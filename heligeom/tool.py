"""
Backend module for the Heligeom calculations
"""

import math
import pathlib
import re
from dataclasses import dataclass, field

import numpy as np
import ptools
from ptools import RigidBody, io, measure, reduce
from ptools import superpose
from ptools.heligeom import chain_intersect, heli_analyze, heli_construct
from ptools.pairlist import PairList
from ptools.superpose import Screw, rmsd

from . import utils
from .forms import InputStructures


# = Custom selections functions are faster than semantic ones i.e rb.select("chain X") ========
def _custom_select_chains(rb, chain):
    indices = np.where(rb.atom_properties["chains"].values == chain)[0]
    return rb[indices]


def _custom_select_residue_range(rb, start, end):
    indices = np.where(
        np.logical_and(
            rb.atom_properties["residue_indices"].values >= start,
            rb.atom_properties["residue_indices"].values <= end,
        )
    )[0]
    return rb[indices]


# = Custom Exceptions ========
class MonomerSizeZeroError(BaseException):
    """Raised when a RigidBody has a size of 0 atoms."""

    pass


class MonomersDifferentSizeError(BaseException):
    """Raised when 2 RigidBodys have differnt sizes."""

    pass


# =================================


@dataclass
class HeligeomMonomer:
    """Define a Heligeom monomer based on the user inputs."""

    input_data: str | RigidBody  #: the input data. Can be a file or a RigidBody structure
    chain: str  #: the chain selected if provided.
    residue_range: str  #: the residue index range in the form `X-Y` if provided.

    #: The RigidBody with the correct atoms.
    rb: RigidBody = field(init=False, repr=False)
    #: the string of the molstar selection of the monomer
    molstar_selection: str = field(init=False, repr=False)

    def __post_init__(self):
        if isinstance(self.input_data, str):
            input_structure = RigidBody.from_pdb(self.input_data)
        else:
            input_structure = self.input_data
        (self.rb, self.molstar_selection) = _create_monomer(
            input_structure, self.chain, self.residue_range
        )

    def size(self):
        """Return the number of atoms of the Rigidbody"""
        return self.rb.size()


def _create_monomer(struct, chain="", res_range=""):
    """Create a monomer from a Rigidbody with the correct extracted atoms
    from either the chains information or the residue range.

    At least one of the 2 optional argument needs to be filled.

    All heteros atoms are discarded.

    It also construct the molstar selection of the monomer, like
    "struct_asym_id: 'B', start_residue_number:1, end_residue_number:100"
    See https://github.com/molstar/pdbe-molstar/wiki/3.-Helper-Methods.

    Parameters
    ----------
    struct : RigidBody
        the global structure
    chain : str, optional
        chain of the monomer.
    res_range : str, optional
        string of the form 'X-Y' where X and Y are the residue ID minimal and maximal.

    Raises
    ------
    TypeError
        if the monomer has a size of 0.

    Returns
    -------
    RigidBody
        the monomer constructed
    string
        the molstar selection of the monomer
    """

    molstar_selection = ""

    if chain:
        monomer = _custom_select_chains(struct, chain)
        molstar_selection = f"struct_asym_id: '{ chain }'"
        if res_range:
            min_resid, max_resid = utils.parse_resrange(res_range)
            monomer = _custom_select_residue_range(monomer, min_resid, max_resid)
            molstar_selection += (
                f", start_residue_number: {min_resid}, end_residue_number: {max_resid}"
            )
    else:
        min_resid, max_resid = utils.parse_resrange(res_range)
        monomer = _custom_select_residue_range(struct, min_resid, max_resid)
        molstar_selection = f"start_residue_number: {min_resid}, end_residue_number: {max_resid}"

    return monomer.copy(), molstar_selection


class HeligeomInterface:
    """Represents a Heligeom interface between 2 defined monomers."""

    #: The 1st monomer of the interface
    monomer1: HeligeomMonomer
    #: The 2nd monomer of the interface
    monomer2: HeligeomMonomer
    #: The RigidBody of the 1st monomer restricted to the core region
    core_monomer1: RigidBody
    #: The RigidBody of the 2nd monomer restricted to the core region
    core_monomer2: RigidBody
    #: The Screw transformation defining the interface
    hp: Screw
    #: The oligomer construction made of monomer1
    oligomer: RigidBody = field(init=False, repr=False)
    #: colors of the 1st monomer ([light, strong])
    colors_monomer1: list = ["#DCFCFE", "#02AAB3"]
    #: colors of the 2nd monomer ([light, strong])
    colors_monomer2: list = ["#FFD1A2", "#F57C00"]

    #: molstar selection of the core monomer 1
    molstar_select_core_monomer1: str
    #: molstar selection of the core monomer 2
    molstar_select_core_monomer2: str

    def __init__(
        self,
        pdb_file,
        chain_id_M1,
        chain_id_M2,
        res_range_M1,
        res_range_M2,
        core_filter,
        core_region1,
        core_region2,
    ):
        # Parse the input file
        global_rb = RigidBody.from_pdb(pdb_file)

        # Add manually the occupency property if it's not read
        if "occupancies" not in global_rb.atom_properties:
            global_rb.add_atom_property("occupancy", "occupancies", [0] * len(global_rb))
        else:
            global_rb.occupancies = [0] * len(global_rb)

        # Discard hetero and water atoms
        protein = global_rb.select("not hetero and not water").copy()
        # Renumber atoms of the protein for the molstar selections
        protein.reset_atom_indices(protein[0].index)

        # Create HeligeomMonomers with the selecting region
        self.monomer1 = HeligeomMonomer(protein, chain_id_M1, res_range_M1)
        self.monomer2 = HeligeomMonomer(protein, chain_id_M2, res_range_M2)

        if self.monomer1.size() == 0:
            raise MonomerSizeZeroError(
                "Monomer 1 defined with chain and/or residues has a size of 0 atoms."
            )
        if self.monomer2.size() == 0:
            raise MonomerSizeZeroError(
                "Monomer 2 defined with chain and/or residues has a size of 0 atoms."
            )

        # Handle missing residues
        # Retrieve the offset between the first resid of the 2 monomers
        try:
            min_res1, _ = utils.parse_resrange(self.monomer1.residue_range)
            min_res2, _ = utils.parse_resrange(self.monomer2.residue_range)
            delta_resid = min_res2 - min_res1
            # Handle missing residues with an intersection
            self.monomer1.rb, self.monomer2.rb = chain_intersect(
                self.monomer1.rb, self.monomer2.rb, delta_resid=delta_resid
            )
        except SyntaxError:
            # No Residue range provided so assume the delta resid = 0
            self.monomer1.rb, self.monomer2.rb = chain_intersect(
                self.monomer1.rb, self.monomer2.rb, delta_resid=0
            )

        if self.monomer1.size() == 0:
            raise MonomerSizeZeroError(
                "Monomer 1 defined with chain and/or residues has a size of 0 atoms."
            )
        if self.monomer2.size() == 0:
            raise MonomerSizeZeroError(
                "Monomer 2 defined with chain and/or residues has a size of 0 atoms."
            )

        # Core region defined?
        if core_filter == "manual":
            self.core_monomer1, self.molstar_select_core_monomer1 = create_core_monomer(  # type: ignore
                self.monomer1.rb, core_region1
            )
            self.core_monomer2, self.molstar_select_core_monomer2 = create_core_monomer(  # type: ignore
                self.monomer2.rb, core_region2
            )

            if self.monomer1.chain:
                self.molstar_select_core_monomer1 = f"struct_asym_id: '{ self.monomer1.chain }', {self.molstar_select_core_monomer1}"  # noqa: E501
            if self.monomer2.chain:
                self.molstar_select_core_monomer2 = f"struct_asym_id: '{ self.monomer2.chain }', {self.molstar_select_core_monomer2}"  # noqa: E501
        else:
            # The core region is the whole monomer
            self.monomer1.rb.occupancies = [1] * len(self.monomer1.rb)
            self.core_monomer1 = self.monomer1.rb
            self.core_monomer2 = self.monomer2.rb
            self.molstar_select_core_monomer1 = self.monomer1.molstar_selection
            self.molstar_select_core_monomer2 = self.monomer2.molstar_selection

        if self.core_monomer1.size() == 0:
            raise MonomerSizeZeroError("Monomer 1 defined with core regions has a size of 0 atoms.")
        if self.core_monomer2.size() == 0:
            raise MonomerSizeZeroError("Monomer 2 defined with core regions has a size of 0 atoms.")

        monomer1_CA = self.core_monomer1.select("name CA")
        monomer2_CA = self.core_monomer2.select("name CA")
        if monomer1_CA.size() != monomer2_CA.size():
            raise MonomersDifferentSizeError(
                f"Monomer 1 & 2 have different number of Carbon-α ({monomer1_CA.size()} vs {monomer2_CA.size()} )."
            )

        self.hp = Screw()

    def compute_screw(self):
        """Compute the screw transformation between the core 2 monomers."""

        # Use CA to compute Screw parameters
        monomer1_CA = self.core_monomer1.select("name CA")
        monomer2_CA = self.core_monomer2.select("name CA")

        self.hp = heli_analyze(monomer1_CA, monomer2_CA)
        rmsd_value = rmsd(monomer1_CA, monomer2_CA, do_fit=True)

        # Retrieve N & Pitch
        rotation_angle_degrees = math.degrees(self.hp.angle)
        if abs(self.hp.angle) > 0.0:
            pitch = abs(self.hp.normtranslation * (360.0 / abs(rotation_angle_degrees)))
            monomers_per_turn = 360.0 / abs(rotation_angle_degrees)
        else:
            # Attention: zeros here indicate NaN (divide by zero)
            pitch = 0.0
            monomers_per_turn = 0.0
        direction = "right-handed" if self.hp.angle * self.hp.normtranslation > 0 else "left-handed"

        dmin, dmax = measure.minmax_distance_to_axis(
            self.monomer1.rb, self.hp.unit, center=self.hp.point
        )

        return (pitch, monomers_per_turn, direction, dmin, dmax, rmsd_value)

    def construct_oligomer(self, fileout, ncopies, z_align=True):
        """Based on the interface, construct an oligomer with
        `ncopies` of the 1st monomer.

        This oligomer will be written as a PDB in `fileout`.

        Parameters
        ----------
        fileout : str
            Path to output file.
        ncopies : int
            number of copies of the 1st monomer to be cosntructed.
        z_align: Bool
            Wether the oligomer should be aligned on the Z axis.
        """

        self.oligomer = heli_construct(self.monomer1.rb, self.hp, N=ncopies, Z=z_align)

        io.write_pdb(self.oligomer, fileout)  # type: ignore

    def heli_ring(self, fileout, ncopies, z_align=True):
        # Create an AttractRigidBody from a RigidBody
        # First need to reduce the structures to Coarse grained structures
        rb1_beads = reduce.reducer.reduce(self.monomer1.rb)
        arec = ptools.AttractRigidBody(rb1_beads)  # type: ignore

        rb2_beads = reduce.reducer.reduce(self.monomer2.rb)
        alig = ptools.AttractRigidBody(rb2_beads)  # type: ignore

        # Compute reference energy
        eref = utils.ener1(arec, alig, 7)

        monomers_per_turn = round(360.0 / abs(math.degrees(self.hp.angle)))

        # Call adjust with Ptarget = 0 to create a ring
        newhp, newARb, bestener, rms, fnat = utils.adjust(arec, self.hp, eref, monomers_per_turn, 0)

        self.oligomer = heli_construct(self.monomer1.rb, newhp, N=ncopies, Z=z_align)
        io.write_pdb(self.oligomer, fileout)  # type: ignore

        return (newhp, rms, fnat)

    def save_monomers(self, fileout):
        """Save both monomers to a PDB file in `fileout`.

        Parameters
        ----------
        fileout : pathlib.Path
            Path to output file.
        """
        concat = self.monomer1.rb + self.monomer2.rb
        io.write_pdb(concat, fileout)  # type: ignore

    def interface_atoms_oligomer(self, cutoff=5):
        """Returns the indexes of the residue atoms of monomer 1 and monomer 1' in contacts from
        the oligomer structure computed.

        In the oligomer structure, monomer 1 will be the chain A and monomer 1' will be the chain B.

        Parameters
        ----------
        cutoff : float, optional
            distance in Angstrom defining a contact between 2 residues, by default 5

        Returns
        -------
        tuple[np.ndarray, np.ndarray]:
            monomer 1 and monomer 1' atom indexes.
        """

        mono1 = self.oligomer.select("chain A")
        mono1prime = self.oligomer.select("chain B")

        mono1_atom_indexes = []
        mono2_atom_indexes = []

        lhs_atom_ids, rhs_atom_ids = PairList(mono1, mono1prime, cutoff).raw_contacts()
        lhs_residue_ids = mono1.residue_indices[lhs_atom_ids]
        rhs_residue_ids = mono1prime.residue_indices[rhs_atom_ids]

        if lhs_atom_ids.size != 0:
            mono1_atom_indexes = mono1.select_by_property("resid", lhs_residue_ids).indices
        if rhs_atom_ids.size != 0:
            mono2_atom_indexes = mono1prime.select_by_property("resid", rhs_residue_ids).indices

        return mono1_atom_indexes, mono2_atom_indexes

    def molstar_selection_interface_oligomer(self):
        """Returns the molstar selection as a string of the atom indexes belonging
        to the interface of the oligomer.

        See `interface_atoms_oligomer()`.
        In the oligomer structure, monomer 1 will be the chain A and monomer 1' will be the chain B.

        Returns
        -------
        str
            the molstar selection as a string of the atom indexes
        """
        mono1_atom_indexes, mono2_atom_indexes = self.interface_atoms_oligomer()

        selection = (
            # selection of the monomer 1
            f"{{ struct_asym_id: 'A', color:'{ self.colors_monomer1[0] }' }},"
            # selection of the monomer 1 atoms at the interface
            f"{{ struct_asym_id: 'A', atom_id: [{", ".join([str(i) for i in mono1_atom_indexes])}], representation:'ball-and-stick', representationColor:'{self.colors_monomer1[1]}', color:'{self.colors_monomer1[0]}', focus:true }},"
            # selection of the monomer 2
            f"{{ struct_asym_id: 'B', color:'{ self.colors_monomer2[0] }' }},"
            # selection of the monomer 2 atoms at the interface
            f"{{ struct_asym_id: 'B', atom_id: [{", ".join([str(i) for i in mono2_atom_indexes])}], representation:'ball-and-stick', representationColor:'{self.colors_monomer2[1]}', color:'{self.colors_monomer2[0]}', focus:true }},"
        )

        return selection

    def save_csv_atom_contacts(self, pl, filout):
        """Save to a file `filout` all contacts stored in `pl`
        with atom information.

        The format will be .csv with those columns:
            Index, Atom1 Id, Atom1 Name, Residue1 Id, Residue1 Name, Atom2 Id, Atom2 Name, Residue2 Id, Residue2 Name, Distance.

        Parameters
        ----------
        pl : Pairlist
            the contacts computed
        filout : str
            the filename to write the contacts results.
        """

        # Retrieve the list of atoms for each partner from the Pairlist
        lhs_atom_ids, rhs_atom_ids = pl.raw_contacts()

        # start at 1
        indexes = np.arange(1, len(lhs_atom_ids) + 1)

        # Retrieve info from the first monomer
        lhs_atom_names = self.core_monomer1[lhs_atom_ids].names
        lhs_atom_resids = self.core_monomer1[lhs_atom_ids].residue_indices
        lhs_atom_resnames = self.core_monomer1[lhs_atom_ids].residue_names

        # Retrieve info from the 2nd monomer
        rhs_atom_names = self.core_monomer2[rhs_atom_ids].names
        rhs_atom_resids = self.core_monomer2[rhs_atom_ids].residue_indices
        rhs_atom_resnames = self.core_monomer2[rhs_atom_ids].residue_names

        # :( cast the  contacts distance in string
        distances = [f"{dist:.2f}" for dist in pl.distances()]

        np.savetxt(
            filout,
            [
                i
                for i in zip(
                    indexes,
                    lhs_atom_ids,
                    lhs_atom_names,
                    lhs_atom_resids,
                    lhs_atom_resnames,
                    rhs_atom_ids,
                    rhs_atom_names,
                    rhs_atom_resids,
                    rhs_atom_resnames,
                    distances,
                    strict=True,
                )
            ],
            delimiter=",",
            fmt="%s",
            header="Index, Atom1 Id, Atom1 Name, Residue1 Id, Residue1 Name, Atom2 Id, Atom2 Name, Residue2 Id, Residue2 Name, Distance",
            comments="",
        )

    def interface_atoms(self, cutoff=5, fileout=None):
        """Returns the indexes of the residue atoms of monomer 1 and monomer 2 in contacts

        Parameters
        ----------
        cutoff : float, optional
            distance in Angstrom defining a contact between 2 residues, by default 5
        filout: str, optional
            If provided, the list of contacts will be save in .csv format in this filename.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]:
            monomer 1 and monomer 2 atom indexes.
        """

        pl = PairList(self.core_monomer1, self.core_monomer2, cutoff)

        if fileout and isinstance(fileout, str | pathlib.Path):
            self.save_csv_atom_contacts(pl, fileout)

        lhs_atom_ids, rhs_atom_ids = pl.raw_contacts()

        lhs_residue_ids = self.core_monomer1.residue_indices[lhs_atom_ids]
        rhs_residue_ids = self.core_monomer2.residue_indices[rhs_atom_ids]

        mono1_atom_indexes = []
        mono2_atom_indexes = []

        if lhs_atom_ids.size != 0:
            mono1_atom_indexes = self.core_monomer1.select_by_property(
                "resid", set(lhs_residue_ids)
            ).indices
        if rhs_atom_ids.size != 0:
            mono2_atom_indexes = self.core_monomer2.select_by_property(
                "resid", set(rhs_residue_ids)
            ).indices

        return mono1_atom_indexes, mono2_atom_indexes

    def molstar_selection_monomers(self, fileout=None):
        """Returns the molstar selection as a string of the 2 monomers
        and the atom indexes belonging to the interface of monomer 1 and 2.

        See `interface_atoms()` for the contacts.

        Returns
        -------
        str
            the molstar selection as a string
        """
        mono1_atom_indexes, mono2_atom_indexes = self.interface_atoms(cutoff=5, fileout=fileout)

        selection = (
            # selection of the monomer 1
            f"{{ { self.monomer1.molstar_selection }, color:'{ self.colors_monomer1[0] }' }},"
            # selection of the monomer 1 atoms at the interface
            f"{{ { self.monomer1.molstar_selection }, atom_id: [{", ".join([str(i) for i in mono1_atom_indexes])}], representation:'ball-and-stick', representationColor:'{self.colors_monomer1[1]}', color:'{self.colors_monomer1[0]}', focus:true }},"
            # selection of the monomer 2
            f"{{ { self.monomer2.molstar_selection }, color:'{ self.colors_monomer2[0] }' }},"
            # selection of the monomer 2 atoms at the interface
            f"{{ { self.monomer2.molstar_selection }, atom_id: [{", ".join([str(i) for i in mono2_atom_indexes])}], representation:'ball-and-stick', representationColor:'{self.colors_monomer2[1]}', color:'{self.colors_monomer2[0]}', focus:true }},"
        )

        return selection

    @classmethod
    def compute_fnat(cls, heli_interface1, heli_interface2, cutoff=5):
        """Compute the FNAT (fraction of native contacts) between 2 interfaces.

        Handle situation where there is an offset in the residue indices between the 2 interfaces.

        Parameters
        ----------
        heli_interface1 : HeligeomInterface
            the 1st heligeom interface
        heli_interface2 : HeligeomInterface
            the 2nd heligeom interface
        cutoff : float, optional
            distance in Angstrom defining a contact between 2 residues, by default 5

        Returns
        -------
        float
            the FNAT value
        """
        rec1 = heli_interface1.monomer1.rb
        lig1 = heli_interface1.monomer2.rb
        rec2 = heli_interface2.monomer1.rb
        lig2 = heli_interface2.monomer2.rb
        res_pair1 = measure.contacts_by_residue(rec1, lig1, cutoff)
        if len(res_pair1) == 0:
            return 0

        res_pair2 = measure.contacts_by_residue(rec2, lig2, cutoff)

        # Compute the offset in the 1st interface, assuming the offset in the 2nd interface is the same
        offset_resid = (
            heli_interface1.monomer2.rb[0].residue_index
            - heli_interface1.monomer1.rb[0].residue_index
        )
        # Inverse the offset if the 1st monomer residue is smaller than the 2nd one
        if offset_resid < 0:
            offset_resid = -offset_resid

        # Apply the offset to the 2nd interface
        new_res_pair2 = set(map(lambda t: (t[0] - offset_resid, t[1] - offset_resid), res_pair2))

        intersect = res_pair1 & new_res_pair2
        return len(intersect) / len(res_pair1)


def create_core_monomer(rb, core_region):
    """Based on the string `core_region`, construct a sub Rigidbody of `rb`
    containing only the residus selected.

    The string should be in the form of "X-Y, Z-A,..." with X,Y,Z,A as residue numbers.


    Parameters
    ----------
    rb : _type_
        _description_
    core_region : String
        litteral string defining a core region, in the form of "X-Y, Z-A,..."

    Returns
    -------
    Rigidbody
        a new RigidBody containing only the core region.
        If no core region found, return the input `rb`.
    """
    # Parse core region input
    res = re.match(InputStructures.cls_regexp_core, core_region)
    if not res:
        return rb

    molstar_selection = ""
    core_rb = RigidBody()
    for res_range in core_region.split(","):
        stripped_res_range = res_range.strip()
        min_res, max_res = utils.parse_resrange(stripped_res_range)
        molstar_selection += f"start_residue_number: {min_res}, end_residue_number: {max_res}"
        tmp_core_rb = rb.select(f"resid {min_res}:{max_res}")
        # Core region = 1 in occupancy
        tmp_core_rb.occupancies = [1] * len(tmp_core_rb)
        core_rb += tmp_core_rb

    return core_rb, molstar_selection
