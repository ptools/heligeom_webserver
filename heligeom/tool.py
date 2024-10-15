"""
Backend module for the Heligeom calculations
"""

import math
from dataclasses import dataclass, field

from ptools import RigidBody, io, measure
from ptools.heligeom import chain_intersect, heli_analyze, heli_construct
from ptools.superpose import Screw

from . import utils


@dataclass
class HeligeomMonomer:
    """Define a Heligeom monomer based on the user inputs."""

    input_file: str  #: the path of the input structure file
    chain: str  #: the chain selected if provided.
    residue_range: str  #: the residue index range in the form `X-Y` if provided.

    #: The RigidBody with the correct atoms.
    rb: RigidBody = field(init=False, repr=False)

    def __post_init__(self):
        input_structure = RigidBody.from_pdb(self.input_file)
        self.rb = _create_monomer(input_structure, self.chain, self.residue_range)


def _create_monomer(struct, chain="", res_range=""):
    """Create a monomer from a Rigidbody with the correct extracted atoms
    from either the chains information or the residue range.

    At least one of the 2 optional argument needs to be filled.

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
    """

    monomer = RigidBody()

    if chain:
        monomer = struct.select_chain(chain)
        if res_range:
            min_resid, max_resid = utils.parse_resrange(res_range)
            monomer = monomer.select_residue_range(min_resid, max_resid)
    else:
        min_resid, max_resid = utils.parse_resrange(res_range)
        monomer = struct.select_residue_range(min_resid, max_resid)

    if monomer.size() == 0:
        raise TypeError("Monomer has a size of 0.")

    return monomer


class HeligeomInterface:
    """Represents a Heligeom interface between 2 defined monomers."""

    #: The 1st monomer of the interface
    monomer1: HeligeomMonomer
    #: The 2nd monomer of the interface
    monomer2: HeligeomMonomer
    #: The Screw transformation defining the interface
    hp: Screw

    def __init__(self, pdb_file, chain_id_M1, chain_id_M2, res_range_M1, res_range_M2):
        self.monomer1 = HeligeomMonomer(pdb_file, chain_id_M1, res_range_M1)
        self.monomer2 = HeligeomMonomer(pdb_file, chain_id_M2, res_range_M2)

        self.hp = Screw()

    def compute_screw(self, force=False):
        """Compute the screw transformation between the 2 monomers.

        Parameters
        ----------
        force : bool, optional
            Force recompute the transformation, by default False
        """

        # Retrieve the offset between the first resid of the 2 monomers
        try:
            min_res1, _ = utils.parse_resrange(self.monomer1.residue_range)
            min_res2, _ = utils.parse_resrange(self.monomer2.residue_range)
            delta_resid = min_res2 - min_res1
            # Handle missing residues with an intersection
            rb1, rb2 = chain_intersect(self.monomer1.rb, self.monomer2.rb, delta_resid=delta_resid)
        except SyntaxError:
            rb1 = self.monomer1.rb
            rb2 = self.monomer2.rb

        if rb1.size() == 0:
            raise ValueError("Monomer 1 has a size of 0.")
        if rb2.size() == 0:
            raise ValueError("Monomer 2 has a size of 0.")

        # Use CA for computating parameters
        monomer1_CA = rb1.select_atom_type("CA")
        monomer2_CA = rb2.select_atom_type("CA")

        self.hp = heli_analyze(monomer1_CA, monomer2_CA)

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

        return (pitch, monomers_per_turn, direction, dmin, dmax)

    def construct_oligomer(self, ncopies, z_align, fileout):
        """Based on the interface, construct an oligomer with
        `ncopies` of the 1st monomer.

        This oligomer will be written as a PDB in `fileout`.

        Parameters
        ----------
        ncopies : int
            number of copies of the 1st monomer to be cosntructed.
        z_align: Bool
            Wether the oligomer should be aligned on the Z axis.
        fileout : str
            Path to output file.
        """

        result = heli_construct(self.monomer1.rb, self.hp, N=ncopies, Z=z_align)

        io.write_pdb(result, fileout)  # type: ignore

    @classmethod
    def compute_fnat(cls, heli_interface1, heli_interface2, cutoff=5):
        """Compute the FNAT (fraction of native contacts) between 2 interfaces.

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
        return measure.fnat(
            heli_interface1.monomer1.rb,
            heli_interface1.monomer2.rb,
            heli_interface2.monomer1.rb,
            heli_interface2.monomer2.rb,
            cutoff,
        )
