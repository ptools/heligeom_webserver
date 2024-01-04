"""
Ptools module for the calculations
"""

from . import utils
import math

from ptools import RigidBody
from ptools import measure
from ptools.heligeom import heli_analyze, heli_construct, chain_intersect
from ptools import io

def get_monomers(pdb_file, chain_id_M1, chain_id_M2, res_range_M1, res_range_M2):
    """Return the 2 monomers exracted from the `pdb_file` with the correct extracted atoms
    from either the chains information or the residue range.

    Parameters
    ----------
    pdb_file : str
        name of the PDB file
    chain_id_M1 : str
        Chain of the 1st monomer
    chain_id_M2 : str
        Chain of the 2nd monomer
    res_range_M1 : str
        string of the form 'X-Y' where X and Y are the residue ID minimal and maximal
    res_range_M2 : str
        string of the form 'X-Y' where X and Y are the residue ID minimal and maximal

    Returns
    -------
    Tuple of RigidBody
        the 2 monomers
    """

    input_structure = RigidBody.from_pdb(pdb_file)

    # Monomer 1:
    if chain_id_M1:
        monomer1 = input_structure.select_chain(chain_id_M1)
        if res_range_M1:
            min_res1, max_res1 = utils.parse_resrange(res_range_M1)
            monomer1 = monomer1.select_residue_range(min_res1, max_res1)
    else:
        min_res1, max_res1 = utils.parse_resrange(res_range_M1)
        monomer1 = input_structure.select_residue_range(min_res1, max_res1)

    if monomer1.size() == 0:
        raise ValueError("Monomer 1 has a size of 0.")

    # Monomer 2:
    if chain_id_M2:
        monomer2 = input_structure.select_chain(chain_id_M2)
        if res_range_M2:
            min_res2, max_res2 = utils.parse_resrange(res_range_M2)
            monomer2 = monomer2.select_residue_range(min_res2, max_res2)
    else:
        min_res2, max_res2 = utils.parse_resrange(res_range_M2)
        monomer2 = input_structure.select_residue_range(min_res2, max_res2)

    if monomer2.size() == 0:
        raise ValueError("Monomer 2 has a size of 0.")

    return (monomer1, monomer2)


def screw_parameters(pdb_file, chain_id_M1, chain_id_M2, res_range_M1, res_range_M2):

    # Create monomers from input data
    monomer1, monomer2 = get_monomers(pdb_file, chain_id_M1, chain_id_M2, res_range_M1, res_range_M2)

    # Retrieve the offset between the first resid of the 2 monomers
    try:
        min_res1, _ = utils.parse_resrange(res_range_M1)
        min_res2, _ = utils.parse_resrange(res_range_M2)
        delta_resid = min_res2 - min_res1
        # Handle missing residues with an intersection
        rb1, rb2 = chain_intersect(monomer1, monomer2, delta_resid=delta_resid)
    except SyntaxError:
        rb1 = monomer1
        rb2 = monomer2

    if rb1.size() == 0:
        raise ValueError("Monomer 1 has a size of 0.")
    if rb2.size() == 0:
        raise ValueError("Monomer 2 has a size of 0.")

    # Use CA for computating parameters
    monomer1_CA = rb1.select_atom_type("CA")
    monomer2_CA = rb2.select_atom_type("CA")

    hp = heli_analyze(monomer1_CA, monomer2_CA)

    #Retrieve N & Pitch
    rotation_angle_degrees = math.degrees(hp.angle)
    if abs(hp.angle) > 0.0:
        pitch = abs(hp.normtranslation * (360. / abs(rotation_angle_degrees)))
        monomers_per_turn = 360. / abs(rotation_angle_degrees)
    else:
        # Attention: zeros here indicate NaN (divide by zero)
        pitch = 0.0
        monomers_per_turn = 0.0
    direction = "right-handed" if hp.angle * hp.normtranslation > 0 else "left-handed"

    dmin, dmax = measure.minmax_distance_to_axis(monomer1, hp.unit, center=hp.point)

    return (hp, pitch, monomers_per_turn, direction, dmin, dmax)


def construct(pdb_file, chain_id_M1, chain_id_M2, res_range_M1, res_range_M2, n_mer, z_align, pdb_out):

    # Create monomers from input data
    monomer1, monomer2 = get_monomers(pdb_file, chain_id_M1, chain_id_M2, res_range_M1, res_range_M2)

    # Retrieve the offset between the first resid of the 2 monomers
    try:
        min_res1, _ = utils.parse_resrange(res_range_M1)
        min_res2, _ = utils.parse_resrange(res_range_M2)
        delta_resid = min_res2 - min_res1
        # Handle missing residues with an intersection
        rb1, rb2 = chain_intersect(monomer1, monomer2, delta_resid=delta_resid)
    except SyntaxError:
        rb1 = monomer1
        rb2 = monomer2

    if rb1.size() == 0:
        raise ValueError("Monomer 1 has a size of 0.")
    if rb2.size() == 0:
        raise ValueError("Monomer 2 has a size of 0.")


    # Use CA for computating parameters
    monomer1_CA = rb1.select_atom_type("CA")
    monomer2_CA = rb2.select_atom_type("CA")

    hp = heli_analyze(monomer1_CA, monomer2_CA)
    result = heli_construct(monomer1, hp, N=n_mer,Z=z_align)

    io.write_pdb(result, pdb_out)

    #Retrieve N & Pitch
    rotation_angle_degrees = math.degrees(hp.angle)
    if abs(hp.angle) > 0.0:
        pitch = abs(hp.normtranslation * (360. / abs(rotation_angle_degrees)))
        monomers_per_turn = 360. / abs(rotation_angle_degrees)
    else:
        # Attention: zeros here indicate NaN (divide by zero)
        pitch = 0.0
        monomers_per_turn = 0.0
    direction = "right-handed" if hp.angle * hp.normtranslation > 0 else "left-handed"

    dmin, dmax = measure.minmax_distance_to_axis(monomer1, hp.unit, center=hp.point)

    return (hp, pitch, monomers_per_turn, direction, dmin, dmax)


def analyze_fnat(pdb_file, chain_id_M1, chain_id_M2, res_range_M1, res_range_M2,
                pdb_file2, chain2_id_M1, chain2_id_M2, res_range2_M1, res_range2_M2):


    monomer1, monomer2 = get_monomers(pdb_file, chain_id_M1, chain_id_M2, res_range_M1, res_range_M2)
    monomer1bis, monomer2bis = get_monomers(pdb_file2, chain2_id_M1, chain2_id_M2, res_range2_M1, res_range2_M2)

    return measure.fnat(monomer1, monomer2, monomer1bis, monomer2bis, cutoff=5)
