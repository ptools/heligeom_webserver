"""
Ptools module for the calculations
"""

from . import utils
import math

from ptools import RigidBody
from ptools.heligeom import heli_analyze, heli_construct

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

    input_structure = RigidBody(pdb_file)

    # Monomer 1:
    if chain_id_M1:
        monomer1 = input_structure.select_chain(chain_id_M1)
        if res_range_M1:
            min_res1, max_res1 = utils.parse_resrange(res_range_M1)
            monomer1 = monomer1.select_residue_range(min_res1, max_res1)
    else:
        min_res1, max_res1 = utils.parse_resrange(res_range_M1)
        monomer1 = input_structure.select_residue_range(min_res1, max_res1)

    # Monomer 2:
    if chain_id_M2:
        monomer2 = input_structure.select_chain(chain_id_M2)
        if res_range_M2:
            min_res2, max_res2 = utils.parse_resrange(res_range_M2)
            monomer2 = monomer2.select_residue_range(min_res2, max_res2)
    else:
        min_res2, max_res2 = utils.parse_resrange(res_range_M2)
        monomer2 = input_structure.select_residue_range(min_res2, max_res2)

    return (monomer1, monomer2)


def screw_parameters(pdb_file, chain_id_M1, chain_id_M2, res_range_M1, res_range_M2):


    monomer1, monomer2 = get_monomers(pdb_file, chain_id_M1, chain_id_M2, res_range_M1, res_range_M2)
    print("monomer1 : ", monomer1.size(), " ; ")
    print(monomer1.writepdb("monomer1.pdb"))
    print("monomer2 : ", monomer2.size(), " ; ")
    print(monomer2.writepdb("monomer2.pdb"))

    hp = heli_analyze(monomer1, monomer2)

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

    return (hp, pitch, monomers_per_turn, direction)


def construct(pdb_file, chain_id_M1, chain_id_M2, res_range_M1, res_range_M2, n_mer, pdb_out):

    monomer1, monomer2 = get_monomers(pdb_file, chain_id_M1, chain_id_M2, res_range_M1, res_range_M2)

    hp = heli_analyze(monomer1, monomer2)
    result = heli_construct(monomer1, hp, N=n_mer)

    result.writepdb(pdb_out)

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

    return (hp, pitch, monomers_per_turn, direction)
