"""
Ptools module for the calculations
"""

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
    res_range_M1 : [type]
        [description]
    res_range_M2 : [type]
        [description]

    Returns
    -------
    Tuple of RigidBody
        the 2 monomers
    """

    input_structure = RigidBody(pdb_file)

    if chain_id_M1 != "" and chain_id_M2 != "":
        monomer1 = input_structure.select_chain(chain_id_M1)
        monomer2 = input_structure.select_chain(chain_id_M2)
    else:
        monomer1 = input_structure.select_residue_range(res_range_M1)
        monomer2 = input_structure.select_residue_range(res_range_M2)


    return (monomer1, monomer2)


def screw_parameters(pdb_file, chain_id_M1, chain_id_M2, res_range_M1, res_range_M2):


    monomer1, monomer2 = get_monomers(pdb_file, chain_id_M1, chain_id_M2, res_range_M1, res_range_M2)

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
