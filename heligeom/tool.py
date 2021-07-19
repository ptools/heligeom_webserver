"""
Ptools module for the calculations
"""

import math

from ptools import RigidBody
from ptools.heligeom import heli_analyze, heli_construct


def run(pdb_file, chain_id_M1, chain_id_M2, res_range_M1, res_range_M2, n_mer, pdb_out):
    input_structure = RigidBody(pdb_file)

    if chain_id_M1 != "" and chain_id_M2 != "":
        monomere1 = input_structure.select_chain(chain_id_M1)
        monomere2 = input_structure.select_chain(chain_id_M2)
    else:
        monomere1 = input_structure.select_residue_range(res_range_M1)
        monomere2 = input_structure.select_residue_range(res_range_M2)

    hp = heli_analyze(monomere1, monomere2)
    print(hp.point, hp.angle, hp.normtranslation)
    result = heli_construct(monomere1, hp, N=n_mer)

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
