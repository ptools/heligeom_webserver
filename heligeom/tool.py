"""
Ptools module for the calculations
"""

import math

from ptools import RigidBody
from ptools.heligeom import heli_analyze, heli_construct


def move_rigidbody(rb, x=0, y=0, z=0):
    out = rb.copy()
    out.moveby([x, y, z])
    return out

def run(pdb1_file, pdb2_file, n_mer, pdb_out):
    mono1 = RigidBody(pdb1_file)
    mono2 = RigidBody(pdb2_file)

    hp = heli_analyze(mono1, mono2)
    print(hp.point, hp.angle, hp.normtranslation)
    result = heli_construct(mono1, hp, N=n_mer)

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
