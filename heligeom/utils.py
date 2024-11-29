"""
Module containing diverses useful functions.
"""

import math
import random
import re

import numpy as np
from ptools import forcefield, measure, superpose, transform
from ptools.heligeom import heli_analyze
from ptools.superpose import rmsd


def parse_resrange(res_range):
    """Parse a string of the form x-y which defines a residue range and return the values.

    To be valid, x < y.

    Parameters
    ----------
    res_range : str
        string containing the residue range

    Returns
    -------
    Tuple:
        Tuple of int of the residue range

    Raises
    ------
    SyntaxError
        If the parsing fails or the value are incorrect.
    """

    res = re.match(r"([0-9]+)-([0-9]+)", res_range)
    if not res:
        raise SyntaxError("Residue range is not in the form of 1-300.")

    min_res = int(res.group(1))
    max_res = int(res.group(2))

    if min_res >= max_res:
        raise SyntaxError("Residue min is superior to residue max")

    return (min_res, max_res)


def check_pair_monomers(chain1, res_range1, chain2, res_range2):
    """Ensure the selection for the pair of monomers is valid.

    A valid selection for one monomer is either a chain attribute
    or/and a residue range.

    Parameters
    ----------
    chain1 : StringField
        Chain value of the 1st monomer
    res_range1 : StringField
        The residue range ("X-Y") for the 1st monomer
    chain2 : StringField
        Chain value of the 2nd monomer
    res_range2 : StringField
        The residue range ("X-Y") for the 2nd monomer

    Returns
    -------
    Bool
        True if both selections are valid. False otherwise
    """
    valid = True
    # Use a zip to test both monomer selections at the same time
    for chain, res_range in zip([chain1, chain2], [res_range1, res_range2], strict=False):
        # Check at least one type of selection is chosen
        if not chain.data and not res_range.data:
            chain.errors = res_range.errors = "No Data specified"
            valid = False
        elif res_range.data:
            # Check if the residue range is consistent
            try:
                parse_resrange(res_range.data)
            except SyntaxError as e:
                res_range.errors = [e.msg]  # Use a list because it's unpacked in HTML.
                valid = False

    return valid


def ener1(rec, lig, cutoff):
    ff = forcefield.AttractForceField1(rec, lig, cutoff)
    return ff.energy()


def adjust(ARb, hpori, enref, Ntarget, Ptarget):
    rec = ARb.copy()
    p = measure.center(rec)

    axe = hpori.unit
    pt0 = hpori.point  # axis point

    ligor1 = rec.copy()
    transform.ab_rotate(ligor1, pt0, pt0 + axe, hpori.angle, degrees=False)
    transform.translate(ligor1, axe * hpori.normtranslation)

    ligor2 = rec.copy()
    transform.ab_rotate(ligor2, pt0, pt0 + axe, -hpori.angle, degrees=False)
    transform.translate(ligor2, -1 * axe * hpori.normtranslation)

    # computes distance from the axis, rayon, radial axis vectn
    # rayori = math.sqrt((pt0.x-p.x)**2 + (pt0.y-p.y)**2 + (pt0.z-p.z)**2) # Hub: unused
    vect = pt0 - p
    vectn = vect - np.dot(vect, axe) * axe
    rayon = np.linalg.norm(vectn)
    vectn = vectn / rayon
    wectn = np.cross(vectn, axe)

    # target geometry  --> angnew is fixed; transnew is fixed
    angnew = 2 * math.pi / Ntarget
    if hpori.angle < 0:
        angnew = -1 * angnew
    transnew = Ptarget / Ntarget

    if hpori.normtranslation < 0:
        transnew = -1 * transnew

    # computes new distance from the axis
    raynew = rayon * math.sin(hpori.angle / 2) / math.sin(angnew / 2)
    pt0new = p + raynew * vectn  # --> axis is laterally translated???
    print("angini,transini", hpori.angle, hpori.normtranslation)
    print("angnew,transnew,rayon", angnew, transnew, raynew)

    ##
    lig = rec.copy()
    transform.ab_rotate(lig, pt0new, pt0new + axe, angnew, degrees=False)
    # transform.ab_rotate(lig, pt0, pt0 + axe, angnew, degrees=False)
    transform.translate(lig, axe * transnew)
    eninit = ener1(rec, lig, 7)

    # debug
    # print "\nAfter modifying ang and trans :\n"
    print(
        "#INIT ",
        Ptarget,
        "/",
        Ntarget,
        f"  Energy ref:  {enref:10.2f}",
        f" Energy init Target:  {eninit:10.2f}",
        f"  -  RMSD ligor ligTarget: {min(rmsd(ligor1, lig), rmsd(ligor2, lig)):7.2f}",
    )
    print(
        f"#                                Fnat rec ligTarget:  {max(measure.fnat(rec, ligor1, rec, lig, 7), measure.fnat(rec, ligor2, rec, lig, 7)):7.2f}"
    )
    # debug

    nbiter = 300
    nbiter2 = 600
    uang = 5.0
    udis = 3.0
    gm = 0.0
    gsda = 2.0
    gsdd = 0.5
    print()
    print("Parameters:   ", nbiter, uang, udis, "\t", nbiter2, gm, gsda, gsdd)
    print()

    # nbiter=500
    # nbaccept=0.
    raytmp = rayon
    enold = eninit

    pt0tmp = pt0.copy()
    dax, day, daz = 0.0, 0.0, 0.0
    # daxbst, daybst, dazbst = 0., 0., 0. # Hub unused
    bestene = 9999
    raybst = raynew
    pt0bst = np.array([0.0, 0.0, 0.0])
    recbst = rec.copy()
    rectmp = rec.copy()
    recnew = rec.copy()
    fnatbst = measure.fnat(rec, ligor1, rec, lig, 7)  # Hub : 7 as cutoff?

    # MC search - constant (N,P) ==> (angnew,transnew)

    ibst = 0
    nbaccept = 0
    nbtry = 0
    for i in range(nbiter):
        if i % 100 == 0:
            print(f"... {i}")
        rectmp = recnew.copy()
        if random.random() > 0.5:  # rotation rec; lig will follow
            dax = math.radians(random.uniform(-1 * uang, uang))  # gaussian rather than uniform???
            day = math.radians(random.uniform(-1 * uang, uang))
            daz = math.radians(random.uniform(-1 * uang, uang))
            transform.ab_rotate(
                rectmp, p, p + vectn, dax
            )  # p = pinit... mais c'est l'axe qui bouge...
            transform.ab_rotate(rectmp, p, p + wectn, day, degrees=False)
            transform.ab_rotate(rectmp, p, p + axe, daz, degrees=False)
        else:  # translation rec; lig will follow
            raytmp = raynew + random.uniform(-1 * udis, udis)  # axis translates along radius
            pt0tmp = p + raytmp * vectn  # rec stays centered on initial p; axis translates

        # apply target transformation to lig
        ligtmp = rectmp.copy()
        transform.ab_rotate(ligtmp, pt0tmp, pt0tmp + axe, angnew, degrees=False)
        transform.translate(ligtmp, axe * transnew)

        # ICI TEST FNAT: si FNAT < 0.25 on passe le tour
        #   attention, rec a bouge, il faut utiliser fnat2
        fn = measure.fnat(rec, rectmp, ligor1, ligtmp, 7)  # Hub : cutoff = 7
        if fn < 0.3:
            print("fnat too low ", i, dax, day, daz, raynew - raytmp, fn)
            continue  # goto next iteration
        nbtry += 1
        enew = ener1(rectmp, ligtmp, 7)

        # MC test
        delta = enew - enold
        if delta < 0:
            delta = 0
            if enew < bestene:
                bestene = enew
                pt0bst = pt0tmp
                recbst = rectmp.copy()
                raybst = raytmp
                fnatbst = fn
                ibst = i
        if random.random() <= math.exp(-delta):  # always true if delta == 0
            enold = enew
            recnew = rectmp
            raynew = raytmp
            nbaccept += 1

    # construct final model
    # ligand
    bestmp = recbst.copy()
    transform.ab_rotate(bestmp, pt0bst, pt0bst + axe, angnew, degrees=False)
    transform.translate(bestmp, axe * transnew)
    # superpose recbst on rec and applymatrix on final and best [CP 12.06.24: superpose the axes instead?]
    ligbst = bestmp.copy()
    matfi = superpose.fit_matrix(rec, recbst)
    transform.move(ligbst, matfi)

    # hpnew = mat_trans_to_screw(superpose(rec,ligbst).matrix)  CORR 4.11.22
    hpbst = heli_analyze(recbst, ligbst)

    # debug
    print()
    print("#ADJ  Nbest raynew raybst-rayinit  bestnrj   rmsd(ligori,ligbst)  fnat")
    print(
        "#ADJ %4d%9.2f%12.2f%12.2f%12.2f%12.2f"
        % (
            ibst,
            raybst,
            raybst - rayon,
            bestene,
            min(rmsd(ligor1, ligbst), rmsd(ligor2, ligbst)),
            fnatbst,
        )
    )
    # debug

    print("acc: ", nbaccept / nbtry, nbtry)

    # new MC run around the best energy
    recnew = recbst
    raynew = raybst
    nbtry2 = 0.0
    nbaccept2 = 0.0
    for i in range(nbiter2):
        if i % 100 == 0:
            print(f"... {i}")
        rectmp = recnew.copy()
        if random.random() > 0.5:  # rotation rec; lig will follow
            dax = math.radians(random.gauss(0.0, gsda))  # gaussian rather than uniform???
            day = math.radians(random.gauss(0.0, gsda))
            daz = math.radians(random.gauss(0.0, gsda))
            transform.ab_rotate(
                rectmp, p, p + vectn, dax
            )  # p = pinit... mais c'est l'axe qui bouge...
            transform.ab_rotate(rectmp, p, p + wectn, day, degrees=False)
            transform.ab_rotate(rectmp, p, p + axe, daz, degrees=False)
        else:  # translation rec; lig will follow
            raytmp = raynew + random.gauss(0.0, gsdd)  # axis translates along radius
            pt0tmp = p + raytmp * vectn  # rec stays centered on initial p; axis translates

        # apply target transformation to lig
        ligtmp = rectmp.copy()
        transform.ab_rotate(ligtmp, pt0tmp, pt0tmp + axe, angnew, degrees=False)
        transform.translate(ligtmp, axe * transnew)

        # ICI TEST FNAT: si FNAT < 0.25 on passe le tour
        #   attention, rec a bouge, il faut utiliser fnat2
        fn = measure.fnat(rec, rectmp, ligor1, ligtmp, 7)  # Hub : cutoff = 7
        if fn < 0.3:
            print("fnat too low ", i, dax, day, daz, raynew - raytmp, fn)
            continue  # goto next iteration
        nbtry2 += 1
        enew = ener1(rectmp, ligtmp, 7)

        # DEGUG
        if enew < -170:
            print("i, enew: ", i, enew)

        # MC test
        delta = enew - enold
        if delta < 0:
            delta = 0
            if enew < bestene:
                bestene = enew
                pt0bst = pt0tmp
                recbst = rectmp.copy()
                raybst = raytmp
                fnatbst = fn
                ibst = i
        if random.random() <= math.exp(-delta):  # always true if delta == 0
            enold = enew
            recnew = rectmp
            raynew = raytmp
            nbaccept2 += 1

    # construct final model
    # ligand
    bestmp = recbst.copy()
    transform.ab_rotate(bestmp, pt0bst, pt0bst + axe, angnew, degrees=False)
    transform.translate(bestmp, axe * transnew)

    # superpose recbst on rec and applymatrix on final and best [CP 12.06.24: superpose the axes instead?]
    ligbst = bestmp.copy()
    matfi = superpose.fit_matrix(rec, recbst)
    transform.move(ligbst, matfi)

    # hpnew = mat_trans_to_screw(superpose(rec,ligbst).matrix)  CORR 4.11.22
    hpbst = heli_analyze(rec, ligbst)

    # debug
    print("#ADJ2  Nbest raynew raybst-rayinit  bestnrj   rmsd(ligori,ligbst)  fnat")
    print(
        "#ADJ2 %3d%10.2f%15.2f%12.2f%15.2f%12.2f"
        % (
            ibst,
            raybst,
            raybst - rayon,
            bestene,
            min(rmsd(ligor1, ligbst), rmsd(ligor2, ligbst)),
            fnatbst,
        )
    )
    # debug

    print("acc2: ", (nbaccept2 * 1.0) / (nbtry2 * 1.0), nbtry2)

    monomers_per_turn = round(360.0 / abs(math.degrees(hpbst.angle)))
    pitch = abs(hpbst.normtranslation * (360.0 / (abs(math.degrees(hpbst.angle)))))

    print(f"hpbst : {hpbst.normtranslation}, {hpbst.angle}")
    print(f"hpbst : mono, pitch: {monomers_per_turn}, {pitch} ")

    return hpbst, ligbst, bestene, min(rmsd(ligor1, ligbst), rmsd(ligor2, ligbst)), fnatbst
