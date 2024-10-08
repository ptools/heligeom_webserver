"""
Module containing diverses useful functions.
"""

import re


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
                res_range.errors = e.msg
                valid = False

    return valid
