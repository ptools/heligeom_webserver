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

    res = re.match(r'([0-9]+)-([0-9]+)', res_range)
    if not res:
        raise SyntaxError("Residue range is not in the form of 1-300.")

    min_res = int(res.group(1))
    max_res = int(res.group(2))

    if min_res >= max_res:
        raise SyntaxError("Residue min is superior to residue max")

    return (min_res, max_res)
