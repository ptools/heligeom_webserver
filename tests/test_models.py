"""
Unit tests for the Heligeom models
"""

from heligeom.models import UserInputs


def test_new_data():
    """
    Test the insertion of user input data
    """

    data = UserInputs("32300a0108a84a26afb3f23eb6bafba2", "1kx2.pdb", "A", "B", "", "40:300")
    assert data.request_id == "32300a0108a84a26afb3f23eb6bafba2"
    assert data.pdb_filename == "1kx2.pdb"
    assert data.chain1_id == "A"
    assert data.chain2_id == "B"
    assert data.res_range1 == ""
    assert data.res_range2 == "40:300"
