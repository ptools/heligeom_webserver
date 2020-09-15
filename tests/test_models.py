"""
Unit tests for the Heligeom models
"""

from heligeom.models import User_Inputs


def test_new_data():
    """
    Test the insertion of user input data
    """

    data = User_Inputs("32300a0108a84a26afb3f23eb6bafba2","1kx2.pdb","1bta.pdb", 5)
    assert data.request_id == "32300a0108a84a26afb3f23eb6bafba2"
    assert data.pdb1_filename =="1kx2.pdb"
    assert data.pdb2_filename == "1bta.pdb"
    assert data.n_mer == 5
    assert not data.z_align
