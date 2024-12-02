from . import db


class UserInputs(db.Model):
    """Class defining the data stored in the Database.

    Matches the fields of the `InputStructures` class.
    """

    # Uniq ID defining the path of the results page
    request_id = db.Column(db.String(32), unique=True, primary_key=True)
    # PDB file obtained from the user input
    pdb_filename = db.Column(db.String(100), nullable=False)
    # Chain ID to define the 1st monomer
    chain1_id = db.Column(db.String(2), default="")
    # Chain ID to define the 2nd monomer
    chain2_id = db.Column(db.String(2), default="")
    # Residue range to define the 1st monomer
    res_range1 = db.Column(db.String(50), default="")
    # Residue range to define the 2nd monomer
    res_range2 = db.Column(db.String(50), default="")
    # Type of core region (monomer 1)
    core_filter1 = db.Column(db.String(10), default="")
    # List of residue range to define the core region (monomer 1)
    core_region1 = db.Column(db.String(50), default="")
    # Type of core region (monomer 2)
    core_filter2 = db.Column(db.String(10), default="")
    # List of residue range to define the core region (monomer 2)
    core_region2 = db.Column(db.String(50), default="")

    # For a second oligomer
    second_oligomer = db.Column(db.Boolean, default=False)
    pdb_filename_2nd = db.Column(db.String(100), default="")
    chain1bis_id = db.Column(db.String(2), default="")
    chain2bis_id = db.Column(db.String(2), default="")
    res_range1bis = db.Column(db.String(50), default="")
    res_range2bis = db.Column(db.String(50), default="")
    core_filter1bis = db.Column(db.String(10), default="")
    core_region1bis = db.Column(db.String(50), default="")
    core_filter2bis = db.Column(db.String(10), default="")
    core_region2bis = db.Column(db.String(50), default="")

    def __init__(
        self,
        request_id,
        pdb_filename="",
        chain1_id="",
        chain2_id="",
        res_range1="",
        res_range2="",
        core_filter1="all",
        core_region1="",
        core_filter2="all",
        core_region2="",
    ):
        self.request_id = request_id
        self.pdb_filename = pdb_filename
        self.chain1_id = chain1_id
        self.chain2_id = chain2_id
        self.res_range1 = res_range1
        self.res_range2 = res_range2
        self.core_filter1 = core_filter1
        self.core_region1 = core_region1
        self.core_filter2 = core_filter2
        self.core_region2 = core_region2

    def add_2nd_oligomer(
        self,
        pdb_filename_2nd="",
        chain1bis_id="",
        chain2bis_id="",
        res_range1bis="",
        res_range2bis="",
        core_filter1bis="all",
        core_region1bis="",
        core_filter2bis="all",
        core_region2bis="",
    ):
        """Fill the values of the 2nd oligomer form data."""
        self.second_oligomer = True
        self.pdb_filename_2nd = pdb_filename_2nd
        self.chain1bis_id = chain1bis_id
        self.chain2bis_id = chain2bis_id
        self.res_range1bis = res_range1bis
        self.res_range2bis = res_range2bis
        self.core_filter1bis = core_filter1bis
        self.core_region1bis = core_region1bis
        self.core_filter2bis = core_filter2bis
        self.core_region2bis = core_region2bis

    def __repr__(self):
        return f"""<Id {self.request_id}, pdb_filename: {self.pdb_filename}"""
