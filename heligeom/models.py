from . import db


class UserInputs(db.Model):
    """Class defining the data stored in the Database.

    Matches the field of the `HeligeomForm`.
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
    # Number of copy/monomers requested to create the filament.
    n_mer = db.Column(db.Integer, default=0)
    # Z alignement of the filament (Yes/No)
    z_align = db.Column(db.Boolean, default=False)

    # For a second oligomer
    second_oligomer = db.Column(db.Boolean, default=False)
    pdb_filename_2nd = db.Column(db.String(100), default="")
    chain1bis_id = db.Column(db.String(2), default="")
    chain2bis_id = db.Column(db.String(2), default="")
    res_range1bis = db.Column(db.String(50), default="")
    res_range2bis = db.Column(db.String(50), default="")

    def __init__(
        self,
        request_id,
        pdb_filename="",
        chain1_id="",
        chain2_id="",
        res_range1="",
        res_range2="",
        n_mer=0,
        z_align=False,
    ):
        self.request_id = request_id
        self.pdb_filename = pdb_filename
        self.chain1_id = chain1_id
        self.chain2_id = chain2_id
        self.res_range1 = res_range1
        self.res_range2 = res_range2
        self.n_mer = n_mer
        self.z_align = z_align

    def add_2nd_oligomer(
        self,
        pdb_filename_2nd="",
        chain1bis_id="",
        chain2bis_id="",
        res_range1bis="",
        res_range2bis="",
    ):
        self.second_oligomer = True
        self.pdb_filename_2nd = pdb_filename_2nd
        self.chain1bis_id = chain1bis_id
        self.chain2bis_id = chain2bis_id
        self.res_range1bis = res_range1bis
        self.res_range2bis = res_range2bis

    def __repr__(self):
        return f"""<Id {self.request_id}, pdb_filename: {self.pdb_filename},
                n_mer: {self.n_mer}, z_align: {self.z_align}>"""
