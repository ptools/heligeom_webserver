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
    chain1_id = db.Column(db.String(2), default="A")
    # Chain ID to define the 2nd monomer
    chain2_id = db.Column(db.String(2), default="B")
    # Residue range to define the 1st monomer
    res_range1 = db.Column(db.String(50), default="")
    # Residue range to define the 2nd monomer
    res_range2 = db.Column(db.String(50), default="")
    # Number of copy/monomers requested to create the filament.
    n_mer = db.Column(db.Integer, default=0)
    # Z alignement of the filament (Yes/No)
    z_align = db.Column(db.Boolean, default=False)


    def __init__(self, request_id, pdb_filename="", chain1_id="A", chain2_id="B",
                 res_range1="", res_range2="", n_mer=0, z_align=False):
        self.request_id = request_id
        self.pdb_filename = pdb_filename
        self.chain1_id = chain1_id
        self.chain2_id = chain2_id
        self.res_range1 = res_range1
        self.res_range2 =res_range2
        self.n_mer = n_mer
        self.z_align = z_align


    def __repr__(self):
        return '<Id {0}, pdb_filename: {1}, n_mer: {2}, z_align: {3}>'.format(self.request_id, self.pdb_filename, self.n_mer, self.z_align)
