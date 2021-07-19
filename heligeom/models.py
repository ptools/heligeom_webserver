
from . import db


class UserInputs(db.Model):
    request_id = db.Column(db.String(32), unique=True, primary_key=True)
    pdb_filename = db.Column(db.String(100), nullable=False)
    chain1_id = db.Column(db.String(2), default="A")
    chain2_id = db.Column(db.String(2), default="B")
    res_range1 = db.Column(db.String(50), default="")
    res_range2 = db.Column(db.String(50), default="")
    n_mer = db.Column(db.Integer, default=0)
    z_align = db.Column(db.Boolean, default=False)

    def __init__(self, request_id, pdb_filename="", chain1_id="A", chain2_id="B", res_range1="", res_range2="",
                 n_mer=0, z_align=False):
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
