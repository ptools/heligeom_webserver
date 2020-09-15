
from . import db


class User_Inputs(db.Model):
    request_id = db.Column(db.String(32), unique=True, primary_key=True)
    pdb1_filename = db.Column(db.String(100), nullable=False)
    pdb2_filename = db.Column(db.String(100), nullable=False)
    n_mer = db.Column(db.Integer, default=0)
    z_align = db.Column(db.Boolean, default=False)

    def __init__(self, request_id, pdb1_filename, pdb2_filename, n_mer=0, z_align=False):
        self.request_id = request_id
        self.pdb1_filename = pdb1_filename
        self.pdb2_filename = pdb2_filename
        self.n_mer = n_mer
        self.z_align = z_align


    def __repr__(self):
        return '<Id {0}, pdb_file: {1}, n_mer: {2}, z_align: {3}>'.format(self.request_id, self.pdb1_filename, self.n_mer, self.z_align)
