
from heligeom import db


class Form_Data(db.Model):
    request_id = db.Column(db.String(32), unique=True, primary_key=True)
    pdb1_filename = db.Column(db.String(100), nullable=False)
    pdb2_filename = db.Column(db.String(100), nullable=False)
    n_mer = db.Column(db.Integer)
    z_align = db.Column(db.Boolean)

    def __repr__(self):
        return '<Id {0}, pdb_file: {1}, n_mer: {2}, z_align: {3}>'.format(self.request_id, self.pdb1_filename, self.z_align)
