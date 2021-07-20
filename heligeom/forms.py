"""
Module to handle the form to run Heligeom.
Contains also functions which verify the field of the form.
"""

import urllib.request

from flask import current_app
from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileAllowed
from wtforms import IntegerField, BooleanField, StringField, validators
from werkzeug.utils import secure_filename


class HeligeomForm(FlaskForm):
    """Class to handle the form to run Heligeom."""

    # Input Structure: either a PDB file or a PDB ID
    input_file = FileField('input_file', validators=[FileAllowed(['pdb'], 'Only a PDB file can be uploaded.')])
    pdb_id = StringField('pdb_id', validators=[validators.Optional(),
                                                 validators.length(min=4, max=4)])

    # Chain ID to select a 1st monomer
    chain1_id = StringField('chain1_id', validators=[validators.Optional(),
                                                 validators.length(min=1, max=2)])
    # Chain ID to select a 2nd monomer
    chain2_id = StringField('chain2_id', validators=[validators.Optional(),
                                                 validators.length(min=1, max=2)])

    # A residue range to select a 1st monomer
    res_range1 = StringField('res_range1', validators=[validators.Optional(),
                                                 validators.length(min=1, max=50)])
    # A residue range to select a 2nd monomer
    res_range2 = StringField('res_range2', validators=[validators.Optional(),
                                                 validators.length(min=1, max=50)])

    # Number of copy/monomers requested to create the filament.
    n_mer = IntegerField("n_mer", validators=[validators.Optional(),
                                              validators.NumberRange(0, 50, message="Only a Number.")])
    # If the filament will be align on Z.
    z_align = BooleanField("z_align")


    def validate_screw(self, extra_validators=None):
        """Custom validate method for a part of the FlaskForm.
        It checks only the field needed to retrieve the screw parameters.

        Fill the `error` attribute of a field with a erro message if it's not valid.

        Returns
        -------
        Bool
            True if the needed field are valid. False otherwise.
        """

        # Start by calling the parent method
        if not super().validate(extra_validators=extra_validators):
            return False

        # Now, check for each structure, if the users filled the upload part or the pdb id
        if self.input_file.data is None and self.pdb_id.data == '':
            self.input_file.errors =  self.pdb_id.errors = "No Data specified"
            return False

        return True


    def validate(self, extra_validators=None):
        """Overload validate() method of the FlaskForm.
        Check all fields.

        Fill the `error` attribute of a field with a erro message if it's not valid.

        Returns
        -------
        Bool
            True if all field forms are valid. False otherwise.
        """
        # Call the custom one
        if not self.validate_screw(extra_validators=extra_validators):
            return False

        # Validate the others inputs
        if self.n_mer.data is None:
            self.n_mer.errors = "Required to construct the filament."
            return False

        return True


def validate_input_structure(pdb_file, pdb_id, path):
    """Check for the input structure if the user submit an uploaded file or a PDB id.
    The user uploaded file or the file downloaded corresponding to the PDB ID is stored inside.

    Parameters
    ----------
    pdb_file : FileField
        File Field of the Flask form which represent the upload PDB file.
    pdb_id : StringField
        String field of the flask form which represent the PDB ID (4 letters)
    folder_name : str
        Name of the folder where to store the PDB file.
    path: pathlib.Path
        Path to save the file.

    Returns
    -------
    str
        the filename of the PDB.
    """

    # Upload ?
    if pdb_file.data:
        pdb_file = pdb_file.data
        pdb_filename = secure_filename(pdb_file.filename)
        if pdb_filename == '':
            return None
        pdb_file.save(path/pdb_filename)
        return pdb_filename

    # PDB ID
    # Generate url
    pdb_filename = pdb_id.data + ".pdb"
    url = current_app.config["PDB_SERVER"] + pdb_filename
    try:
        urllib.request.urlretrieve(url, path/pdb_filename)
        return pdb_filename
    except urllib.error.URLError:
        pdb_id.errors = "PDB ID unknown"
        return None
