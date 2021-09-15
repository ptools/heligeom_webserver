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

from . import utils

class HeligeomForm(FlaskForm):
    """Class to handle the form to run Heligeom."""

    # Input Structure: either a PDB file or a PDB ID
    input_file = FileField('input_file', validators=[
        FileAllowed(['pdb'], 'Only a PDB file can be uploaded.')])
    pdb_id = StringField('pdb_id', validators=[
        validators.Optional(),
        validators.length(min=4, max=4)])

    # Input Structure: either a PDB file or a PDB ID
    input_file_2nd = FileField('input_file', validators=[
        FileAllowed(['pdb'], 'Only a PDB file can be uploaded.')])
    pdb_id_2nd = StringField('pdb_id', validators=[
        validators.Optional(),
        validators.length(min=4, max=4)])

    # Chain ID to select a 1st monomer
    chain1_id = StringField('chain1_id', validators=[
        validators.Optional(),
        validators.Regexp('^[A-Za-z]+$', message="Chain ID must contain only letters."),
        validators.length(min=1, max=2, message="Chain ID must be between 1 or 2 letters.")])
    # Chain ID to select a 2nd monomer
    chain2_id = StringField('chain2_id', validators=[
        validators.Optional(),
        validators.Regexp('^[A-Za-z]+$', message="Chain ID must contain only letters."),
        validators.length(min=1, max=2, message="Chain ID must be between 1 or 2 letters.")])

    # Chain ID to select an optional variant of 1st monomer
    chain1bis_id = StringField('chain1bis_id', validators=[
        validators.Optional(),
        validators.Regexp('^[A-Za-z]+$', message="Chain ID must contain only letters."),
        validators.length(min=1, max=2, message="Chain ID must be between 1 or 2 letters.")])
    # Chain ID to select an optional variant of 2nd monomer
    chain2bis_id = StringField('chain2bis_id', validators=[
        validators.Optional(),
        validators.Regexp('^[A-Za-z]+$', message="Chain ID must contain only letters."),
        validators.length(min=1, max=2, message="Chain ID must be between 1 or 2 letters.")])

    # A residue range to select a 1st monomer
    res_range1 = StringField('res_range1', validators=[
        validators.Optional(),
        validators.length(min=1, max=50)])
    # A residue range to select a 2nd monomer
    res_range2 = StringField('res_range2', validators=[
        validators.Optional(),
        validators.length(min=1, max=50)])

    # A residue range to select an optional variant of 1st monomer
    res_range1bis = StringField('res_range1bis', validators=[
        validators.Optional(),
        validators.length(min=1, max=50)])
    # A residue range to select an optional variant of 2nd monomer
    res_range2bis = StringField('res_range2bis', validators=[
        validators.Optional(),
        validators.length(min=1, max=50)])


    # Number of copy/monomers requested to create the filament.
    n_mer = IntegerField("n_mer", validators=[
        validators.Optional(),
        validators.NumberRange(0, 50, message="Only number is accepted.")])
    # If the filament will be align on Z.
    z_align = BooleanField("z_align")


    def validate_screw(self, extra_validators=None):
        """Custom validate method for a part of the HeligeomForm.
        It checks only the field needed to retrieve the screw parameters.

        Fill the `error` attribute of a field with an error message if it's not valid.

        Returns
        -------
        Bool
            True if the needed fields are valid. False otherwise.
        """

        # Start by calling the parent method
        if not super().validate(extra_validators=extra_validators):
            return False

        # Now, check for each structure, if the users filled the upload part or the pdb id
        if self.input_file.data is None and self.pdb_id.data == '':
            self.input_file.errors =  self.pdb_id.errors = "No Data specified"
            return False

        # Validation of the monomer selections
        # For each monomer, a valid selection is either a chain attribute or/and a residue range
        monomer_validation = True # Use a boolean to test both monomer selections at the same time
        for chain_id, res_range in zip([self.chain1_id, self.chain2_id], [self.res_range1, self.res_range2]):
            # Check at least one type of selection is chosen
            if not chain_id.data and not res_range.data:
                chain_id.errors = res_range.errors = "No Data specified"
                monomer_validation = False
            elif res_range.data:
                # Check if the residue range is consistent
                try:
                    utils.parse_resrange(res_range.data)
                except SyntaxError as e:
                    res_range.errors = e.msg
                    monomer_validation = False

        return monomer_validation


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

        if not self.validate_2nd_oligomer(extra_validators=extra_validators):
            return False

        # Validate the others inputs
        if self.n_mer.data is None:
            self.n_mer.errors = "Required to construct the filament."
            return False

        return True

    def validate_2nd_oligomer(self, extra_validators=None):
        """Check if the form inputs for a 2nd assembly (a different selection of 2 monomers) is valide.
        """

        # Start by calling the parent method
        if not super().validate(extra_validators=extra_validators):
            return False

        # Validation of the monomer selections
        # For each monomer, a valid selection is either a chain attribute or/and a residue range
        monomer_validation = True # Use a boolean to test both monomer selections at the same time
        monomer_data = {1: True, 2: True}
        number_monomer = 1
        for chain_id, res_range in zip([self.chain1bis_id, self.chain2bis_id], [self.res_range1bis, self.res_range2bis]):
            # Check at least one type of selection is chosen
            if not chain_id.data and not res_range.data:
                monomer_data[number_monomer] = False
            elif res_range.data:
                # Check if the residue range is consistent
                try:
                    utils.parse_resrange(res_range.data)
                except SyntaxError as e:
                    res_range.errors = e.msg
                    monomer_validation = False
            number_monomer += 1

        #Determine if data has been entered
        if monomer_data[1] == monomer_data[2] == False:
            monomer_validation = True
        elif not monomer_data[1]:
            self.chain1bis_id.errors = self.res_range1bis.errors =  "No Data specified"
            monomer_validation = False
        elif not monomer_data[2]:
            self.chain2bis_id.errors = self.res_range2bis.errors =  "No Data specified"
            monomer_validation = False

        return monomer_validation


    def has_2nd_oligomer(self):
        """Check if the form fields for a 2nd assembly (a different selection of 2 monomers) is present.

        Returns
        -------
        Bool
            True if at least one field is not None.
        """

        fields = [self.chain1bis_id.data, self.chain2bis_id.data, self.res_range1bis.data, self.res_range2bis.data, self.pdb_id_2nd.data, self.input_file_2nd.data]

        for field in fields:
            if field:
                return True
        return False



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

    # If there is no data in the fields.
    if not pdb_file.data and not pdb_id.data:
        return None

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
