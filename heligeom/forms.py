
from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileAllowed
from wtforms import IntegerField, BooleanField, StringField, validators


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
