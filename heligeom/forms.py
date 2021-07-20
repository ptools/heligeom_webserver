
from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileAllowed
from wtforms import IntegerField, BooleanField, StringField, validators


class TestForm(FlaskForm):

    input_file = FileField('input_file', validators=[FileAllowed(['pdb'], 'Only a PDB file can be uploaded')])
    pdb_id = StringField('pdb_id', validators=[validators.Optional(),
                                                 validators.length(min=4, max=4)])
    chain1_id = StringField('chain1_id', validators=[validators.Optional(),
                                                 validators.length(min=1, max=2)])
    chain2_id = StringField('chain2_id', validators=[validators.Optional(),
                                                 validators.length(min=1, max=2)])
    res_range1 = StringField('res_range1', validators=[validators.Optional(),
                                                 validators.length(min=1, max=50)])
    res_range2 = StringField('res_range2', validators=[validators.Optional(),
                                                 validators.length(min=1, max=50)])

    n_mer = IntegerField("n_mer", validators=[validators.Optional(),
                                              validators.NumberRange(0, 50, message="Only a Number")])
    z_align = BooleanField("z_align")


    # Custom validate method for a part of the FlaskForm
    def validate_screw(self, extra_validators=None):

        # Start by calling the parent method
        if not super().validate(extra_validators=extra_validators):
            return False

        # Now, check for each structure, if the users filled the upload part or the pdb id
        if self.input_file.data is None and self.pdb_id.data == '':
            self.input_file.errors =  self.pdb_id.errors = "No Data specified"
            return False

        return True

    # Overload validate() method of the FlaskForm
    def validate(self, extra_validators=None):

        # Call the custom one
        if not self.validate_screw(extra_validators=extra_validators):
            return False

        # Validate the others inputs
        if self.n_mer.data is None:
            self.n_mer.errors = "Required to construct the filament."
            return False

        return True
