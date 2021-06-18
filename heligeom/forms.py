
from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileAllowed
from wtforms import IntegerField, BooleanField, StringField, validators


class TestForm(FlaskForm):

    pdb1_file = FileField('pdb1_file', validators=[FileAllowed(['pdb'], 'Only a PDB file')])
    pdb2_file = FileField('pdb2_file', validators=[FileAllowed(['pdb'], 'Only a PDB file')])
    pdb1_id = StringField('pdb1_id', validators=[validators.Optional(),
                                                 validators.length(min=4, max=4)])
    pdb2_id = StringField('pdb1_id', validators=[validators.Optional(),
                                                 validators.length(min=4, max=4)])
    n_mer = IntegerField("n_mer", validators=[validators.InputRequired(message=""),
                                              validators.NumberRange(0, message="Only a Number")])
    z_align = BooleanField("z_align")


    #Overload validate() method of the FlaskForm
    def validate(self, extra_validators=None):

        # Start by calling the parent method
        if not super().validate(extra_validators=extra_validators):
            return False

        # Now, check for each structure, if the users filled the upload part or the pdb id
        if self.pdb1_file.data is None and self.pdb1_id.data == '':
            self.pdb1_file.errors =  self.pdb1_id.errors = "No Data specified"
            return False
        if self.pdb2_file.data is None and self.pdb2_id.data == '':
            self.pdb2_file.errors =  self.pdb2_id.errors = "No Data specified"
            return False

        return True
