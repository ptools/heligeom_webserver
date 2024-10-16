"""
Module to handle the form to run Heligeom.
Contains also functions which verify the field of the form.
"""

import urllib.error
import urllib.request

from flask import current_app
from flask_wtf import FlaskForm
from flask_wtf.file import FileAllowed, FileField
from werkzeug.utils import secure_filename
from wtforms import BooleanField, IntegerField, RadioField, StringField, validators

from . import utils


class InputStructures(FlaskForm):
    """Class to handle the input structures form."""

    # Input Structure: either a PDB file or a PDB ID
    input_file = FileField(
        "input_file",
        validators=[FileAllowed(["pdb"], "Only a PDB file can be uploaded.")],
    )
    pdb_id = StringField(
        "pdb_id",
        validators=[validators.Optional(), validators.length(min=4, max=4)],
        render_kw={"placeholder": "2GLS"},
    )

    # Input Structure: either a PDB file or a PDB ID
    input_file_2nd = FileField(
        "input_file",
        validators=[FileAllowed(["pdb"], "Only a PDB file can be uploaded.")],
    )
    pdb_id_2nd = StringField(
        "pdb_id",
        validators=[validators.Optional(), validators.length(min=4, max=4)],
        render_kw={"placeholder": "2GLS"},
    )

    # Chain ID to select a 1st monomer
    chain1_id = StringField(
        "chain1_id",
        validators=[
            validators.Optional(),
            validators.Regexp("^[A-Za-z]+$", message="Chain ID must contain only letters."),
            validators.length(min=1, max=2, message="Chain ID must be between 1 or 2 letters."),
        ],
        render_kw={"placeholder": "A"},
    )
    # Chain ID to select a 2nd monomer
    chain2_id = StringField(
        "chain2_id",
        validators=[
            validators.Optional(),
            validators.Regexp("^[A-Za-z]+$", message="Chain ID must contain only letters."),
            validators.length(min=1, max=2, message="Chain ID must be between 1 or 2 letters."),
        ],
        render_kw={"placeholder": "B"},
    )

    # Chain ID to select an optional variant of 1st monomer
    chain1bis_id = StringField(
        "chain1bis_id",
        validators=[
            validators.Optional(),
            validators.Regexp("^[A-Za-z]+$", message="Chain ID must contain only letters."),
            validators.length(min=1, max=2, message="Chain ID must be between 1 or 2 letters."),
        ],
        render_kw={"placeholder": "C"},
    )
    # Chain ID to select an optional variant of 2nd monomer
    chain2bis_id = StringField(
        "chain2bis_id",
        validators=[
            validators.Optional(),
            validators.Regexp("^[A-Za-z]+$", message="Chain ID must contain only letters."),
            validators.length(min=1, max=2, message="Chain ID must be between 1 or 2 letters."),
        ],
        render_kw={"placeholder": "D"},
    )

    # A residue range to select a 1st monomer
    res_range1 = StringField(
        "res_range1",
        validators=[validators.Optional(), validators.length(min=1, max=50)],
        render_kw={"placeholder": "1-300"},
    )
    # A residue range to select a 2nd monomer
    res_range2 = StringField(
        "res_range2",
        validators=[validators.Optional(), validators.length(min=1, max=50)],
        render_kw={"placeholder": "301-600"},
    )

    # A residue range to select an optional variant of 1st monomer
    res_range1bis = StringField(
        "res_range1bis",
        validators=[validators.Optional(), validators.length(min=1, max=50)],
        render_kw={"placeholder": "20-200"},
    )
    # A residue range to select an optional variant of 2nd monomer
    res_range2bis = StringField(
        "res_range2bis",
        validators=[validators.Optional(), validators.length(min=1, max=50)],
        render_kw={"placeholder": "220-400"},
    )

    core_regions = RadioField(
        "core_regions",
        choices=[("all", "toto"), ("Select core residues", "toto")],
    )

    def validate_1st_oligomer(self):
        """Custom validate method for InputStructures form.
        It checks the field needed to retrieve the screw parameters.

        Fill the `error` attribute of a field with an error message if it's not valid.

        Check only fields for the first oligomer.

        Returns
        -------
        Bool
            True if the needed fields are valid. False otherwise.
        """
        # Now, check for each structure, if the users filled the upload part or the pdb id
        if self.input_file.data is None and self.pdb_id.data == "":
            self.input_file.errors = self.pdb_id.errors = "No Data specified"
            return False

        # Ensure monomer selections are correct
        return utils.check_pair_monomers(
            self.chain1_id, self.res_range1, self.chain2_id, self.res_range2
        )

    def validate(self, extra_validators=None):
        """Overload validate() method of the FlaskForm.
        Check all fields.

        Fill the `error` attribute of a field with a error message if it's not valid.

        Returns
        -------
        Bool
            True if all field forms are valid. False otherwise.
        """
        if not self.validate_1st_oligomer():
            return False

        if self.has_2nd_oligomer():
            if not self.validate_2nd_oligomer():
                return False

        return True

    def validate_2nd_oligomer(self):
        """Check if the form inputs for a 2nd assembly
        (a different selection of 2 monomers) is valid.

        It doesn't test the PDB entry (file or ID) because blank
        means the first entry is reused.

        Returns
        -------
        Bool
            True if 2nd oligomer fields forms are valid. False otherwise.
        """
        # Ensure monomer selections are correct
        return utils.check_pair_monomers(
            self.chain1bis_id,
            self.res_range1bis,
            self.chain2bis_id,
            self.res_range2bis,
        )

    def has_2nd_oligomer(self):
        """Check if the form fields for a 2nd assembly
           (a different selection of 2 monomers) is present.

        Returns
        -------
        Bool
            True if at least one field is not None.
        """

        fields = [
            self.chain1bis_id.data,
            self.chain2bis_id.data,
            self.res_range1bis.data,
            self.res_range2bis.data,
            self.pdb_id_2nd.data,
            self.input_file_2nd.data,
        ]

        for field in fields:
            if field:
                return True
        return False


def validate_input_structure(pdb_file, pdb_id, path):
    """Check for the input structure if the user submit an uploaded file or a PDB id.
    The user uploaded file or the file downloaded corresponding to the PDB ID is stored.

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
        the filename of the PDB. None otherwise
    """

    # If there is no data in the fields.
    if not pdb_file.data and not pdb_id.data:
        return None

    # Upload ?
    if pdb_file.data:
        pdb_file = pdb_file.data
        pdb_filename = secure_filename(pdb_file.filename)
        if pdb_filename == "":
            return None
        pdb_file.save(path / pdb_filename)
        return pdb_filename

    # PDB ID
    # Generate url
    pdb_filename = pdb_id.data + ".pdb"
    url = current_app.config["PDB_SERVER"] + pdb_filename
    try:
        urllib.request.urlretrieve(url, path / pdb_filename)
        return pdb_filename
    except urllib.error.URLError:
        pdb_id.errors = "PDB ID unknown"
        return None


class Construction(FlaskForm):
    """Class to handle the construction parameters form."""

    # Number of copy/monomers requested to create the filament.
    n_mer = IntegerField(
        "n_mer",
        validators=[
            validators.NumberRange(0, 100, message="Only a number is accepted."),
        ],
        render_kw={"placeholder": "10"},
    )

    # If the filament will be align on Z.
    z_align = BooleanField("z_align", validators=[validators.Optional()])

    def validate(self, extra_validators=None):
        """Overload validate() method of the Construction FlaskForm.
        Check the n_mer field.

        Fill the `error` attribute of a field with a erro message if it's not valid.

        Returns
        -------
        Bool
            True if the fields forms are valid. False otherwise.
        """

        if self.n_mer.data is None:
            self.n_mer.errors = "Required to construct the filament."
            return False

        return True
