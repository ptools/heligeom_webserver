"""
Module to handle the form to run Heligeom.
Contains also functions which verify the field of the form.
"""

import re
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

    #: Regular epxression for the chain
    cls_regexp_chain = r"^[A-Za-z]+$"
    #: Regular expression for selecting the residue range
    cls_regexp_resrange = r"^\d+-\d+$"
    #: Regular expression for selecting a core region
    cls_regexp_core = r"^(\d+-\d+)(\s*,\s*\d+-\d+)*\s*$"

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
        validators=[
            validators.Optional(),
            validators.length(min=1, max=50),
            validators.Regexp("^[0-9]+-[0-9]+$", message="The range is invalid."),
        ],
        render_kw={"placeholder": "1-300"},
    )
    # A residue range to select a 2nd monomer
    res_range2 = StringField(
        "res_range2",
        validators=[
            validators.Optional(),
            validators.length(min=1, max=50),
            validators.Regexp("^[0-9]+-[0-9]+$", message="The range is invalid."),
        ],
        render_kw={"placeholder": "301-600"},
    )

    # A residue range to select an optional variant of 1st monomer
    res_range1bis = StringField(
        "res_range1bis",
        validators=[
            validators.Optional(),
            validators.length(min=1, max=50),
            validators.Regexp("^[0-9]+-[0-9]+$", message="The range is invalid."),
        ],
        render_kw={"placeholder": "20-200"},
    )
    # A residue range to select an optional variant of 2nd monomer
    res_range2bis = StringField(
        "res_range2bis",
        validators=[
            validators.Optional(),
            validators.length(min=1, max=50),
            validators.Regexp("^[0-9]+-[0-9]+$", message="The range is invalid."),
        ],
        render_kw={"placeholder": "220-400"},
    )

    core_filter1 = RadioField(
        "core_filter1",
        choices=[("all", "All"), ("manual", "Select core residues:")],
        default="all",
    )
    core_region1 = StringField(
        "core_region1",
        validators=[
            validators.Optional(),
            validators.length(min=1, max=50),
        ],
        render_kw={"placeholder": "1-30, 60-100"},
    )

    core_filter2 = RadioField(
        "core_filter2",
        choices=[("all", "All"), ("manual", "Select core residues:")],
        default="all",
    )
    core_region2 = StringField(
        "core_region2",
        validators=[
            validators.Optional(),
            validators.length(min=1, max=50),
        ],
        render_kw={"placeholder": "100-130, 160-200"},
    )

    core_filter1bis = RadioField(
        "core_filter1bis",
        choices=[("all", "All"), ("manual", "Select core residues:")],
        default="all",
    )
    core_region1bis = StringField(
        "core_region1bis",
        validators=[
            validators.Optional(),
            validators.length(min=1, max=50),
        ],
        render_kw={"placeholder": "1-30, 60-100"},
    )

    core_filter2bis = RadioField(
        "core_regions_choice2bis",
        choices=[("all", "All"), ("manual", "Select core residues:")],
        default="all",
    )
    core_region2bis = StringField(
        "core_region2bis",
        validators=[
            validators.Optional(),
            validators.length(min=1, max=50),
        ],
        render_kw={"placeholder": "100-130, 160-200"},
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
        if not utils.check_pair_monomers(
            self.chain1_id, self.res_range1, self.chain2_id, self.res_range2
        ):
            return False

        # Core regions options need to be identical between the 2 monomers
        if self.core_filter1.data != self.core_filter2.data:
            self.core_region1.errors = self.core_region2.errors = [
                "Both core region selections needs to be filled."
            ]
            return False

        # Ensure the core region selections for both monomers are correct
        return self.check_core_regions()

    def check_core_regions(self):
        return check_core_region_monomer(
            self.core_filter1, self.core_region1
        ) and check_core_region_monomer(self.core_filter2, self.core_region2)

    def check_core_regions_bis(self):
        return check_core_region_monomer(
            self.core_filter1bis, self.core_region1bis
        ) and check_core_region_monomer(self.core_filter2bis, self.core_region2bis)

    def validate(self, extra_validators=None):
        """Overload validate() method of the FlaskForm.
        Check all fields.

        Fill the `error` attribute of a field with a error message if it's not valid.

        Returns
        -------
        Bool
            True if all field forms are valid. False otherwise.
        """

        # Start by calling the parent method
        if not super().validate(extra_validators=extra_validators):
            return False

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
        if not utils.check_pair_monomers(
            self.chain1bis_id,
            self.res_range1bis,
            self.chain2bis_id,
            self.res_range2bis,
        ):
            return False

        # Core regions options need to be identical between the 2 monomers
        if self.core_filter1bis.data != self.core_filter2bis.data:
            self.core_region1bis.errors = self.core_region2bis.errors = [
                "Both core region selections needs to be filled."
            ]
            return False

        return self.check_core_regions_bis()

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


def check_core_region_monomer(filter, string_range):
    # No core region defined
    if filter.data == "all":
        return True

    # The radio button is set to manual
    if not string_range.data:
        string_range.errors = ["This field cannot be empty."]
        return False

    res = re.match(InputStructures.cls_regexp_core, string_range.data)
    if not res:
        string_range.errors = ["The range is invalid."]
        return False

    for res_range in string_range.data.split(","):
        try:
            utils.parse_resrange(res_range)
        except SyntaxError as e:
            string_range.errors = [e.msg]  # Use a list because it's unpacked in HTML.
            return False

    return True


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
