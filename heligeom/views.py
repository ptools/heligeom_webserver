
import os
import pathlib
import uuid
import math
import urllib.request

from flask import Blueprint, current_app, render_template, redirect, url_for, request
from werkzeug.utils import secure_filename

from .forms import TestForm
from .models import db, UserInputs
from .tool import run

# Blueprint Configuration
heligeom_bp = Blueprint(
    'heligeom_bp', __name__,
    template_folder='templates',
    static_folder='static'
)

def validate_input_structure(pdb_file, pdb_id, folder_name):
    """
    Check for the input structure if the user submit an uploaded file or a PDB id.
    """

    # Upload
    if pdb_file.data:
        pdb_file = pdb_file.data
        pdb_filename = secure_filename(pdb_file.filename)
        if pdb_filename == '':
            return None
        pdb_file.save(os.path.join(current_app.config["DATA_UPLOADS"], folder_name, pdb_filename))
        return pdb_filename
    # or PDB id ?
    else:
        pdb_filename = pdb_id.data + ".pdb"
        url = current_app.config["PDB_SERVER"] + pdb_filename
        path = os.path.join(current_app.config["DATA_UPLOADS"], folder_name, pdb_filename)
        try:
            urllib.request.urlretrieve(url, path)
            return pdb_filename
        except urllib.error.URLError:
            pdb_id.errors = "PDB id unknown"
            return None


@heligeom_bp.route('/')
def homepage():
    return render_template('index.html')

@heligeom_bp.route('/run', methods=['GET', 'POST'])
def runpage():

    # Initialize submission form
    form = TestForm()
    # Validation part
    if form.validate():
        if request.method == 'POST':

            #Generate UUID for page results
            uniq_id = uuid.uuid4().hex
            # Create result folder

            print("path :", current_app.config['DATA_UPLOADS'], uniq_id)
            pathlib.Path(current_app.config['DATA_UPLOADS'], uniq_id).mkdir(exist_ok=True)


            pdb1_filename = validate_input_structure(form.pdb1_file, form.pdb1_id, uniq_id)
            pdb2_filename = validate_input_structure(form.pdb2_file, form.pdb2_id, uniq_id)
            if pdb1_filename and pdb2_filename:

                # Save to the database the form inputs
                # Only way for now to pass the form data to another page.
                # We could use session or flash messages but neither seems to fit the need.
                user_inputs = UserInputs(request_id=uniq_id,
                                         pdb1_filename=pdb1_filename, pdb2_filename=pdb2_filename,
                                         n_mer=form.n_mer.data, z_align=form.z_align.data)
                db.session.add(user_inputs)
                db.session.commit()

                #Redirect to the results page with correct id
                return redirect(url_for('heligeom_bp.results', results_id=uniq_id))
    return render_template('run.html', form=form)

@heligeom_bp.route("/results/<results_id>", methods=['GET', 'POST'])
def results(results_id):

    # Query the database to retrieve the form inputs
    query_result = db.session.query(UserInputs).filter(UserInputs.request_id == results_id).first()
    pdb1_filename = query_result.pdb1_filename
    pdb2_filename = query_result.pdb2_filename
    n_mer = query_result.n_mer
    z_align = query_result.z_align

    #Ptools part
    pdb_out_name = "out.pdb"
    pdb_out_abs_path = os.path.join(current_app.config["DATA_UPLOADS"], results_id, pdb_out_name)

    pdb1_abs_path = os.path.join(current_app.config["DATA_UPLOADS"], results_id, pdb1_filename)
    pdb2_abs_path = os.path.join(current_app.config["DATA_UPLOADS"], results_id, pdb2_filename)
    #Run the Heligeom calculations and write the PDB result in pdb_out_abs_path
    # Return the helicoidal parameters
    hp, pitch, nb_monomers, direction = run(pdb1_abs_path, pdb2_abs_path, n_mer, pdb_out_abs_path)

    # Construct a relative path for the pdb result to be givent to the liteMol plugin on the html page.
    pdb_out_webpath = pathlib.Path('static/data', results_id, pdb_out_name)

    # Create dict of data to pass to render_template
    data = {
        "pdb_path" : pdb_out_webpath,
        "n_mer": n_mer,
        "z_align": z_align,
        "pitch": f"{pitch:3.2f}",
        "nb_monomers": f"{nb_monomers:3.2f}",
        "direction": direction,
        "rotation_angle": f"{hp.angle:3.2f}",
        "rotation_angle_deg": f"{math.degrees(hp.angle):3.2f}",
        "axis_point": [f"{coord:3.2f}" for coord in hp.point],
        "translation": f"{hp.normtranslation:3.2f}"
    }

    return render_template('results.html',data=data)

@heligeom_bp.route('/help')
def help():
    return render_template('help.html')

@heligeom_bp.route('/contact')
def contact():
    return render_template('contact.html')

