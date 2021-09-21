
import pathlib
import uuid
import math
import traceback

from flask import Blueprint, current_app, render_template, redirect, url_for, request, send_from_directory

from .forms import HeligeomForm, validate_input_structure
from .models import db, UserInputs
from .tool import construct, screw_parameters, analyze_fnat

# Blueprint Configuration
heligeom_bp = Blueprint(
    'heligeom_bp', __name__,
    template_folder='templates',
    static_folder='static'
)


@heligeom_bp.route('/')
def homepage():
    return render_template('index.html')

@heligeom_bp.route('/run', methods=['GET', 'POST'])
def runpage():

    # Initialize submission form
    form = HeligeomForm()

    if request.method == 'POST':
        # Check wether only the screw parameters or a construction has been requested

        # Screw parameters asked
        if request.form.get("action") == "screw" and form.validate_screw() and form.validate_2nd_oligomer():

            #Generate UUID for storing files
            uniq_id = uuid.uuid4().hex
            # Create result folder
            result_path = pathlib.Path(current_app.config['DATA_UPLOADS'], uniq_id)
            result_path.mkdir(exist_ok=True)

            pdb_filename = validate_input_structure(form.input_file, form.pdb_id, result_path)
            if pdb_filename:

                try:
                    pdb_abs_path = result_path / pdb_filename
                    # Compute Helicoidal parameters
                    hp, pitch, monomers_per_turn, direction, dmin, dmax = screw_parameters(pdb_abs_path, form.chain1_id.data, form.chain2_id.data,
                                                                            form.res_range1.data, form.res_range2.data)

                    screw_data = {
                        "pitch": f"{pitch:3.2f}",
                        "nb_monomers": f"{monomers_per_turn:3.2f}",
                        "direction": direction,
                        "rotation_angle": f"{hp.angle:3.2f}",
                        "rotation_angle_deg": f"{math.degrees(hp.angle):3.2f}",
                        "axis_point": [f"{coord:3.2f}" for coord in hp.point],
                        "translation": f"{hp.normtranslation:3.2f}",
                        "dmin"      : f"{dmin:3.2f}",
                        "dmax"       : f"{dmax:3.2f}"
                    }


                    # Does the user input for a 2nd assembly ?
                    if form.has_2nd_oligomer():

                        # A second structure has been provided?
                        # If not, use the pdbfile from the first oligomeric form
                        pdb_filename_2nd = validate_input_structure(form.input_file_2nd, form.pdb_id_2nd, result_path)
                        if(pdb_filename_2nd):
                            pdb_abs_path = result_path / pdb_filename_2nd

                        # Compute Helicoidal parameters for 2nd assembly
                        hp_bis, pitch_bis, monomers_per_turn_bis, direction_bis, dminbis, dmaxbis = screw_parameters(pdb_abs_path, form.chain1bis_id.data, form.chain2bis_id.data,
                                                                                                   form.res_range1bis.data, form.res_range2bis.data)

                        data_bis = {
                            "pitch": f"{pitch_bis:3.2f}",
                            "nb_monomers": f"{monomers_per_turn_bis:3.2f}",
                            "direction": direction_bis,
                            "rotation_angle": f"{hp_bis.angle:3.2f}",
                            "rotation_angle_deg": f"{math.degrees(hp_bis.angle):3.2f}",
                            "axis_point": [f"{coord:3.2f}" for coord in hp_bis.point],
                            "translation": f"{hp_bis.normtranslation:3.2f}",
                            "dmin"      : f"{dminbis:3.2f}",
                            "dmax"       : f"{dmaxbis:3.2f}"
                        }
                    else:
                        data_bis = None

                    return render_template('run.html', form=form, screw_data=screw_data, screw_data_bis=data_bis)
                except Exception as e:
                    return render_template('error.html',summary=str(e), error_log=traceback.format_exc())

        # Construction requested
        elif request.form.get("action") == "construct" and form.validate():

            #Generate UUID for page results
            uniq_id = uuid.uuid4().hex
            # Create result folder
            result_path = pathlib.Path(current_app.config['DATA_UPLOADS'], uniq_id)
            result_path.mkdir(exist_ok=True)

            pdb_filename = validate_input_structure(form.input_file, form.pdb_id, result_path)
            if pdb_filename:

                # Save to the database the form inputs
                # Only way for now to pass the form data to another page.
                # We could use session or flash messages but neither seems to fit the need.
                user_inputs = UserInputs(request_id=uniq_id,
                                         pdb_filename=pdb_filename,
                                         chain1_id=form.chain1_id.data,
                                         chain2_id=form.chain2_id.data,
                                         res_range1=form.res_range1.data,
                                         res_range2=form.res_range2.data,
                                         n_mer=form.n_mer.data, z_align=form.z_align.data)

                # Does the user input for a 2nd assembly ?
                if form.has_2nd_oligomer():
                    pdb_filename_2nd = validate_input_structure(form.input_file_2nd, form.pdb_id_2nd, result_path)
                    user_inputs.add_2nd_oligomer(pdb_filename_2nd=pdb_filename_2nd,
                                                 chain1bis_id=form.chain1bis_id.data,
                                                 chain2bis_id=form.chain2bis_id.data,
                                                 res_range1bis=form.res_range1bis.data,
                                                 res_range2bis=form.res_range2bis.data)

                db.session.add(user_inputs)
                db.session.commit()

                #Redirect to the results page with correct id
                return redirect(url_for('heligeom_bp.results', results_id=uniq_id))
    return render_template('run.html', form=form)


@heligeom_bp.route("/results/<results_id>", methods=['GET', 'POST'])
def results(results_id):

    try:
        query_result = db.session.query(UserInputs).filter(UserInputs.request_id == results_id).first()

        # Query the database to retrieve the form inputs
        query_result = db.session.query(UserInputs).filter(UserInputs.request_id == results_id).first()
        pdb_filename = query_result.pdb_filename
        chain1_id    = query_result.chain1_id
        chain2_id    = query_result.chain2_id
        res_range1   = query_result.res_range1
        res_range2   = query_result.res_range2
        n_mer        = query_result.n_mer
        z_align      = query_result.z_align

        #Ptools part
        # Name of the constructed PDB (used also in download function)
        pdb_out_name = f"Construct_N{n_mer}_{pdb_filename}"
        pdb_out_abs_path = pathlib.Path(current_app.config['DATA_UPLOADS'], results_id, pdb_out_name)

        pdb_abs_path = pathlib.Path(current_app.config['DATA_UPLOADS'], results_id, pdb_filename)
        #Run the Heligeom calculations and write the PDB result in pdb_out_abs_path
        hp, pitch, nb_monomers, direction, dmin, dmax = construct(pdb_abs_path, chain1_id, chain2_id, res_range1, res_range2,
                                                                  n_mer, pdb_out_abs_path)



        # Create dict of data to pass to render_template
        data = {
            "results_id": results_id,
            "pdb_input": pdb_filename,
            "chain1_id": chain1_id,
            "chain2_id": chain2_id,
            "res_range1": res_range1,
            "res_range2": res_range2,
            "pdb_out_name": pdb_out_name,
            "n_mer": n_mer,
            "z_align": z_align,
            "pitch": f"{pitch:3.2f}",
            "nb_monomers": f"{nb_monomers:3.2f}",
            "direction": direction,
            "rotation_angle": f"{hp.angle:3.2f}",
            "rotation_angle_deg": f"{math.degrees(hp.angle):3.2f}",
            "axis_point": [f"{coord:3.2f}" for coord in hp.point],
            "translation": f"{hp.normtranslation:3.2f}",
            "dmin"      : f"{dmin:3.2f}",
            "dmax"       : f"{dmax:3.2f}"
        }

        if not query_result.second_oligomer:
            return render_template('results.html',data=data)
        else:
            pdb_filename2 = query_result.pdb_filename_2nd
            if not pdb_filename2:
                pdb_filename2 = pdb_filename
            chain1bis_id    = query_result.chain1bis_id
            chain2bis_id    = query_result.chain2bis_id
            res_range1bis   = query_result.res_range1bis
            res_range2bis   = query_result.res_range2bis

            #Ptools part
            # Name of the constructed PDB (used also in download function)
            pdb_out_name2 = f"Construct_2nd_N{n_mer}_{pdb_filename2}"
            pdb_out_abs_path2 = pathlib.Path(current_app.config['DATA_UPLOADS'], results_id, pdb_out_name2)

            pdb_abs_path2 = pathlib.Path(current_app.config['DATA_UPLOADS'], results_id, pdb_filename2)
            #Run the Heligeom calculations and write the PDB result in pdb_out_abs_path
            hp2, pitch2, nb_monomers2, direction2, dminbis, dmaxbis = construct(pdb_abs_path2, chain1bis_id, chain2bis_id,
                                                                                res_range1bis, res_range2bis,
                                                                                n_mer, pdb_out_abs_path2)

            # compute FNAT
            fnat = analyze_fnat(pdb_abs_path, chain1_id, chain2_id, res_range1, res_range2,
                                pdb_abs_path2, chain1bis_id, chain2bis_id, res_range1bis, res_range2bis)

            # Create dict of data to pass to render_template
            data_bis = {
                "pdb_input": pdb_filename2,
                "chain1_id": chain1bis_id,
                "chain2_id": chain2bis_id,
                "res_range1": res_range1bis,
                "res_range2": res_range2bis,
                "pdb_out_name": pdb_out_name2,
                "n_mer": n_mer,
                "z_align": z_align,
                "pitch": f"{pitch2:3.2f}",
                "nb_monomers": f"{nb_monomers2:3.2f}",
                "direction": direction2,
                "rotation_angle": f"{hp2.angle:3.2f}",
                "rotation_angle_deg": f"{math.degrees(hp2.angle):3.2f}",
                "axis_point": [f"{coord:3.2f}" for coord in hp2.point],
                "translation": f"{hp2.normtranslation:3.2f}",
                "dmin"      : f"{dminbis:3.2f}",
                "dmax"       : f"{dmaxbis:3.2f}",
                "fnat"       : f"{fnat:3.4f}"
            }

            return render_template('results_2_oligomers.html',data=data, data_bis=data_bis)

    except Exception:
        return render_template('error.html',error_log=traceback.format_exc())


# Create a route for the generated PDB. Use to download it and in the LiteMol plugin
@heligeom_bp.route('/results/<results_id>"/<path:filename>',  methods=['GET', 'POST'])
def download(results_id, filename):
    return send_from_directory(pathlib.Path(current_app.config["DATA_UPLOADS"], results_id), filename, as_attachment=True)


@heligeom_bp.route('/help')
def help():
    return render_template('help.html')

@heligeom_bp.route('/contact')
def contact():
    return render_template('contact.html')

@heligeom_bp.route('/error')
def error():
    return render_template('error.html')
