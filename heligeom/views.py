import math
import pathlib
import traceback
import uuid

from flask import (
    Blueprint,
    current_app,
    redirect,
    render_template,
    request,
    send_from_directory,
    url_for,
)

from .forms import Construction, InputStructures, validate_input_structure
from .models import UserInputs, db
from .tool import analyze_fnat, construct, screw_parameters

# Blueprint Configuration
heligeom_bp = Blueprint(
    "heligeom_bp", __name__, template_folder="templates", static_folder="static"
)


@heligeom_bp.route("/")
def homepage():
    return render_template("index.html")


@heligeom_bp.route("/run", methods=["GET", "POST"])
def runpage():
    # Initialize submission form
    form = InputStructures()

    if request.method == "POST" and form.validate():
        # Generate UUID for storing files
        uniq_id = uuid.uuid4().hex
        # Create result folder
        result_path = pathlib.Path(current_app.config["DATA_UPLOADS"], uniq_id)
        result_path.mkdir(exist_ok=True)

        pdb_filename = validate_input_structure(form.input_file, form.pdb_id, result_path)
        if pdb_filename:
            # Save to the database the form inputs
            # Only way for now to pass the form data to another page.
            # We could use session or flash messages but neither seems to fit the need.
            user_inputs = UserInputs(
                request_id=uniq_id,
                pdb_filename=pdb_filename,
                chain1_id=form.chain1_id.data,
                chain2_id=form.chain2_id.data,
                res_range1=form.res_range1.data,
                res_range2=form.res_range2.data,
            )

            # Does the user input for a 2nd assembly ?
            if form.has_2nd_oligomer():
                pdb_filename_2nd = validate_input_structure(
                    form.input_file_2nd, form.pdb_id_2nd, result_path
                )
                # If no new structure is provided, use the first one
                if pdb_filename_2nd is None:
                    pdb_filename_2nd = pdb_filename

                user_inputs.add_2nd_oligomer(
                    pdb_filename_2nd=pdb_filename_2nd,
                    chain1bis_id=form.chain1bis_id.data,
                    chain2bis_id=form.chain2bis_id.data,
                    res_range1bis=form.res_range1bis.data,
                    res_range2bis=form.res_range2bis.data,
                )

            db.session.add(user_inputs)
            db.session.commit()

            # Redirect to the results page with correct id
            return redirect(url_for("heligeom_bp.results", results_id=uniq_id))
    return render_template("run.html", form=form)


@heligeom_bp.route("/results/<results_id>", methods=["GET", "POST"])
def results(results_id):
    try:
        # Query the database to retrieve the form inputs
        query_result = UserInputs.query.filter_by(request_id=results_id).first()

        if not query_result:
            raise AttributeError("Could not find the results id.")

        pdb_filename = query_result.pdb_filename
        chain1_id = query_result.chain1_id
        chain2_id = query_result.chain2_id
        res_range1 = query_result.res_range1
        res_range2 = query_result.res_range2

        pdb_abs_path = pathlib.Path(current_app.config["DATA_UPLOADS"], results_id, pdb_filename)

        # Compute Helicoidal parameters
        hp, pitch, monomers_per_turn, direction, dmin, dmax = screw_parameters(
            pdb_abs_path,
            chain1_id,
            chain2_id,
            res_range1,
            res_range2,
        )

        screw_data = {
            "results_id": results_id,
            "pdb_input": pdb_filename,
            "chain1_id": chain1_id,
            "chain2_id": chain2_id,
            "res_range1": res_range1,
            "res_range2": res_range2,
            "pitch": f"{pitch:3.2f}",
            "nb_monomers": f"{monomers_per_turn:3.2f}",
            "direction": direction,
            "rotation_angle": f"{hp.angle:3.2f}",
            "rotation_angle_deg": f"{math.degrees(hp.angle):3.2f}",
            "axis_point": [f"{coord:3.2f}" for coord in hp.point],
            "translation": f"{hp.normtranslation:3.2f}",
            "dmin": f"{dmin:3.2f}",
            "dmax": f"{dmax:3.2f}",
        }

        # Initialize construction form
        construct_form = Construction()

        if construct_form.validate():
            # Retrieve form values
            n_mer = construct_form.n_mer.data
            z_align = construct_form.z_align

            # Ptools part
            # Name of the constructed PDB (used also in download function)
            if z_align:
                pdb_out_name = f"Construct_N{n_mer}_Z_{pdb_filename}"
            else:
                pdb_out_name = f"Construct_N{n_mer}_{pdb_filename}"

            pdb_out_abs_path = pathlib.Path(
                current_app.config["DATA_UPLOADS"], results_id, pdb_out_name
            )
            # Run the Heligeom calculations and write the PDB result in pdb_out_abs_path
            construct(
                pdb_abs_path,
                chain1_id,
                chain2_id,
                res_range1,
                res_range2,
                n_mer,
                z_align,
                pdb_out_abs_path,
            )

            # Create dict of construction details to pass to render_template
            construct_data = {
                "pdb_out_name": pdb_out_name,
                "n_mer": n_mer,
                "z_align": z_align,
            }

            return render_template(
                "results.html",
                construct_form=construct_form,
                screw_data=screw_data,
                construct_data=construct_data,
            )

        return render_template("results.html", construct_form=construct_form, screw_data=screw_data)

    except Exception:
        return render_template("error.html", error_log=traceback.format_exc())


# Create a route for the generated PDB. Use to download it and in the LiteMol plugin
@heligeom_bp.route('/results/<results_id>"/<path:filename>', methods=["GET", "POST"])
def download(results_id, filename):
    return send_from_directory(
        pathlib.Path(current_app.config["DATA_UPLOADS"], results_id), filename, as_attachment=True
    )


@heligeom_bp.route("/help")
def help():
    return render_template("help.html")


@heligeom_bp.route("/contact")
def contact():
    return render_template("contact.html")


@heligeom_bp.route("/error")
def error():
    return render_template("error.html")
