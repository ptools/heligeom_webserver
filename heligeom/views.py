import math
import pathlib
import shutil
import time
import traceback
import uuid

from flask import (
    Blueprint,
    abort,
    current_app,
    g,
    redirect,
    render_template,
    request,
    send_from_directory,
    url_for,
)
from ptools.io.exceptions import InvalidPDBFormatError

from .forms import Construction, InputStructures, validate_input_structure
from .models import UserInputs, db
from .tool import HeligeomInterface, MonomersDifferentSizeError, MonomerSizeZeroError

# Blueprint Configuration
heligeom_bp = Blueprint(
    "heligeom_bp", __name__, template_folder="templates", static_folder="static"
)

VERSION = "1.0.0"


# Functions used to time spent on each page
# TODO: Remove in production
@heligeom_bp.before_request
def log_route_start():
    g.start_time = time.time()


@heligeom_bp.after_request
def log_route_end(response):
    route = request.endpoint
    print(f"{route} ended after {time.time() - g.pop('start_time', None)}")
    return response


# Allow a global variable to passed by default to render_template()
@heligeom_bp.context_processor
def set_version():
    return dict(version=VERSION)


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
            # Handle cases when monomers have different sizes or one is equal to 0
            # For that, we need to create real RigidBodys
            try:
                HeligeomInterface(
                    result_path / pathlib.Path(pdb_filename),
                    form.chain1_id.data,
                    form.chain2_id.data,
                    form.res_range1.data,
                    form.res_range2.data,
                    form.core_filter1.data,
                    form.core_region1.data,
                    form.core_region2.data,
                )
            except InvalidPDBFormatError as e:
                modalError = {
                    "title": f"File Parsing Error for {pdb_filename}",
                    "message": str(e),
                }
                shutil.rmtree(result_path)  # Removed newly created folder due to the errors
                return render_template("run.html", form=form, modalError=modalError)
            except (MonomerSizeZeroError, MonomersDifferentSizeError) as e:
                modalError = {
                    "title": "Invalid Monomer Selection for the 1st Oligomeric form",
                    "message": str(e),
                }
                shutil.rmtree(result_path)  # Removed newly created folder due to the errors
                return render_template("run.html", form=form, modalError=modalError)

            # Save to the database the form inputs
            user_inputs = UserInputs(
                request_id=uniq_id,
                pdb_filename=pdb_filename,
                chain1_id=form.chain1_id.data,
                chain2_id=form.chain2_id.data,
                res_range1=form.res_range1.data,
                res_range2=form.res_range2.data,
                core_filter1=form.core_filter1.data,
                core_region1=form.core_region1.data,
                core_filter2=form.core_filter2.data,
                core_region2=form.core_region2.data,
            )

            # Does the user input for a 2nd assembly ?
            if form.has_2nd_oligomer():
                pdb_filename_2nd = validate_input_structure(
                    form.input_file_2nd, form.pdb_id_2nd, result_path
                )
                # If no new structure is provided, use the first one
                if pdb_filename_2nd is None:
                    pdb_filename_2nd = pdb_filename

                try:
                    HeligeomInterface(
                        result_path / pathlib.Path(pdb_filename_2nd),
                        form.chain1bis_id.data,
                        form.chain2bis_id.data,
                        form.res_range1bis.data,
                        form.res_range2bis.data,
                        form.core_filter1bis.data,
                        form.core_region1bis.data,
                        form.core_region2bis.data,
                    )
                except InvalidPDBFormatError as e:
                    modalError = {
                        "title": f"File Parsing Error for {pdb_filename_2nd}",
                        "message": str(e),
                    }
                    shutil.rmtree(result_path)  # Removed newly created folder due to the errors
                    return render_template("run.html", form=form, modalError=modalError)
                except (MonomerSizeZeroError, MonomersDifferentSizeError) as e:
                    modalError = {
                        "title": "Invalid Monomer Selection for the 2nd Oligomeric form",
                        "message": str(e),
                    }
                    shutil.rmtree(result_path)  # Removed newly created folder due to the errors
                    return render_template("run.html", form=form, modalError=modalError)

                user_inputs.add_2nd_oligomer(
                    pdb_filename_2nd=pdb_filename_2nd,
                    chain1bis_id=form.chain1bis_id.data,
                    chain2bis_id=form.chain2bis_id.data,
                    res_range1bis=form.res_range1bis.data,
                    res_range2bis=form.res_range2bis.data,
                    core_filter1bis=form.core_filter1bis.data,
                    core_region1bis=form.core_region1bis.data,
                    core_filter2bis=form.core_filter2bis.data,
                    core_region2bis=form.core_region2bis.data,
                )

            # Store the input values in the database
            db.session.add(user_inputs)
            db.session.commit()

            # Redirect to the results page with correct id
            return redirect(url_for("heligeom_bp.results", results_id=uniq_id))
    return render_template("run.html", form=form)


@heligeom_bp.route("/results/<results_id>", methods=["GET", "POST"])
def results(results_id):
    # Absolute path to the result folder
    path_to_result = pathlib.Path(current_app.config["DATA_UPLOADS"], results_id)

    # Create this folder not matter what.
    # It should be done after the form submission but results are also directly accessible with the URL.
    # Useful to store the error.log file
    path_to_result.mkdir(exist_ok=True)

    # Query the database to retrieve the form inputs
    query_result = UserInputs.query.filter_by(request_id=results_id).first()

    # Redirect to 404 page when the query didn't return anything
    if not query_result:
        abort(404)

    # Path for all output files accessible to the user
    pdb_monomers = "mono1-2.pdb"
    contacts_csv = "contacts_interface1.csv"

    try:
        # Extract data from the database
        pdb_filename = query_result.pdb_filename
        chain1_id = query_result.chain1_id
        chain2_id = query_result.chain2_id
        res_range1 = query_result.res_range1
        res_range2 = query_result.res_range2
        core_filter1 = query_result.core_filter1
        core_region1 = query_result.core_region1
        core_filter2 = query_result.core_filter2
        core_region2 = query_result.core_region2

        pdb_abs_path = path_to_result / pdb_filename

        # Object to handle all the computation
        heli_interface1 = HeligeomInterface(
            pdb_abs_path,
            chain1_id,
            chain2_id,
            res_range1,
            res_range2,
            core_filter1,
            core_region1,
            core_region2,
        )

        # Save monomers only for the 1st molstar visualisation, depicting the interface
        heli_interface1.save_monomers(path_to_result / pdb_monomers)

        # Compute Helicoidal parameters
        pitch, monomers_per_turn, direction, dmin, dmax, rmsd = heli_interface1.compute_screw()

        # Information passed to the jinja2 template
        screw_data = {
            "results_id": results_id,
            "pdb_input": pdb_filename,
            "monos_file": pdb_monomers,
            "chain1_id": chain1_id,
            "chain2_id": chain2_id,
            "res_range1": res_range1,
            "res_range2": res_range2,
            "core_filter1": core_filter1,
            "core_region1": core_region1,
            "core_filter2": core_filter2,
            "core_region2": core_region2,
            "pitch": f"{pitch:3.2f}",
            "nb_monomers": f"{monomers_per_turn:3.2f}",
            "direction": direction,
            "rotation_angle": f"{heli_interface1.hp.angle:3.2f}",
            "rotation_angle_deg": f"{math.degrees(heli_interface1.hp.angle):3.2f}",
            "axis_point": [f"{coord:3.2f}" for coord in heli_interface1.hp.point],
            "translation": f"{heli_interface1.hp.normtranslation:3.2f}",
            "dmin": f"{dmin:3.2f}",
            "dmax": f"{dmax:3.2f}",
            "rmsd": f"{rmsd:3.2f}",
            "select_interface_mono": heli_interface1.molstar_selection_monomers(
                path_to_result / contacts_csv
            ),
            "select_core1": heli_interface1.molstar_select_core_monomer1,
            "select_core2": heli_interface1.molstar_select_core_monomer2,
            "contacts_csv": contacts_csv,
        }

        # 2nd Oligomer
        has_2nd_oligo = False
        screw_data_bis = dict()
        heli_interface2 = None

        pdb_filename2 = query_result.pdb_filename_2nd

        if pdb_filename2:
            has_2nd_oligo = True
            chain1bis_id = query_result.chain1bis_id
            chain2bis_id = query_result.chain2bis_id
            res_range1bis = query_result.res_range1bis
            res_range2bis = query_result.res_range2bis
            core_filter1bis = query_result.core_filter1bis
            core_region1bis = query_result.core_region1bis
            core_filter2bis = query_result.core_filter2bis
            core_region2bis = query_result.core_region2bis

            # Different input structure ?
            if pdb_filename2 != pdb_filename:
                pdb_abs_path2 = path_to_result / pdb_filename2
            else:
                pdb_abs_path2 = pdb_abs_path

            heli_interface2 = HeligeomInterface(
                pdb_abs_path2,
                chain1bis_id,
                chain2bis_id,
                res_range1bis,
                res_range2bis,
                core_filter1bis,
                core_region1bis,
                core_region2bis,
            )

            heli_interface2.save_monomers(path_to_result / "mono1-2_bis.pdb")

            # Compute Helicoidal parameters
            pitch2, monomers_per_turn2, direction2, dmin2, dmax2, rmsd2 = (
                heli_interface2.compute_screw()
            )

            contactsbis_csv = "contacts_interface2.csv"

            # Compute FNAT
            fnat = HeligeomInterface.compute_fnat(heli_interface1, heli_interface2)

            # Create dict of data to pass to render_template
            screw_data_bis = {
                "results_id": results_id,
                "pdb_input": pdb_filename2,
                "monos_file": "mono1-2_bis.pdb",
                "chain1_id": chain1bis_id,
                "chain2_id": chain2bis_id,
                "res_range1": res_range1bis,
                "res_range2": res_range2bis,
                "core_filter1": core_filter1bis,
                "core_region1": core_region1bis,
                "core_filter2": core_filter2bis,
                "core_region2": core_region2bis,
                "pitch": f"{pitch2:3.2f}",
                "nb_monomers": f"{monomers_per_turn2:3.2f}",
                "direction": direction2,
                "rotation_angle": f"{heli_interface2.hp.angle:3.2f}",
                "rotation_angle_deg": f"{math.degrees(heli_interface2.hp.angle):3.2f}",
                "axis_point": [f"{coord:3.2f}" for coord in heli_interface2.hp.point],
                "translation": f"{heli_interface2.hp.normtranslation:3.2f}",
                "dmin": f"{dmin2:3.2f}",
                "dmax": f"{dmax2:3.2f}",
                "rmsd": f"{rmsd2:3.2f}",
                "fnat": f"{fnat:3.4f}",
                "select_interface_mono": heli_interface2.molstar_selection_monomers(
                    path_to_result / contactsbis_csv
                ),
                "select_core1": heli_interface2.molstar_select_core_monomer1,
                "select_core2": heli_interface2.molstar_select_core_monomer2,
                "contacts_csv": contactsbis_csv,
            }
        # End if 2nd Oligomer

        # Initialize construction form
        construct_form = Construction()
        if request.method == "POST" and construct_form.validate():
            # Retrieve form values
            n_mer = construct_form.n_mer.data
            z_align = construct_form.z_align.data
            flatten = construct_form.flatten.data

            # Name of the constructed PDB (used also in download function)
            filename = pathlib.Path(pdb_filename).stem
            pdb_out_name = f'Construct_{filename}_N{n_mer}{"_Z" if z_align else ""}{"_flatten" if flatten else ""}.pdb'  # noqa: E501

            pdb_out_abs_path = path_to_result / pdb_out_name

            if flatten:
                # Try to construct a ring oligomer and write the PDB result in pdb_out_abs_path
                ring_hp, ring_rmsd, ring_fnat = heli_interface1.heli_ring(
                    pdb_out_abs_path, n_mer, z_align
                )
                # Retrieve N & Pitch
                ring_angle_degrees = math.degrees(ring_hp.angle)
                ring_pitch = abs(ring_hp.normtranslation * (360.0 / abs(ring_angle_degrees)))
                ring_monomers_per_turn = 360.0 / abs(ring_angle_degrees)
                ring_data = {
                    "ring_angle": f"{ring_angle_degrees:3.2f}",
                    "ring_pitch": f"{ring_pitch:3.2f}",
                    "ring_nb_monomers": f"{ring_monomers_per_turn:3.2f}",
                    "ring_rmsd": f"{ring_rmsd:3.2f}",
                    "ring_fnat": f"{ring_fnat:3.4f}",
                }
            else:
                # Construct oligomer and write the PDB result in pdb_out_abs_path
                heli_interface1.construct_oligomer(pdb_out_abs_path, n_mer, z_align)

            # Create dict of construction details to pass to render_template
            # In the oligomer structure, monomer 1 will be the chain A and monomer 1' will be the chain B.
            construct_data = {
                "pdb_out_name": pdb_out_name,
                "select_monomer1": f"struct_asym_id: 'A', color:'{heli_interface1.colors_monomer1[1]}'",
                "select_monomer2": f"struct_asym_id: 'B', color:'{heli_interface1.colors_monomer2[1]}'",
                "select_interface": heli_interface1.molstar_selection_interface_oligomer(),
                "flatten": flatten,
            }

            # Add ring data if selected
            if flatten:
                construct_data.update(ring_data)  # type: ignore

            # 2nd Oligomer
            if has_2nd_oligo and heli_interface2:
                # Different input structure ?
                if pdb_filename2 != pdb_filename:
                    pdb_abs_path2 = path_to_result / pdb_filename2
                else:
                    pdb_abs_path2 = pdb_abs_path

                # Name of the constructed PDB (used also in download function)
                filename2 = pathlib.Path(pdb_filename2).stem
                pdb_out_name2 = f'Construct_{filename2}_2nd_N{n_mer}{"_Z" if z_align else ""}{"_flatten" if flatten else ""}.pdb'  # noqa: E501

                pdb_out_abs_path2 = path_to_result / pdb_out_name2

                if flatten:
                    # Try to construct a ring oligomer and write the PDB result in pdb_out_abs_path
                    ring_hp, ring_rmsd, ring_fnat = heli_interface2.heli_ring(
                        pdb_out_abs_path2, n_mer, z_align
                    )
                    # Retrieve N & Pitch
                    ring_angle_degrees = math.degrees(ring_hp.angle)
                    ring_pitch = abs(ring_hp.normtranslation * (360.0 / abs(ring_angle_degrees)))
                    ring_monomers_per_turn = 360.0 / abs(ring_angle_degrees)
                    ring_data_bis = {
                        "ring_angle": f"{ring_angle_degrees:3.2f}",
                        "ring_pitch": f"{ring_pitch:3.2f}",
                        "ring_nb_monomers": f"{ring_monomers_per_turn:3.2f}",
                        "ring_rmsd": f"{ring_rmsd:3.2f}",
                        "ring_fnat": f"{ring_fnat:3.4f}",
                    }
                else:
                    # Construct oligomer and write the PDB result in pdb_out_abs_path
                    heli_interface2.construct_oligomer(pdb_out_abs_path2, n_mer, z_align)  # type: ignore

                # Create dict of construction details to pass to render_template
                construct_data_bis = {
                    "pdb_out_name": pdb_out_name2,
                    "select_monomer1": f"struct_asym_id: 'A', color:'{heli_interface2.colors_monomer1[1]}'",
                    "select_monomer2": f"struct_asym_id: 'B', color:'{heli_interface2.colors_monomer2[1]}'",
                    "select_interface": heli_interface2.molstar_selection_interface_oligomer(),
                    "flatten": flatten,
                }
                # Add ring data if selected
                if flatten:
                    construct_data_bis.update(ring_data_bis)  # type: ignore

                return render_template(
                    "results_2_oligomers.html",
                    construct_form=construct_form,
                    screw_data=screw_data,
                    screw_data_bis=screw_data_bis,
                    construct_data=construct_data,
                    construct_data_bis=construct_data_bis,
                )
            # End if 2nd Oligomer
            else:
                return render_template(
                    "results.html",
                    construct_form=construct_form,
                    screw_data=screw_data,
                    construct_data=construct_data,
                )

        if has_2nd_oligo:
            return render_template(
                "results_2_oligomers.html",
                construct_form=construct_form,
                screw_data=screw_data,
                screw_data_bis=screw_data_bis,
            )
        else:
            return render_template(
                "results.html", construct_form=construct_form, screw_data=screw_data
            )

    except Exception:
        # Save the error log into a file and redirect to a 500 error
        error_file = path_to_result / "error.log"
        with open(error_file, "w") as log:
            log.write(traceback.format_exc())
        # TODO: switch to "abort(500)" in prod
        return render_template("error.html", error_log=traceback.format_exc())


# Create a route for the generated PDB. Use to download it and in the Molstar plugin
@heligeom_bp.route("/results/<results_id>/<path:filename>", methods=["GET", "POST"])
def download(results_id, filename):
    return send_from_directory(
        pathlib.Path(current_app.config["DATA_UPLOADS"], results_id), filename, as_attachment=True
    )


@heligeom_bp.route("/help")
def help():
    return render_template("help.html")


@heligeom_bp.route("/about")
def about():
    return render_template("about.html")


@heligeom_bp.route("/contact")
def contact():
    return render_template("contact.html")


@heligeom_bp.route("/error")
def error():
    return render_template("error.html")


@heligeom_bp.app_errorhandler(404)
def page_not_found(e):
    return render_template("404.html"), 404


@heligeom_bp.app_errorhandler(500)
def internal_server_error(e):
    return render_template("500.html"), 500
