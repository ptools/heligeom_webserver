from flask import Flask
from flask_cors import CORS
from flask_sqlalchemy import SQLAlchemy

db = SQLAlchemy()


def add_examples(db):
    """Manually add examples to the database

    Parameters
    ----------
    db : SQLAlchemy
        the database file
    """
    from .models import UserInputs

    gls = UserInputs(
        request_id="2GLS",
        pdb_filename="2GLS.pdb",
        chain1_id="A",
        chain2_id="B",
    )
    db.session.add(gls)
    db.session.commit()

    esv_core = UserInputs(
        request_id="4ESV_core",
        pdb_filename="4ESV.pdb",
        chain1_id="B",
        chain2_id="C",
    )
    esv_core.add_2nd_oligomer(
        pdb_filename_2nd="4ESV.pdb",
        chain1bis_id="B",
        chain2bis_id="C",
        core_filter1bis="manual",
        core_region1bis="184-365",
        core_filter2bis="manual",
        core_region2bis="184-365",
    )
    db.session.add(esv_core)
    db.session.commit()

    ice_AB_EF = UserInputs(
        request_id="3ICE_AB_EF",
        pdb_filename="3ICE.pdb",
        chain1_id="A",
        chain2_id="B",
    )
    ice_AB_EF.add_2nd_oligomer(pdb_filename_2nd="3ICE.pdb", chain1bis_id="E", chain2bis_id="F")
    db.session.add(ice_AB_EF)
    db.session.commit()


def create_app(config_object=None, clear_database=True):
    # create and configure the app
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_object(config_object)
    CORS(app)

    db.init_app(app)

    with app.app_context():
        # Clear database
        from .models import UserInputs

        db.create_all()
        # Clear database
        if clear_database:
            db.session.query(UserInputs).delete()
            db.session.commit()

            # Add examples to the DB
            add_examples(db)

        # Include Routes
        from . import views

        # Register the blueprint
        app.register_blueprint(views.heligeom_bp)

        return app
