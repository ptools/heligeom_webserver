from flask import Flask
from flask_cors import CORS
from flask_sqlalchemy import SQLAlchemy

db = SQLAlchemy()


def create_app(config_object=None):
    # create and configure the app
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_object(config_object)
    CORS(app)

    db.init_app(app)

    with app.app_context():
        # Clear database
        db.drop_all()
        from .models import UserInputs

        db.create_all()

        # Add examples to the DB
        # gls = UserInputs(
        #     request_id="b83edb40aa934f76b866cd24e3aa",
        #     pdb_filename="2GLS.pdb",
        #     chain1_id="A",
        #     chain2_id="B",
        # )
        # db.session.add(gls)
        # db.session.commit()

        # Include Routes
        from . import views

        # Register the blueprint
        app.register_blueprint(views.heligeom_bp)

        return app
