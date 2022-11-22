
import os
from flask_cors import CORS
from flask import Flask
from flask_sqlalchemy import SQLAlchemy

db = SQLAlchemy()

def create_app(config_object=None):

    # create and configure the app
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_object(config_object)
    cors = CORS(app)

    db.init_app(app)

    with app.app_context():

        #Clear database
        db.drop_all()
        from .models import UserInputs
        db.create_all()

        # Include Routes
        from . import views

        #Register the blueprint
        app.register_blueprint(views.heligeom_bp)

        return app
