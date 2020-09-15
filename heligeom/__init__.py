
import os
from flask import Flask
from flask_sqlalchemy import SQLAlchemy

db = SQLAlchemy()

def create_app(config_filename=None):

    basedir = os.path.abspath(os.path.dirname(__file__))

    # create and configure the app
    app = Flask(__name__, instance_relative_config=True)
    # app.config.from_mapping(
    #     SECRET_KEY='you-will-never-guess',
    #     DATABASE=os.path.join(app.instance_path, 'flaskr.sqlite'),
    # )

    app.config['SECRET_KEY'] = "you-will-never-guess"
    app.config['IMAGE_UPLOADS'] =  os.path.join(basedir, "static/data/")
    app.config['MAX_CONTENT_LENGTH'] = 10000 * 10000
    # Database
    app.config["SQLALCHEMY_DATABASE_URI"] = 'sqlite:///' + os.path.join(basedir, 'heligeom.db')
    app.config["SQLALCHEMY_TRACK_MODIFICATIONS"] = False

    app.config["PDB_SERVER"] = "https://files.rcsb.org/download/"

    # app.config.from_pyfile(config_filename, silent=True)

    db.init_app(app)

    with app.app_context():

        #Clear database
        db.drop_all(app=app)
        db.create_all(app=app)

        # Include Routes
        from . import views

        #Register the blueprint
        app.register_blueprint(views.heligeom_bp)

        return app
