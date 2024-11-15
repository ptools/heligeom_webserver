from flask import Flask
from flask_cors import CORS
from flask_sqlalchemy import SQLAlchemy

import flask_debugtoolbar

db = SQLAlchemy()


def create_app(config_object=None):
    # create and configure the app
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_object(config_object)
    CORS(app)

    db.init_app(app)

    # Specify the debug panels you want
    app.config['DEBUG_TB_PANELS'] = [
        'flask_debugtoolbar.panels.versions.VersionDebugPanel',
        'flask_debugtoolbar.panels.timer.TimerDebugPanel',
        'flask_debugtoolbar.panels.headers.HeaderDebugPanel',
        'flask_debugtoolbar.panels.request_vars.RequestVarsDebugPanel',
        'flask_debugtoolbar.panels.template.TemplateDebugPanel',
        'flask_debugtoolbar.panels.sqlalchemy.SQLAlchemyDebugPanel',
        'flask_debugtoolbar.panels.logger.LoggingPanel',
        'flask_debugtoolbar.panels.profiler.ProfilerDebugPanel',
        # Add the line profiling
        'flask_debugtoolbar_lineprofilerpanel.panels.LineProfilerPanel'
    ]
    toolbar = flask_debugtoolbar.DebugToolbarExtension(app)

    with app.app_context():
        # Clear database
        db.drop_all()
        from .models import UserInputs

        db.create_all()

        # Include Routes
        from . import views

        # Register the blueprint
        app.register_blueprint(views.heligeom_bp)

        return app
