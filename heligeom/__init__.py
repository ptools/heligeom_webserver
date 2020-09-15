
import os
from flask import Flask
from flask_sqlalchemy import SQLAlchemy


basedir = os.path.abspath(os.path.dirname(__file__))

app = Flask(__name__)
app.config['SECRET_KEY'] = "you-will-never-guess"
app.config['IMAGE_UPLOADS'] =  os.path.join(basedir, "static/data/")
app.config['MAX_CONTENT_LENGTH'] = 10000 * 10000
# Database
app.config["SQLALCHEMY_DATABASE_URI"] = 'sqlite:///' + os.path.join(basedir, 'heligeom.db')
app.config["SQLALCHEMY_TRACK_MODIFICATIONS"] = False

app.config["PDB_SERVER"] = "https://files.rcsb.org/download/"

db = SQLAlchemy(app)


from . import views
