"""Flask config."""
import os

basedir = os.path.abspath(os.path.dirname(__file__))


class Config:
    """Base config."""
    SECRET_KEY = os.environ.get('SECRET_KEY')
    STATIC_FOLDER = 'static'
    DATA_UPLOADS =  os.path.join(basedir, "heligeom/static/data/")
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    MAX_CONTENT_LENGTH = 10000 * 10000
    PDB_SERVER = "https://files.rcsb.org/download/"



class DevConfig(Config):
    SECRET_KEY = "dev-key"
    FLASK_ENV = 'development'
    DEBUG = True
    TESTING = False
    SQLALCHEMY_DATABASE_URI = 'sqlite:///' + os.path.join(basedir, 'heligeom.db')

class TestConfig(Config):
    SECRET_KEY = "pouet"
    FLASK_ENV = 'test'
    DEBUG = True
    TESTING = True
    WTF_CSRF_ENABLED = False
    SQLALCHEMY_TRACK_MODIFICATIONS = True
    SQLALCHEMY_DATABASE_URI = 'sqlite:///test.db'
