"""
Configuration of the tests.
Instantiate app and db variables through fixtures
"""

import pytest

from heligeom import create_app, db
from heligeom.models import User_Inputs


@pytest.fixture(scope='module')
def client():
    """
    Configure the application for testing
    """
    flask_app = create_app()
    flask_app.config['TESTING'] = True
    flask_app.config['WTF_CSRF_ENABLED'] = False
    flask_app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///test.db'

    with flask_app.test_client() as client:
        yield client

@pytest.fixture(scope='module')
def init_database():
    """
    Configure the database for testing
    """

    # Create the database and the database table
    db.create_all()

    # Insert user data
    data = User_Inputs("32300a0108a84a26afb3f23eb6bafba2","1kx2.pdb","1bta.pdb", 5)
    db.session.add(data)

    # Commit the changes for the users
    db.session.commit()

    yield db  # this is where the testing happens!

    db.drop_all()
