"""
Configuration of the tests.
Instantiate app and db variables through fixtures
"""

import pytest
from heligeom import create_app, db
from heligeom.models import UserInputs


@pytest.fixture(scope="module")
def client():
    """
    Configure the application for testing
    """
    flask_app = create_app("config.TestConfig")

    with flask_app.test_client() as c:
        yield c


@pytest.fixture(scope="module")
def init_database():
    """
    Configure the database for testing
    """

    # Create the database and the database table
    db.create_all()

    # Insert user data
    data = UserInputs("32300a0108a84a26afb3f23eb6bafba2", "1kx2.pdb", "1bta.pdb")
    db.session.add(data)

    # Commit the changes for the users
    db.session.commit()

    yield db  # this is where the testing happens!

    db.drop_all()
