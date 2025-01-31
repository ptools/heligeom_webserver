"""
Tests for the Heligeom simple pages
"""


def test_index(client):
    """
    GIVEN a Flask application
    WHEN the '/' page is requested (GET)
    THEN check the response is valid
    """
    response = client.get("/")
    assert response.status_code == 200
    assert b"Heligeom" in response.data
    assert b"Relating interfaces in 3D architectures" in response.data


def test_help(client):
    """
    GIVEN a Flask application
    WHEN the '/help' page is requested (GET)
    THEN check the response is valid
    """
    response = client.get("/help")
    assert response.status_code == 200
    assert b"Documentation" in response.data


def test_contact(client):
    """
    GIVEN a Flask application
    WHEN the '/contact' page is requested (GET)
    THEN check the response is valid
    """
    response = client.get("/contact")
    assert response.status_code == 200
    assert b"Created" in response.data
    assert b"Address" in response.data
    assert b"UPR 9080 CNRS" in response.data
    assert b"Citations" in response.data
