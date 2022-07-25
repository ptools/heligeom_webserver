===============================
Heligeom WebServer
===============================

This is the implementation of a Web Server for the tool **Heligeom**.

**Heligeom** is a tool, part of `Ptools <https://github.com/ptools/ptools>`_, to characterizing, manipulating and assembling structural units with a screw organization, and in which the structural units may be individual proteins or protein hetero-multimers.
Heligeom relies on the structures of monomer-monomer interfaces both for deriving the transformations and for filament construction.
For the latter it is thus complementary to other packages that apply known space group symmetries to obtain the structures of supra-assemblies [1].

The Web server is based on `Flask <https://flask.palletsprojects.com/en/1.1.x/>`_ for the backend and `Materialize <https://materializecss.com/>`_ for the frontend.


Requirements
------------

This Web server requires:

* Python >= 3.9
* Python modules: Flask, flask-sqlalchemy, flask-wtf
* the new Ptools version, rewritten in full Python.

Installation
------------

Clone the depot and use pipenv to create a virtual environnement with all dependencies (thanks to Pipfile).

.. code-block:: bash

    $ pipenv install


Then, install Ptools (documentation to come). Once installed, just run the test webserver from Flask:

.. code-block:: bash

    $ pipenv run python run.py




.. [1] Boyer B, Ezelin J, Poulain P, Saladin A, Zacharias M, et al. (2015) An Integrative Approach to the Study of Filamentous Oligomeric Assemblies, with Application to RecA. PLOS ONE 10(3): e0116414.

