==============================
Heligeom WebServer
==============================

This is the implementation of a Web Server for the tool **Heligeom**. The webserver address is `https://heligeom.galaxy.ibpc.fr/ <https://heligeom.galaxy.ibpc.fr/>`_

**Heligeom** is a tool, part of `Ptools <https://github.com/ptools/ptools>`_, for characterizing, manipulating and assembling structural units with a screw organization. Each structural unit, or monomer, may be an individual protein or a protein hetero-multimer.
Heligeom relies on the structures of monomer-monomer interfaces both for deriving the transformations and for constructing a filament of specified length.
For the latter it is thus complementary to other packages that apply known space group symmetries to obtain the structures of supra-assemblies [1].

The Web server is based on `Flask v2.2 <https://flask.palletsprojects.com/en/2.2.x/>`_ for the backend and `Materialize v1.1.0 <https://materializecss.github.io/materialize/>`_ for the frontend.


Requirements
------------

This Web server requires:

* Python >= 3.11
* Python modules Flask, flask-sqlalchemy, flask-wtf, scipy, numpy
* the new Ptools version written in 100% Python (`https://github.com/ptools/ptools`).


Installation
------------

Clone the depot and use pipenv to create a virtual environnement with all dependencies (thanks to Pipfile).

.. code-block:: bash

    $ pipenv install


Then, install Ptools (documentation to come). Once installed, just run the test webserver from Flask:

.. code-block:: bash

    $ pipenv run python run.py


A `requirements.txt` file is also provided for pip.

Citation
------------

To cite the Heligeom webserver, please refer to the following publications:

References
------------

.. [1] Santuz H, Laurent B, Robert CH, and Prévost C (2025) Heligeom: A web resource to generate, analyze, and visualize filament architectures based on pairwise association geometries of biological macromolecules. (submitted).

.. [2] Boyer B, Ezelin J, Poulain P, Saladin A, Zacharias M, Robert CH, and Prévost C (2015) An Integrative Approach to the Study of Filamentous Oligomeric Assemblies, with Application to RecA. PLOS ONE 10(3): e0116414.

