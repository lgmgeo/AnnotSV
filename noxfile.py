import os

import nox

os.environ.update(PDM_IGNORE_SAVED_PYTHON="1", PDM_USE_VENV="1")


@nox.session(python=("3.7", "3.8", "3.9", "3.10"))
def test(session):
    session.run("pdm", "install", "-Gtest", external=True)
    session.run("pytest", "tests/")
