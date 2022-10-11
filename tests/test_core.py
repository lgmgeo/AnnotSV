import annotsv


def test_import_package():
    assert isinstance(annotsv.__all__, list)
