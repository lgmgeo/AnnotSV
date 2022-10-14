# AnnotSv

[![Tests](https://github.com/lgmgeo/AnnotSV/workflows/Tests/badge.svg)](https://github.com/lgmgeo/AnnotSV/actions?query=workflow%3Aci)
[![pypi version](https://img.shields.io/pypi/v/AnnotSV.svg)](https://pypi.org/project/AnnotSV/)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![pdm-managed](https://img.shields.io/badge/pdm-managed-blueviolet)](https://pdm.fming.dev)

Annotation and Ranking of Structural Variation

## Requirements

AnnotSv requires Python >=3.7

## Installation

It is recommended to install with `pipx`, if `pipx` haven't been installed yet, refer to the [pipx's docs](https://github.com/pipxproject/pipx)

```bash
$ pipx install annotsv
```

Alternatively, install with `pip` to the user site:

```bash
$ python -m pip install --user annotsv
```

## Contribute

### Prepare

```
pdm install
```

### Run

```
pdm run AnnotSV
```

### Test

```
pdm run test
```

### Lint

```
pdm run black
pdm run isort
pdm run mypy
```

### Install from source

```
python -m pip install --user .
```
