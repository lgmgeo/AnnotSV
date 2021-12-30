from __future__ import annotations
from pathlib import Path

import re
from typing import List
from annotsv.constants import exomiser_dir
from annotsv.context import Context
from annotsv.enums import Organisms


def start_rest_service(app: Context, properties_file: Path, port: int, service_file: Path):
    ...


def check_exomiser_installation(app: Context):
    results_file = app.config.ncbi_dir / "results.txt"
    if app.config.organism is not Organisms.Human:
        app.log.debug(f"Disabling HPO for non-human organism {app.config.organism}")
        app.config.hpo = None
    elif not results_file.exists():
        app.log.warning(f"No Exomiser annotations available. {results_file} does not exist")
        app.config.hpo = None

    hpo_subdirs = list(app.config.exomiser_dir.glob("*"))
    for subd in hpo_subdirs[:]:
        if not re.match(r"\d+$", subd.name):
            hpo_subdirs.remove(subd)

    if hpo_subdirs:
        hpo_version = hpo_subdirs[-1].name
        app.log.info(
            "AnnotSV takes use of Exomiser (Smedley et al., 2015) for the phenotype-driven analysis."
        )
        app.log.info(
            f"AnnotSV is using the Human Phenotype Ontology (version {hpo_version}). "
            "Find out more at http://www.human-phenotype-ontology.org"
        )
    else:
        app.log.warning(f"No Exomiser annotations available in {app.config.exomiser_dir}")


def search_for_gene_id(app: Context, gene_name: str):
    ...


def run_exomiser(app: Context, gene_list: List[str], hpo_list: List[str]):
    ...


def exomiser_annotation(app: Context, gene_name: str, what: str):
    ...
