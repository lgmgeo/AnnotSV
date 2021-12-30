from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from annotsv.context import Context
from annotsv.util import ymd

OMIM_COLS = {
    "gene1": "Gene Symbols",
    "gene2": "Approved Gene Symbol",
    "mim": "MIM Number",
    "pheno": "Phenotypes",
}

PHENO_MAP = {
    "Autosomal dominant": "AD",
    "Autosomal recessive": "AR",
    "X-linked dominant": "XLD",
    "X-linked recessive": "XLR",
    "Y-linked dominant": "YLD",
    "Y-linked recessive": "YLR",
    "X-linked": "XL",
    "Y-linked": "YL",
}


def check_omim_file(app: Context):
    file_downloaded = app.config.omim_dir / "genemap2.txt"
    omim1_files = list(app.config.omim_dir.glob("*_OMIM-1-annotations.tsv.gz"))
    omim2_files = list(app.config.omim_dir.glob("*_OMIM-2-annotations.tsv.gz"))

    if not file_downloaded.exists():
        app.log.debug("No new OMIM annotation to format")
    elif len(omim1_files) > 1:
        app.keep_last_file("OMIM-1", omim1_files)
    elif len(omim2_files) > 1:
        app.keep_last_file("OMIM-2", omim2_files)
    elif not omim1_files or not omim2_files:
        ## - Create the 'date'_OMIM-1-annotations.tsv and 'date'_OMIM-2-annotations.tsv files.
        ##   Header1: genes, OMIM_ID
        ##   Header2: genes, OMIM_phenotype, OMIM_inheritance

        omim1_file = app.config.omim_dir / f"{ymd()}_OMIM-1-annotations.tsv.gz"
        omim2_file = app.config.omim_dir / f"{ymd()}_OMIM-2-annotations.tsv.gz"

        app.log.info(f"creating {omim1_file.name}, {omim2_file.name} in {app.config.omim_dir}")
        with file_downloaded.open("rt") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if not reader.fieldnames:
                app.abort(f"No header found in {file_downloaded}")
            for fname in OMIM_COLS.values():
                if fname not in reader.fieldnames:
                    app.abort(f"Missing required field {fname} in {file_downloaded}")

            for row in reader:
                all_genes = set(row[OMIM_COLS["gene1"]].split(","))
                all_genes.add(row[OMIM_COLS["gene2"]])

                pheno_tmp = row[OMIM_COLS["pheno"]]
                for pat, repl in PHENO_MAP.items():
                    pheno_tmp = re.sub(pat, repl, pheno_tmp, re.I)

                for p in pheno_tmp.split(";"):
                    match = re.search(r" *(.*?\(\d+\)),? *(.*)", p)
                    if match:
                        lpheno, linherit = match.groups()
                        if linherit:
                            linherit = re.sub(r" *, +", ",", linherit)
                            # TODO: finish later. file creation unnecessary if using downloaded annotation

        raise NotImplementedError()


def check_morbid_file(app: Context):
    file_downloaded = app.config.omim_dir / "morbidmap.txt"
    formatted_files = list(app.config.omim_dir.glob("*_morbid.tsv.gz"))
    candidate_files = list(app.config.omim_dir.glob("*_morbidCandidate.tsv.gz"))

    if not file_downloaded.exists() and not formatted_files:
        app.log.debug("No Morbid annotation")
    elif len(formatted_files) > 1:
        app.keep_last_file("OMIM Morbid Genes", formatted_files)
    elif len(candidate_files) > 1:
        app.keep_last_file("OMIM Morbid Genes candidate", candidate_files)
    elif not formatted_files:
        raise NotImplementedError()
    else:
        app.log.debug(f"No new Morbid annotation to format")


def omim2phenotype(app: Context, omim_id: str):
    ...


def omim_gene2phenotype(app: Context, omim_gene: str):
    ...


def is_morbid(app: Context, gene: str):
    ...
