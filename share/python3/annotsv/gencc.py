from __future__ import annotations

import csv
import gzip
from pathlib import Path
import re
from dataclasses import dataclass
from typing import Dict, Set, TYPE_CHECKING
from annotsv.util import ymd
from annotsv.schemas import AnnotationValidator, ResolvedFiles

if TYPE_CHECKING:
    from annotsv.context import Context

GENCC_COLS = {
    "gene": "gene_symbol",
    "disease": "disease_title",
    "moi": "moi_title",
    "classif": "classification_title",
    "pmid": "submitted_as_pmids",
}


class GenCCValidator(AnnotationValidator):
    def __init__(self, app: Context):
        downloaded_rf = ResolvedFiles(app.config.gencc_dir, "submissions-export-tsv")
        formatted_rf = ResolvedFiles(app.config.gencc_dir, "*_GenCC.tsv", "*_GenCC.tsv.gz")
        super().__init__(
            app,
            label="GenCC",
            downloaded=downloaded_rf,
            formatted=formatted_rf,
        )

    ### GenCC
    #####################################
    # GenCC downloaded file: https://search.thegencc.org/download/action/submissions-export-tsv
    #
    # Header:
    # "uuid","gene_curie","gene_symbol","disease_curie","disease_title","disease_original_curie","disease_original_title","classification_curie","classification_title","moi_curie","moi_title","submitter_curie","submitter_title","submitted_as_hgnc_id","submitted_as_hgnc_symbol","submitted_as_disease_id","submitted_as_disease_name","submitted_as_moi_id","submitted_as_moi_name","submitted_as_submitter_id","submitted_as_submitter_name","submitted_as_classification_id","submitted_as_classification_name","submitted_as_date","submitted_as_public_report_url","submitted_as_notes","submitted_as_pmids","submitted_as_assertion_criteria_url","submitted_as_submission_id","submitted_run_date"
    ## - Check if the following GenCC file has been downloaded:
    #    - submissions-export-tsv
    #
    ## - Check and create if necessary the following file:
    #    - 'date'_GenCC.sorted.tsv.gz
    def update(self):
        downloaded_file = self.downloaded_path()
        formatted_file = self.formatted_path().with_name(
            self.formatted[0].patterns[1].replace("*", ymd())
        )

        self._app.log.info(
            f"GenCC configuration - creation of {formatted_file} (done only once during first GenCC annotation)"
        )

        # GenCC says it's tsv, but it's still just a csv
        with downloaded_file.open("rt") as csv_in:
            rdr = csv.DictReader(csv_in)
            if rdr.fieldnames is None:
                self._app.abort(f"No header found in {downloaded_file}")

            # stupid type checker
            assert rdr.fieldnames
            for fname in GENCC_COLS.values():
                if fname not in rdr.fieldnames:
                    self._app.abort(f"{fname} not found in {downloaded_file} header")

            genes: Dict[str, Gene] = {}
            for row in rdr:
                gene_name = row[GENCC_COLS["gene"]]

                if gene_name not in genes:
                    genes[gene_name] = Gene(row)
                else:
                    genes[gene_name].update(row)

        with gzip.open(formatted_file, "wt") as tsv_out:
            tsv_out.write("genes\tGenCC_disease\tGenCC_moi\tGenCC_classification\tGenCC_pmid\n")
            for g in sorted(genes):
                tsv_out.write(f"{genes[g]}\n")

        downloaded_file.unlink()


@dataclass
class Gene:
    __slots__ = ["name", "disease", "moi", "classif", "pmid"]
    name: str
    disease: Set[str]
    moi: Set[str]
    classif: Set[str]
    pmid: Set[str]

    def __init__(self, row: Dict[str, str]) -> None:
        self.name = row[GENCC_COLS["gene"]]
        self.disease = set()
        self.moi = set()
        self.classif = set()
        self.pmid = set()
        self.update(row)

    def __str__(self):
        d_str = ";".join(sorted(self.disease))
        m_str = ";".join(sorted(self.moi))
        c_str = ";".join(sorted(self.classif))
        p_str = ";".join(sorted(self.pmid))
        return "\t".join([self.name, d_str, m_str, c_str, p_str])

    def update(self, row: Dict[str, str]):
        if row[GENCC_COLS["disease"]]:
            self.disease.add(row[GENCC_COLS["disease"]])

        if row[GENCC_COLS["moi"]]:
            moi = row[GENCC_COLS["moi"]].lower()
            new_moi = None
            if "autosomal dominant" in moi:
                if "with maternal imprinting" in moi:
                    new_moi = "ADm"
                elif "with paternal imprinting" in moi:
                    new_moi = "ADm"
                else:
                    new_moi = "AD"
            elif "autosomal recessive" in moi:
                new_moi = "AR"
            elif "digenic" in moi:
                new_moi = "2G"
            elif "mitochondrial" in moi:
                new_moi = "MT"
            elif "semidominant" in moi:
                new_moi = "sD"
            elif "somatic mosaicism" in moi:
                new_moi = "SOM"
            elif "x-linked" in moi:
                if "dominant" in moi:
                    new_moi = "XLD"
                elif "recessive" in moi:
                    new_moi = "XLR"
                else:
                    new_moi = "XL"
            elif "y-linked" in moi:
                if "dominant" in moi:
                    new_moi = "YLD"
                elif "recessive" in moi:
                    new_moi = "YLR"
                else:
                    new_moi = "YL"

            if new_moi:
                self.moi.add(new_moi)
                line = ",".join(row.values()).lower()
                if "incomplete penetrance" in line or "reduced penetrance" in line:
                    self.moi.add("IPVE")

        if row[GENCC_COLS["classif"]]:
            self.classif.add(row[GENCC_COLS["classif"]])

        pmid = row[GENCC_COLS["pmid"]]
        if pmid:
            pmid_list = re.split(r"[,;]", pmid.replace(" ", ""))
            self.pmid.update([p for p in pmid_list if p not in ["", "0"]])
