from __future__ import annotations
from pathlib import Path

from logging import Logger
from typing import Dict, List, Optional, Set, Tuple

from annotsv.config import AnnotSVConfig


class Context:
    config: AnnotSVConfig
    log: Logger
    sv_ident: Set[str]
    id_map: Dict[Tuple[str, str, str], str]
    sv_lens: Dict[str, int]
    gccontent_ann: bool
    repeat_ann: bool
    segdup_ann: bool
    gap_ann: bool
    blacklist_ann: bool
    tad_ann: bool
    cosmic_ann: bool = False  # enabled in annotsv.cosmic
    mirna_ann: bool = False  # enabled in annotsv.regulartory_elements
    gh_ann: bool = False  # enabled in annotsv.regulartory_elements
    ea_ann: bool = False  # enabled in annotsv.regulartory_elements
    promoter_ann: bool = False  # enabled in annotsv.regulartory_elements
    cytoband_ann: bool = False  # enabled in annotsv.cytoband
    vcf_header: Optional[List[str]] = None
    bed_header: Optional[Path] = None
    genes_file: Optional[Path] = None
    refseq_genes: Set[str]
    ensembl_genes: Set[str]

    def __init__(self, config: AnnotSVConfig, log: Logger) -> None:
        self.config = config
        self.log = log
        self.sv_ident = set()
        self.id_map = {}
        self.sv_lens = {}
        self.refseq_genes = set()
        self.ensembl_genes = set()
        self.gccontent_ann = any(
            f"GC_content_{x}" in self.config.output_columns for x in ["left", "right"]
        )
        self.repeat_ann = any(
            f"Repeat_{x}_{y}" in self.config.output_columns
            for x in ["coord", "type"]
            for y in ["left", "right"]
        )
        self.segdup_ann = any(
            f"SegDup_{x}" in self.config.output_columns for x in ["left", "right"]
        )
        self.gap_ann = any(f"Gap_{x}" in self.config.output_columns for x in ["left", "right"])
        self.blacklist_ann = any(
            f"ENCODE_blacklist_{x}{y}" in self.config.output_columns
            for x in ["", "characteristics_"]
            for y in ["left", "right"]
        )
        self.tad_ann = any(
            x in self.config.output_columns for x in ["TAD_coordinate", "ENCODE_experiment"]
        )

    def abort(self, msg: str):
        self.log.error(msg)
        exit(2)

    def get_id(self, var_ident: str, ref: str, alt: str):
        svid = self.id_map.get((var_ident, ref, alt))
        if svid is None:
            svid = self.gen_id(var_ident)
            self.id_map[(var_ident, ref, alt)] = svid
        return svid

    def gen_id(self, var_ident: str):
        i = 1
        while True:
            svid = f"{var_ident}_{i}"
            if svid in self.sv_ident:
                i += 1
            else:
                break
        self.sv_ident.add(svid)
        return svid

    def keep_last_file(self, file_type: str, files: List[Path]):
        self.log.info(f"Found multiple {file_type} files, only keep last: {files}")
        for fpath in files[:-1]:
            fpath.rename(f"{fpath}.notused")
