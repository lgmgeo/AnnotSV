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
    vcf_header: Optional[List[str]] = None
    bed_header: Optional[Path] = None
    genes_file: Optional[Path] = None

    def __init__(self, config: AnnotSVConfig, log: Logger) -> None:
        self.config = config
        self.log = log
        self.sv_ident = set()
        self.id_map = {}
        self.sv_lens = {}

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
