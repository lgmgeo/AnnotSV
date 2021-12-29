from __future__ import annotations

from typing import Dict, List, Optional, Set, Tuple

from annotsv.config import AnnotSVConfig


class Context:
    config: AnnotSVConfig
    sv_ident: Set[str]
    id_map: Dict[Tuple[str, str, str], str]
    vcf_header: Optional[List[str]]
    sv_lens: Dict[str, int]

    def __init__(self, config: AnnotSVConfig) -> None:
        self.config = config
        self.sv_ident = set()
        self.id_map = {}
        self.vcf_header = None
        self.sv_lens = {}

    def get_id(self, sv_key: str, ref: str, alt: str):
        """sv_key: f"{chrom}_{pos}_{end}_{svtype}" """
        svid = self.id_map.get((sv_key, ref, alt))
        if svid is None:
            svid = self.gen_id(sv_key)
            self.id_map[(sv_key, ref, alt)] = svid
        return svid

    def gen_id(self, sv_key: str):
        i = 1
        while True:
            svid = f"{sv_key}_{i}"
            if svid in self.sv_ident:
                i += 1
            else:
                break
        self.sv_ident.add(svid)
        return svid
