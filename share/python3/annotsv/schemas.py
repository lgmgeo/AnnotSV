import re
from dataclasses import dataclass
from typing import Optional, Tuple

from cyvcf2 import Variant
from pydantic import BaseModel


def _strip_chr(val: str):
    return re.sub("chr", "", val, re.I)


@dataclass(init=False)
class Variant:
    __slots__ = ["chrom", "pos", "ref", "alt", "end", "svtype", "svlen", "cipos", "ciend"]
    chrom: str
    pos: int
    ref: str
    alt: str
    end: Optional[int]
    svtype: Optional[str]
    svlen: Optional[int]
    cipos: Optional[Tuple[int, int]]
    ciend: Optional[Tuple[int, int]]

    def __init__(self, var: Variant):
        self.chrom = _strip_chr(var.CHROM)
        # self.pos = var.POS
        self.pos = var.start + 1
        self.ref = _strip_chr(var.REF)
        if len(var.ALT) == 0:
            self.alt = ""
        else:
            self.alt = _strip_chr(",".join(var.ALT))
        self.svtype = var.INFO.get("SVTYPE")
        self.svlen = var.INFO.get("SVLEN")
        # cyvcf2 uses 0-indexed .end and .start
        self.end = var.INFO.get("END")
        self.cipos = var.INFO.get("CIPOS")
        self.ciend = var.INFO.get("CIEND")

    def __str__(self):
        return self.ident

    @property
    def ident(self):
        return f"{self.chrom}_{self.pos}_{self.end if self.end else ''}_{self.svtype if self.svtype else ''}"
