# -*- coding: utf-8 -*-
"""
@Author: Samuel Nicaise
@Version: v1.2.1
"""

import hashlib


class Variant:
    def __init__(self):
        self.chrom = ""
        self.pos = ""
        self.id = ""
        self.ref = ""
        self.alt = ""
        self.qual = ""
        self.filter = ""
        self.info = {}

    def set_column(self, vcf_column, value):
        if vcf_column == "#CHROM":
            self.chrom = value
        elif vcf_column == "POS":
            self.pos = value
        elif vcf_column == "ID":
            self.id = value
        elif vcf_column == "REF":
            self.ref = value
        elif vcf_column == "ALT":
            self.alt = value
        elif vcf_column == "QUAL":
            self.qual = value
        elif vcf_column == "FILTER":
            self.filter = value
        elif vcf_column in ("INFO", "FORMAT"):
            raise NotImplementedError(
                "Variant.set_column(args) not implemented for INFO or FORMAT fields yet"
            )
        else:
            raise ValueError("Unexpected vcf_column: " + str(vcf_column))

    def set_hash(self, length=7):
        required_properties = [self.chrom, self.pos, self.ref, self.alt]
        if any([v == "" for v in required_properties]):
            raise RuntimeError(
                "CHROM, POS, REF, ALT are required for set_hash(). Current required_properties:"
                + ",".join([str(v) for v in required_properties])
            )
        hash_source = "".join([str(v) for v in required_properties])

        # based on the INFO keys reserved to encode structural variants
        optional_properties = ["SVTYPE", "LEN", "NOVEL", "IMPRECISE"]
        for prop in optional_properties:
            hash_source += self.info.get(prop, "")
        self.hash = hashlib.md5(hash_source.encode("utf-8")).hexdigest()[:length]

    def get_hash(self):
        if not hasattr(self, "hash"):
            self.set_hash()
        return self.hash
