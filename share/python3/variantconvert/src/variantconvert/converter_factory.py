# -*- coding: utf-8 -*-
"""
@Author: Samuel Nicaise
@Version: v2.0.0
"""

from variantconvert.converters.vcf_from_annotsv import VcfFromAnnotsv
from variantconvert.converters.vcf_from_bed import VcfFromBed
from variantconvert.converters.vcf_from_bedpe import VcfFromBedpe
from variantconvert.converters.vcf_from_breakpoints import VcfFromBreakpoints
from variantconvert.converters.vcf_from_tsv import VcfFromTsv
from variantconvert.converters.vcf_from_varank import VcfFromVarank
from variantconvert.converters.vcf_from_snp import VcfFromSnp

import logging as log


class ConverterFactory:
    """
    Factory pattern implementation
    To add new converters, use register_converter() or add them directly in __init__()
    """

    def __init__(self):
        self._converters = {}
        self._converters["varank>vcf"] = VcfFromVarank
        self._converters["annotsv>vcf"] = VcfFromAnnotsv
        self._converters["bed>vcf"] = VcfFromBed
        self._converters["tsv>vcf"] = VcfFromTsv
        self._converters["breakpoints>vcf"] = VcfFromBreakpoints
        self._converters["snp>vcf"] = VcfFromSnp
        self._converters["bedpe>vcf"] = VcfFromBedpe

    def register_converter(self, source_format, dest_format, converter):
        self._converters[source_format + ">" + dest_format] = converter

    def get_converter(self, source_format, dest_format, config):
        converter = self._converters.get(source_format + ">" + dest_format)
        if not converter:
            raise ValueError(
                f"Unknown converter: {source_format}>{dest_format} - probably a mistake in <yourconfig.yml> >> [GENERAL] >> input_format or output_format"
            )
        return converter(config)
