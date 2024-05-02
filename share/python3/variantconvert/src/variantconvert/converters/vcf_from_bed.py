# -*- coding: utf-8 -*-
"""
"""

import logging as log
import os
import sys

from os.path import join as osj
import pandas as pd

from variantconvert.converters.abstract_converter import AbstractConverter
from variantconvert.helpers.helper_functions import HelperFunctions


class VcfFromBed:
    """
    TODO: insert related celine's code here
    """

    def __init__(*args, **kwargs):
        raise ValueError(
            "Not implemented yet. In most cases you should be able to use the TSV to VCF converter with an appropriate config. One exists for CANOES"
        )

    def convert(self, bed, output_path):
        log.debug("Converting to vcf from bed using config: " + self.config_filepath)
        raise ValueError(
            "Not implemented yet. In most cases you should be able to use the TSV to VCF converter with an appropriate config. One exists for CANOES"
        )
