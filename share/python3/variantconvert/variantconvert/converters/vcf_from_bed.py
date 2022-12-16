# -*- coding: utf-8 -*-
"""
"""
from __future__ import division
from __future__ import print_function

import logging as log
import os
import sys

from os.path import join as osj
import pandas as pd

from converters.abstract_converter import AbstractConverter

sys.path.append("..")
from helper_functions import HelperFunctions


class VcfFromBed:
    """
    TODO: insert related celine's code here
    """

    def convert(self, bed, output_path):
        log.info("Converting to vcf from bed using config: " + self.config_filepath)
        raise ValueError(
            "Not implemented yet. In most cases you should be able to use the TSV to VCF converter with an appropriate config. One exists for CANOES"
        )


if __name__ == "__main__":
    pass
