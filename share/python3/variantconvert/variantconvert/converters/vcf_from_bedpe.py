# -*- coding: utf-8 -*-

import logging as log
import os
import pandas as pd
import sys

from converters.vcf_from_breakpoints import VcfFromBreakpoints

sys.path.append("..")
from commons import create_vcf_header, is_helper_func, clean_string
from helper_functions import HelperFunctions
from variant import Variant


class VcfFromBedpe(VcfFromBreakpoints):
    """Made for the BEDPE format
    https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format
    Each input line will result in two VCF lines, one for each side of the breakpoint.
    """

    def get_input_columns(self):
        """identify if bed has a header
        return the given or predefined columns"""
        bedpe_columns = [
            "chrom1",
            "start1",
            "end1",
            "chrom2",
            "start2",
            "end2",
            "name",
            "score",
            "strand1",
            "strand2",
        ]
        minimum_header = "\t".join(bedpe_columns[0:6])

        column_names = None
        with open(self.filepath, "r") as f:
            for l in f:
                for i in range(self.config["GENERAL"]["skip_rows"]):
                    next(l)
                if l.startswith("#"):
                    if not l.startswith(minimum_header):
                        raise ValueError(
                            f"{self.filepath} does not start with expected header.\nExpected at least: {minimum_header}\nGot instead:{l}"
                        )
                elif "\t" in l:
                    total_cols = len(l.split("\t"))
                    if total_cols < len(minimum_header.split("\t")):
                        raise ValueError(
                            f"{self.filepath} needs a minimum of 6 columns: {minimum_header}"
                        )
                    column_names = bedpe_columns[0:total_cols]
                else:
                    raise ValueError(f"{self.filepath} has to be tab-separated")

                return column_names

    def _init_dataframe(self):
        self.df = pd.read_csv(
            self.filepath,
            skiprows=self.config["GENERAL"]["skip_rows"],
            sep="\t",
            low_memory=False,
            names=self.get_input_columns(),
        )
        self.df.reset_index(drop=True, inplace=True)
        self.df.fillna(".", inplace=True)

        log.debug(self.df)
        self.df[self.unique_variant_id_column] = self.df.apply(
            lambda row: self._get_unique_variant_id(row), axis=1
        )
        log.debug(self.df)

        if "strand1" not in self.df.columns:
            self.df["strand1"] = self.config["GENERAL"]["defaut_strand_1"]
        if "strand2" not in self.df.columns:
            self.df["strand2"] = self.config["GENERAL"]["defaut_strand_1"]
