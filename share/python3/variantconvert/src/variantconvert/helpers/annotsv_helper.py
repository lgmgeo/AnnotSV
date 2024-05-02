"""
Helper functions are functions that can be called from config files.

See the docstring of helper_functions.HelperFunctions() for more info.
"""

import logging as log
import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from commons import get_genome
from variantconvert.helpers.helper_functions import HelperFunctions


def remove_decimal(func):
    """
    pandas considers NaN a float. Hence, all integer columns containing empty cells are cast as float and a '.0' is added to every value.
    This decorator removes them.
    """

    def wrapper(*args, **kwargs):
        result = func(*args, **kwargs)
        if result == None:
            log.warning(f"Helper with params '{args}' returned None (should be '.')")
            return result
        if result.endswith(".0"):
            return result[:-2]
        return result

    return wrapper


class AnnotSvHelper(HelperFunctions):
    def __init__(self, *args, **kwargs):
        super(AnnotSvHelper, self).__init__(*args, **kwargs)
        self.FULL = "full"
        self.SPLIT = "split"
        annotsv_functions = {
            "get_pos_annotsv": self.get_pos_annotsv,
            "get_ref_annotsv": self.get_ref_annotsv,
            "get_alt_for_bed_based_annotsv": self.get_alt_for_bed_based_annotsv,
            "get_end_annotsv": self.get_end_annotsv,
            "get_svlen_annotsv": self.get_svlen_annotsv,
        }
        self.dispatcher = {**self.dispatcher, **annotsv_functions}

    def _check_annotation_mode(self, mode):
        if mode not in (self.FULL, self.SPLIT):
            raise ValueError(
                f"All values in {self.config['VCF_COLUMNS']['INFO']['Annotation_mode']} should be either {self.FULL} or {self.SPLIT}"
            )

    def get_pos_annotsv(self, sv_start, tx_start, annotation_mode):
        """Unsure if this is a good idea. For now I'll return the original value and let users use tx_start by themselves"""
        return sv_start
        # self._check_annotation_mode(annotation_mode)
        # # tx_start won't cast as int directly. Compare as floats to do 1 less cast.
        # if annotation_mode == self.SPLIT and float(tx_start) > float(sv_start):
        #     return str(int(float(tx_start)))
        # return sv_start

    def get_ref_annotsv(self, chrom, start):
        f = get_genome(self.config["GENOME"]["path"])
        if self.config["GENOME"]["vcf_header"][0].startswith(
            "##contig=<ID=chr"
        ) and not chrom.startswith("chr"):
            chrom = "chr" + str(chrom)
        return f[chrom][int(start) - 1].seq.upper()

    @staticmethod
    def get_alt_for_bed_based_annotsv(sv_type):
        return "<" + sv_type + ">"

    @remove_decimal
    def get_end_annotsv(self, annotation_mode, sv_end, tx_end):
        """Unsure if this is a good idea. For now I'll return the original value and let users use tx_end by themselves"""
        return sv_end
        # if annotation_mode == self.FULL:
        #     return sv_end
        # elif annotation_mode == self.SPLIT:
        #     if sv_end < tx_end:
        #         return sv_end
        #     return tx_end

    @remove_decimal
    def get_svlen_annotsv(self, annotation_mode, sv_length, overlapped_tx_length):
        """Unsure if this is a good idea. For now I'll return the original value and let users use overlapped_tx_length by themselves"""
        return sv_length
        # if annotation_mode == self.FULL:
        #     return sv_length
        # elif annotation_mode == self.SPLIT:
        #     if sv_length == ".":
        #         return sv_length
        #     elif float(sv_length) < 0:
        #         return str(-int(float(overlapped_tx_length)))
        #     return overlapped_tx_length
