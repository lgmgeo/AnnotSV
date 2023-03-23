# -*- coding: utf-8 -*-

import logging as log
import natsort
import os
import pandas as pd
import sys

from converters.abstract_converter import AbstractConverter

sys.path.append("..")
from commons import create_vcf_header, is_helper_func, clean_string
from helper_functions import HelperFunctions
from variant import Variant


class VcfFromBreakpoints(AbstractConverter):
    """Made for file formats such as the TSV outputs of STAR-Fusion and ARRIBA.
    Other converters (for now) are not able to generate a VCF containing breakpoints.

    Each input line will result in two VCF lines, one for each side of the breakpoint.
    """

    def __init__(self, *args, **kwargs):
        super(VcfFromBreakpoints, self).__init__(*args, **kwargs)
        self.unique_variant_id_column = "__!UNIQUE_VARIANT_ID!__"

    def _init_dataframe(self):
        self.df = pd.read_csv(
            self.filepath,
            skiprows=self.config["GENERAL"]["skip_rows"],
            sep="\t",
            low_memory=False,
        )
        self.df.reset_index(drop=True, inplace=True)
        self.df.fillna(".", inplace=True)
        log.debug(self.df)
        self.df[self.unique_variant_id_column] = self.df.apply(
            lambda row: self._get_unique_variant_id(row), axis=1
        )
        log.debug(self.df)

    def _get_sample_list(self):
        # is the file multisample?
        if self.config["VCF_COLUMNS"]["SAMPLE"] != "":
            sample_list = []
            for sample in self.df[self.config["VCF_COLUMNS"]["SAMPLE"]].unique():
                sample_list.append(sample)
            return sample_list
        else:
            default_name = os.path.basename(self.output_path)
            if default_name.endswith(".vcf"):
                default_name = default_name[:-4]
            return [default_name]

    def _get_unique_variant_id(self, row):
        variant_id = []
        for col in self.config["GENERAL"]["unique_variant_id"]:
            variant_id.append(str(row[col]))
        return "_".join(variant_id)

    def _get_unique_id_to_index_list(self, data):
        id_dic = {}
        for k, v in data[self.unique_variant_id_column].items():
            if v not in id_dic:
                id_dic[v] = [k]
            else:
                id_dic[v].append(k)
        return id_dic

    def convert(self, tsv, output_path):
        log.debug("Converting to vcf from annotSV using config: " + self.config_filepath)

        self.filepath = tsv
        self.output_path = output_path
        self._init_dataframe()
        sample_list = self._get_sample_list()
        helper = HelperFunctions(self.config)

        with open(output_path, "w") as vcf:
            vcf_header = create_vcf_header(tsv, self.config, sample_list, True)
            for l in vcf_header:
                vcf.write(l + "\n")

            data = self.df.astype(str).to_dict()
            # In some variant callers, output files contain a list of variant-sample associations
            # so the same variant can be on multiple lines
            # __!UNIQUE_VARIANT_ID!__ allows to identify such variants and only add them to the VCF once
            already_seen_variants = set()
            unique_id_to_index_list = self._get_unique_id_to_index_list(data)

            final_variants_list = []
            for i in range(len(data[self.unique_variant_id_column])):
                if data[self.unique_variant_id_column][i] in already_seen_variants:
                    continue

                lines = [[], []]  # left side of the breakpoint, right side of the breakpoint
                left_var = Variant()
                right_var = Variant()

                for vcf_col in ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL"]:
                    col = self.config["VCF_COLUMNS"][vcf_col]

                    if is_helper_func(col):
                        # col[1] is a function name, col[2] its list of args
                        # the function named in col[1] has to be callable from this module
                        func = helper.get(col[1])
                        args = [data[c][i] for c in col[2:]]
                        result = func(*args)
                        if len(result) != 2:
                            raise ValueError(
                                "HELPER_FUNCTIONS used with vcf_from_breakpoints.py are expected to return a tuple of len 2. Got instead:"
                                + str(result)
                            )
                        left_var.set_column(vcf_col, result[0])
                        right_var.set_column(vcf_col, result[1])

                    elif col == "":
                        left_var.set_column(vcf_col, ".")
                        right_var.set_column(vcf_col, ".")
                    else:
                        left_var.set_column(vcf_col, data[col][i])
                        right_var.set_column(vcf_col, data[col][i])

                # doing this after defining chr/pos/ref/alt because they're required fields to make a unique hash
                if self.config["VCF_COLUMNS"]["ID"] == "":
                    left_var.id = left_var.get_hash()
                    right_var.id = right_var.get_hash()

                # Variant class not implementing INFO field yet: extract data into text lines now
                lines[0] = [
                    left_var.chrom,
                    left_var.pos,
                    left_var.get_hash(),
                    left_var.ref,
                    left_var.alt,
                    left_var.qual,
                ]
                lines[1] = [
                    right_var.chrom,
                    right_var.pos,
                    right_var.get_hash(),
                    right_var.ref,
                    right_var.alt,
                    right_var.qual,
                ]

                # Cutting-edge FILTER implementation
                lines[0].append("PASS")
                lines[1].append("PASS")

                left_info_field = ["SVTYPE=BND", "MATEID=" + right_var.get_hash()]
                right_info_field = ["SVTYPE=BND", "MATEID=" + left_var.get_hash()]
                for vcf_col, tsv_col in self.config["VCF_COLUMNS"]["INFO"].items():
                    if is_helper_func(tsv_col):
                        func = helper.get(tsv_col[1])
                        args = [data[c][i] for c in tsv_col[2:]]
                        s = vcf_col + "=" + func(*args)
                    else:
                        s = vcf_col + "=" + data[tsv_col][i]
                    left_info_field.append(clean_string(s))
                    right_info_field.append(clean_string(s))
                lines[0].append(";".join(left_info_field))
                lines[1].append(";".join(right_info_field))

                vcf_format_fields = []
                tsv_format_fields = []
                for vcf_col, tsv_col in self.config["VCF_COLUMNS"]["FORMAT"].items():
                    vcf_format_fields.append(vcf_col)
                    tsv_format_fields.append(tsv_col)
                lines[0].append(":".join(vcf_format_fields))
                lines[1].append(":".join(vcf_format_fields))

                # monosample input
                if len(sample_list) == 1:
                    sample_field = []
                    for key, val in self.config["VCF_COLUMNS"]["FORMAT"].items():
                        if key == "GT" and val == "":
                            sample_field.append(self.config["GENERAL"]["default_genotype"])
                            continue
                        sample_field.append(data[val][index])
                    lines[0].append(":".join(sample_field))
                    lines[1].append(":".join(sample_field))
                # multisample input
                else:
                    sample_field_dic = {}
                    # If the variant exists in other lines in the source file, fetch their sample data now
                    for index in unique_id_to_index_list[data[self.unique_variant_id_column][i]]:
                        sample_field = []
                        for key, val in self.config["VCF_COLUMNS"]["FORMAT"].items():
                            if key == "GT" and val == "":
                                sample_field.append(self.config["GENERAL"]["default_genotype"])
                                continue
                            sample_field.append(data[val][index])
                        sample_field_dic[
                            data[self.config["VCF_COLUMNS"]["SAMPLE"]][index]
                        ] = ":".join(sample_field)

                    for sample in sample_list:
                        if sample in sample_field_dic:
                            lines[0].append(sample_field_dic[sample])
                            lines[1].append(sample_field_dic[sample])
                        else:
                            if "GT" in self.config["VCF_COLUMNS"]["FORMAT"]:
                                if len(self.config["VCF_COLUMNS"]["FORMAT"]) == 1:
                                    # there's only GT. Avoid adding a trailing ":"
                                    empty = "./."
                                else:
                                    empty = "./.:" + ":".join(
                                        ["."] * (len(self.config["VCF_COLUMNS"]["FORMAT"]) - 1)
                                    )
                            else:
                                empty = ":".join(
                                    ["."] * (len(self.config["VCF_COLUMNS"]["FORMAT"]) - 1)
                                )
                            lines[0].append(empty)
                            lines[1].append(empty)
                    already_seen_variants.add(data[self.unique_variant_id_column][i])

                for line in lines:
                    final_variants_list.append(line)

            # finally, sort variants by chr/pos
            final_variants_list = natsort.natsorted(final_variants_list)
            for line in final_variants_list:
                vcf.write("\t".join(line) + "\n")
