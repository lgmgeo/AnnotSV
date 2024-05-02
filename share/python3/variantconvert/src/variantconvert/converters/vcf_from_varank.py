# -*- coding: utf-8 -*-

import logging as log
import numpy
import os
import pandas as pd
import re
import sys
import time

from converters.abstract_converter import AbstractConverter

sys.path.append("..")
from commons import clean_string, is_helper_func, rename_duplicates_in_list, varank_to_vcf_coords
from variantconvert.helpers.helper_functions import HelperFunctions


class VcfFromVarank(AbstractConverter):
    """
    TODO: update vcffromvarank.py to fit the interface and import it instead
    """

    def _init_dataframe(self, filepath):
        self.filepath = filepath
        self.df = pd.read_csv(
            filepath,
            skiprows=self.config["GENERAL"]["skip_rows"],
            sep="\t",
            low_memory=False,
        )

        self.df.sort_values(
            [self.config["VCF_COLUMNS"]["#CHROM"], self.config["VCF_COLUMNS"]["POS"]],
            inplace=True,
        )
        self.df = self.df.drop_duplicates()  # varank files have duplicate lines!
        self.df.reset_index(drop=True, inplace=True)
        self.df.columns = rename_duplicates_in_list(self.df.columns)

        # convert french commas to dot in floats
        for field in ["INFO", "FORMAT"]:
            for col in self.config["COLUMNS_DESCRIPTION"][field]:
                if self.config["COLUMNS_DESCRIPTION"][field][col]["Type"] == "Float":
                    if col in self.df.columns:
                        self.df[col] = self.df.apply(
                            lambda row: self.french_commas_to_dots(row[col]), axis=1
                        )

        # request from Jean: remove the transcript part in cNomen columns
        if "cNomen" in self.df.columns:
            self.df["cNomen"] = self.df["cNomen"].apply(
                lambda row: self.remove_transcript_from_cnomen(row)
            )

        if "HI_percent" in self.df.columns:
            self.df["HI_percent"] = self.df["HI_percent"].apply(
                lambda row: self.remove_percent(row)
            )

        # homemade annotation
        self.add_gene_counts_to_df()
        log.debug(self.df)

    def remove_percent(self, val):
        if not isinstance(val, float):  # dirty way to check if value is not nan
            if val.endswith("%"):
                return val.split("%")[0]

    def remove_transcript_from_cnomen(self, val):
        if not isinstance(val, float):  # dirty way to check if value is not nan
            if ":" in val:
                return val.split(":")[1]

    def french_commas_to_dots(self, val):
        if isinstance(val, str):
            if "," in val:
                return val.replace(",", ".")
        return val

    def add_gene_counts_to_df(self):
        """
        Count how many times a gene appears muted in the sample and add it to the dataframe.
        This allows to create cutevariant filters selecting genes that have at least two separate heterozygote recessive mutations.

        No need to check for genotype because Varank TSV files are a list of variants contained in one sample.
        There are no "0/0" or "./." in the output VCF made from a Varank TSV file.
        """
        self.df["gene_mut_counts"] = self.df.groupby("genes")["genes"].transform("size")
        self.df["gene_mut_counts"] = self.df["gene_mut_counts"].fillna(-1)
        # pd.set_option('display.max_rows', None)
        # print(self.df["variantID"])
        # print(self.df["gene_mut_counts"])

    def get_sample_name(self, varank_tsv):
        # with open(varank_file, 'r') as f:
        # next(f)
        # name = f.readline()
        # if not name.startswith("## FamilyBarcode: "):
        # raise ValueError("Couldn't find FamilyBarcode in 2nd line of header. File causing issue: " + varank_file)
        # name = name.split("## FamilyBarcode: ")[1].strip()
        # return name
        name = os.path.basename(varank_tsv)
        name = re.sub("^fam[0-9]*_", "", name)
        for end in self.config["GENERAL"]["varank_filename_ends"]:
            if name.endswith(end):
                return name.split(end)[0]
        raise ValueError("Couldn't determine sample name from varank filename:" + varank_tsv)

    def set_coord_conversion_file(self, coord_conversion_file):
        self.coord_conversion_file = coord_conversion_file

    def get_known_columns(self):
        """
        TODO: load self.config[VCF_COLUMNS] instead and flatten it with https://stackoverflow.com/a/31439438
        """
        known = ["chr", "start", "end", "ref", "alt"]
        known.append("QUALphred")  # QUAL
        known.append("zygosity")  # GT
        known.append("totalReadDepth")  # DP
        known.append("varReadDepth")  # AD[1] (AD[0] computed through a HELPER_FUNCTION)
        known.append("varReadPercent")  # VAF
        known.append("gene_mut_counts")  # GMC ; see add_gene_counts_to_df(self)

        # Optionak STARK annotations not to be included
        known.append("FindByPipelines")
        known.append("BARCODE_2")  # Pool barcode
        known.append("POOL_F_base_counts")
        known.append("POOL_F_Depth")
        known.append("POOL_M_base_counts")
        known.append("POOL_M_Depth")
        known.append("trio_variant_type")
        # No GQ, no PL, and apparently no multi allelic variants
        return known

    def convert(self, varank_tsv, output_path):
        log.debug("Converting to vcf from varank using config: " + self.config_filepath)
        id_to_coords = varank_to_vcf_coords(self.coord_conversion_file)
        self.sample_name = self.get_sample_name(varank_tsv)
        self._init_dataframe(varank_tsv)
        helper = HelperFunctions(self.config)

        with open(output_path, "w") as vcf:
            vcf_header = self.create_vcf_header()
            for l in vcf_header:
                vcf.write(l + "\n")

            data = self.df.fillna(".").astype(str).to_dict()
            # "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", self.sample_name]
            for i in range(len(data["variantID"])):
                line = id_to_coords[data["variantID"][i]]["#CHROM"] + "\t"
                line += id_to_coords[data["variantID"][i]]["POS"] + "\t"
                line += data[self.config["VCF_COLUMNS"]["ID"]][i] + "\t"
                line += id_to_coords[data["variantID"][i]]["REF"] + "\t"
                line += id_to_coords[data["variantID"][i]]["ALT"] + "\t"
                line += data[self.config["VCF_COLUMNS"]["QUAL"]][i] + "\t"
                line += "PASS\t"

                # INFO
                info_field = []
                for key in data.keys():
                    if key not in self.get_known_columns():
                        s = key + "=" + data[key][i]
                        s = clean_string(s)
                        info_field.append(s)
                line += ";".join(info_field) + "\t"

                # FORMAT
                vcf_format_fields = []
                tsv_format_fields = []
                for vcf_col, tsv_col in self.config["VCF_COLUMNS"]["FORMAT"].items():
                    vcf_format_fields.append(vcf_col)
                    tsv_format_fields.append(tsv_col)
                line += ":".join(vcf_format_fields) + "\t"

                # Sample (Varank files are assumed to be monosample)
                sample_field = []
                for vcf_col, tsv_col in self.config["VCF_COLUMNS"]["FORMAT"].items():
                    if is_helper_func(tsv_col):
                        func = helper.get(tsv_col[1])
                        args = [data[c][i] for c in tsv_col[2:]]
                        sample_field.append(func(*args))
                    elif vcf_col == "GMC" and tsv_col == "":
                        sample_field.append(str(data["gene_mut_counts"][i]))
                    else:
                        sample_field.append(data[tsv_col][i])
                line += ":".join(sample_field)

                vcf.write(line + "\n")

        log.debug("Wrote: " + output_path)

    def create_vcf_header(self):
        header = []
        # basics
        header.append("##fileformat=VCFv4.3")
        header.append("##fileDate=%s" % time.strftime("%d/%m/%Y"))
        header.append("##source=" + self.config["GENERAL"]["origin"])
        header.append("##InputFile=%s" % os.path.abspath(self.filepath))

        # FILTER is not present in Varank, so all variants are set to PASS
        header.append('##FILTER=<ID=PASS,Description="Passed filter">')

        # INFO
        for key in self.config["COLUMNS_DESCRIPTION"]["INFO"]:
            if key in self.df.columns or key in self.config["VCF_COLUMNS"]:
                number = "1"  # TODO: update config to include number
                description = self.config["COLUMNS_DESCRIPTION"]["INFO"][key]["Description"]
                info_type = self.config["COLUMNS_DESCRIPTION"]["INFO"][key]["Type"]
                header.append(
                    f'##INFO=<ID={key},Number={number},Type={info_type},Description="{description}">'
                )

        # Undefined TSV fields also go to INFO field
        for key in self.df.columns:
            if (
                key not in self.get_known_columns()
                and key not in self.config["COLUMNS_DESCRIPTION"]["INFO"]
                and key not in self.config["COLUMNS_DESCRIPTION"]["FORMAT"]
            ):
                number = "1"  # TODO: length inference
                description = "Extracted from " + self.config["GENERAL"]["origin"]
                info_type = "String"  # No type inference is safer. Add known int/float annotations in config instead
                header.append(
                    f'##INFO=<ID={key},Number={number},Type={info_type},Description="{description}">'
                )

        # FORMAT
        for key in self.config["COLUMNS_DESCRIPTION"]["FORMAT"]:
            if key in self.df.columns or key in self.config["VCF_COLUMNS"]["FORMAT"]:
                number = self.config["COLUMNS_DESCRIPTION"]["FORMAT"][key]["Number"]
                description = self.config["COLUMNS_DESCRIPTION"]["FORMAT"][key]["Description"]
                format_type = self.config["COLUMNS_DESCRIPTION"]["FORMAT"][key]["Type"]
                header.append(
                    f'##FORMAT=<ID={key},Number={number},Type={format_type},Description="{description}">'
                )

        # GENOME
        header += self.config["GENOME"]["vcf_header"]
        header.append(
            "\t".join(
                [
                    "#CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "QUAL",
                    "FILTER",
                    "INFO",
                    "FORMAT",
                    self.sample_name,
                ]
            )
        )
        return header
