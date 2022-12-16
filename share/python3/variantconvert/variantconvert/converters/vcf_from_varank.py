# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

import logging as log
import numpy
import os
import pandas as pd
import re
import sys
import time

from converters.abstract_converter import AbstractConverter

sys.path.append("..")
from commons import clean_string, rename_duplicates_in_list, varank_to_vcf_coords
from helper_functions import HelperFunctions


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
        for col in self.config["COLUMNS_DESCRIPTION"]:
            if self.config["COLUMNS_DESCRIPTION"][col]["Type"] == "Float":
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
        if not isinstance(val, float): #dirty way to check if value is not nan
            if val.endswith("%"):
                return val.split("%")[0]

    def remove_transcript_from_cnomen(self, val):
        if not isinstance(val, float): #dirty way to check if value is not nan
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
        raise ValueError(
            "Couldn't determine sample name from varank filename:" + varank_tsv
        )

    def set_coord_conversion_file(self, coord_conversion_file):
        self.coord_conversion_file = coord_conversion_file

    def get_known_columns(self):
        """
        TODO: load self.config[VCF_COLUMNS] instead and flatten it with https://stackoverflow.com/a/31439438
        """
        known = ["chr", "start", "end", "ref", "alt"]
        known.append("QUALphread")  # QUAL
        known.append("zygosity")  # GT
        known.append("totalRead")  # DP
        known.append("varReadDepth")  # AD[1] ; AD[0] = DP - AD[1]
        known.append("varReadPercent")  # VAF
        known.append("gene_mut_counts")  # GMC ; see add_gene_counts_to_df(self)
        # No GQ, no PL, and apparently no multi allelic variants
        return known

    def convert(self, varank_tsv, output_path):
        log.info("Converting to vcf from varank using config: " + self.config_filepath)
        id_to_coords = varank_to_vcf_coords(self.coord_conversion_file)
        self.sample_name = self.get_sample_name(varank_tsv)
        self._init_dataframe(varank_tsv)

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

                info_field = []
                for key in data.keys():
                    if key not in self.get_known_columns():
                        s = key + "=" + data[key][i]
                        s = clean_string(s)
                        info_field.append(s)
                line += ";".join(info_field) + "\t"

                line += "GT:DP:AD:VAF:GMC\t"
                gt_dic = {"hom": "1/1", "het": "0/1"}
                sample_field = (
                    gt_dic[data[self.config["VCF_COLUMNS"]["FORMAT"]["GT"]][i]] + ":"
                )
                sample_field += (
                    data[self.config["VCF_COLUMNS"]["FORMAT"]["DP"]][i] + ":"
                )
                sample_field += (
                    str(int(data["totalReadDepth"][i]) - int(data["varReadDepth"][i]))
                    + ","
                    + data["varReadDepth"][i]
                    + ":"
                )
                vaf = data[self.config["VCF_COLUMNS"]["FORMAT"]["VAF"]][i]
                if vaf != ".":
                    vaf = str(float(vaf) / 100)
                sample_field += vaf + ":"
                sample_field += str(data["gene_mut_counts"][i])
                line += sample_field

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

        # INFO contains all columns that are not used anywhere specific
        # log.debug(dict(self.df.dtypes))
        for key in self.df.columns:
            if key in self.get_known_columns():
                continue
            if str(self.df[key].dtypes) in ("object", "O", "bool"):
                info_type = "String"
            elif str(self.df[key].dtypes) == "float64":
                info_type = "Float"
            elif str(self.df[key].dtypes) in ("int64", "Int64"):
                info_type = "Integer"
            else:
                raise ValueError(
                    "Unrecognized type in Varank dataframe. Column causing issue: "
                    + key
                )

            if key in self.config["COLUMNS_DESCRIPTION"]:
                description = self.config["COLUMNS_DESCRIPTION"][key]["Description"]
                info_type = self.config["COLUMNS_DESCRIPTION"][key]["Type"]
            else:
                description = "Extracted from " + self.config["GENERAL"]["origin"]
                info_type = "String"  # ugly fix of that bug where bcftools change POOL_ columns to Float (--> cutevariant crash)
            header.append(
                "##INFO=<ID="
                + key
                + ",Number=1,Type="
                + info_type
                + ',Description="'
                + description
                + '">'
            )
        # FORMAT
        header.append(
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">'
        )
        header.append(
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">'
        )
        header.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        header.append(
            '##FORMAT=<ID=VAF,Number=1,Type=Float,Description="VAF Variant Frequency">'
        )
        header.append(
            '##FORMAT=<ID=GMC,Number=1,Type=String,Description="Gene Mutations count: number of variants occuring in the same gene based on <genes> column. Computed when Varank files are converted to VCF">'
        )
        # genome stuff
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


if __name__ == "__main__":
    pass
