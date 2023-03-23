# -*- coding: utf-8 -*-
"""
@Author: Samuel Nicaise
@Version: v1.2.1
"""

from commons import get_genome


class HelperFunctions:
    """
    For when you can't just convert columns by changing column names
    Steps needed:
    - define the helper function
    - update self.dispatcher so this class can redirect to the proper function
    - tell the config file you want to use a HELPER_FUNCTION with the following pattern:
            [HELPER_FUNCTION, <name of your function>, <arg1>, <arg2>, ..., <argN>]

    Example: I need a LENGTH value in my destination format,
    but my source file only has START and END columns.
    You would need:
    # somewhere in the module
            def get_length(start, end):
                    return str(end - start)
    #in this class __init__():
            self.dispatcher["get_length_from_special_format"]: get_length
    # in the JSON configfile
            LENGTH: [HELPER_FUNCTION, "get_length_from_special_format", START, END]
    """

    def __init__(self, config):
        self.config = config
        self.dispatcher = {
            "get_ref_from_decon": self.get_ref_from_decon,
            "get_alt_from_decon": self.get_alt_from_decon,
            "get_svlen_from_decon": self.get_svlen_from_decon,
            "get_info_from_annotsv": self.get_info_from_annotsv,
            "get_alt_for_bed_based_annotsv": self.get_alt_for_bed_based_annotsv,
            "get_ref_from_canoes_bed": self.get_ref_from_canoes_bed,
            "get_alt_from_canoes_bed": self.get_alt_from_canoes_bed,
            "get_chr_from_breakpoint": self.get_chr_from_breakpoint,
            "get_pos_from_breakpoint": self.get_pos_from_breakpoint,
            "get_ref_from_breakpoint": self.get_ref_from_breakpoint,
            "get_alt_from_star_breakpoint": self.get_alt_from_star_breakpoint,
            "get_alt_from_arriba_breakpoint": self.get_alt_from_arriba_breakpoint,
            "readable_starfusion_annots": self.readable_starfusion_annots,
            "get_undefined_value": self.get_undefined_value,
            "get_gt_from_varank": self.get_gt_from_varank,
            "get_ad_from_varank": self.get_ad_from_varank,
            "get_vaf_from_varank": self.get_vaf_from_varank,
            "get_chr_from_bedpe": self.get_chr_from_bedpe,
            "get_pos_from_bedpe": self.get_pos_from_bedpe,
            "get_ref_from_bedpe": self.get_ref_from_bedpe,
            "get_alt_from_bedpe": self.get_alt_from_bedpe,
        }

    def get(self, func_name):
        return self.dispatcher[func_name]

    def get_ref_from_decon(self, chrom, start):
        f = get_genome(self.config["GENOME"]["path"])
        if self.config["GENOME"]["vcf_header"][0].startswith(
            "##contig=<ID=chr"
        ) and not chrom.startswith("chr"):
            chrom = "chr" + str(chrom)
        return f[chrom][int(start) - 1].seq

    def get_ref_from_canoes_bed(self, chr, start):
        f = get_genome(self.config["GENOME"]["path"])
        return f["chr" + str(chr)][int(start) - 1].seq

    def get_ref_from_breakpoint(self, left_breakpoint, right_breakpoint):
        f = get_genome(self.config["GENOME"]["path"])

        left_chr = left_breakpoint.split(":")[0]
        if not left_chr.startswith("chr"):
            left_chr = "chr" + left_chr
        left_start = left_breakpoint.split(":")[1]

        right_chr = right_breakpoint.split(":")[0]
        if not right_chr.startswith("chr"):
            right_chr = "chr" + right_chr
        right_start = right_breakpoint.split(":")[1]

        return (
            f[left_chr][int(left_start) - 1].seq,
            f[right_chr][int(right_start) - 1].seq,
        )

    def get_alt_with_breakpoints(self, chr1, pos1, strand1, ref1, chr2, pos2, strand2, ref2):
        """
        See VCF 4.3 specification, section 5.4 "Specifying complex rearrangements with breakends"

        Returns a tuple of two strings: left ALT and right ALT fields
        """
        for chrom in (chr1, chr2):
            if not chrom.startswith("chr"):
                chrom = "chr" + str(chrom)

        if strand1 == "+" and strand2 == "+":
            alt1 = f"{ref1}[{chr2}:{pos2}["
            alt2 = f"]{chr1}:{pos1}]{ref2}"
        elif strand1 == "+" and strand2 == "-":
            alt1 = f"{ref1}]{chr2}:{pos2}]"
            alt2 = f"{ref2}]{chr1}:{pos1}]"
        elif strand1 == "-" and strand2 == "+":
            alt1 = f"[{chr2}:{pos2}[{ref1}"
            alt2 = f"[{chr1}:{pos1}[{ref2}"
        elif strand1 == "-" and strand2 == "-":
            alt1 = f"]{chr2}:{pos2}]{ref1}"
            alt2 = f"{ref2}[{chr1}:{pos1}["
        else:
            raise ValueError(
                "Strand should be + or -. Got values: strand1:"
                + str(strand1)
                + " ; strand2:"
                + str(strand2)
            )

        return (alt1, alt2)

    def get_alt_from_star_breakpoint(self, left_breakpoint, right_breakpoint):
        left_ref, right_ref = self.get_ref_from_breakpoint(left_breakpoint, right_breakpoint)
        left_chr, left_pos, left_strand = left_breakpoint.split(":")
        right_chr, right_pos, right_strand = right_breakpoint.split(":")

        return self.get_alt_with_breakpoints(
            left_chr,
            left_pos,
            left_strand,
            left_ref,
            right_chr,
            right_pos,
            right_strand,
            right_ref,
        )

    def get_alt_from_arriba_breakpoint(self, breakpoint1, breakpoint2, direction1, direction2):
        """
        In Arriba, unlike with STAR-Fusion, breakpoints aren't directly given as "left" and "right" breakpoints.
        This information can be inferred from the "direction" column.
        The inferrence method here should be equivalent to using the second strand from each breakpoint's strand column, when it is present.
        See Arriba documentation: https://arriba.readthedocs.io/en/latest/output-files/
        """
        ref1, ref2 = self.get_ref_from_breakpoint(breakpoint1, breakpoint2)
        chr1, pos1 = breakpoint1.split(":")
        chr2, pos2 = breakpoint2.split(":")

        if direction1 == "upstream" and direction2 == "upstream":
            strand1 = "-"
            strand2 = "+"
        elif direction1 == "upstream" and direction2 == "downstream":
            strand1 = "-"
            strand2 = "-"
        elif direction1 == "downstream" and direction2 == "upstream":
            strand1 = "+"
            strand2 = "+"
        elif direction1 == "downstream" and direction2 == "downstream":
            strand1 = "+"
            strand2 = "-"
        else:
            raise ValueError(
                "Direction should be upstream or downstream. Got values: direction1:"
                + str(direction1)
                + " ; direction2:"
                + str(direction2)
            )

        return self.get_alt_with_breakpoints(chr1, pos1, strand1, ref1, chr2, pos2, strand2, ref2)

    @staticmethod
    def get_alt_from_decon(cnv_type_field):
        if cnv_type_field == "deletion":
            return "<DEL>"
        if cnv_type_field == "duplication":
            return "<DUP>"
        raise ValueError("Unexpected CNV.type value:" + str(cnv_type_field))

    @staticmethod
    def get_alt_from_canoes_bed(cnv_type_field):
        if cnv_type_field == "DEL":
            return "<DEL>"
        if cnv_type_field == "DUP":
            return "<DUP>"
        raise ValueError("Unexpected CNV.type value:" + str(cnv_type_field))

    @staticmethod
    def get_svlen_from_decon(start, end):
        return str(int(end) - int(start))

    @staticmethod
    def get_info_from_annotsv(info):
        """
        only used in attempts to convert annotsv files
        as if they were generic TSV (not recommended)
        """
        return "."

    @staticmethod
    def get_alt_for_bed_based_annotsv(sv_type):
        return "<" + sv_type + ">"

    @staticmethod
    def get_chr_from_breakpoint(left_breakpoint, right_breakpoint):
        return (left_breakpoint.split(":")[0], right_breakpoint.split(":")[0])

    @staticmethod
    def get_pos_from_breakpoint(left_breakpoint, right_breakpoint):
        return (left_breakpoint.split(":")[1], right_breakpoint.split(":")[1])

    @staticmethod
    def readable_starfusion_annots(annots):
        """
        input: '["Mitelman","ChimerKB","GUO2018CR_TCGA","DEEPEST2019","HGNC_GENEFAM","Cosmic","ChimerSeq","INTERCHROMOSOMAL[chr11--chr10]"]'
        output: 'Mitelman,ChimerKB,GUO2018CR_TCGA,DEEPEST2019,HGNC_GENEFAM,Cosmic,ChimerSeq,INTERCHROMOSOMAL[chr11--chr10]'
        """
        return ",".join([v[1:-1] for v in annots[1:-1].split(",")])

    @staticmethod
    def get_undefined_value():
        return "."

    @staticmethod
    def get_gt_from_varank(zigosity):
        gt_dic = {"hom": "1/1", "het": "0/1"}
        return gt_dic[zigosity]

    @staticmethod
    def get_ad_from_varank(total_read_depth, var_read_depth):
        return str(int(total_read_depth) - int(var_read_depth)) + "," + var_read_depth

    @staticmethod
    def get_vaf_from_varank(var_read_percent):
        if var_read_percent == ".":
            return var_read_percent
        return str(float(var_read_percent) / 100)

    @staticmethod
    def get_chr_from_bedpe(chr1, chr2):
        return (chr1, chr2)

    @staticmethod
    def get_pos_from_bedpe(end1, end2):
        return (end1, end2)

    def get_ref_from_bedpe(self, chr1, end1, chr2, end2):
        f = get_genome(self.config["GENOME"]["path"])

        return (
            f[chr1][int(end1) - 1].seq,
            f[chr2][int(end2) - 1].seq,
        )

    def get_alt_from_bedpe(self, chr1, end1, strand1, chr2, end2, strand2):
        ref1, ref2 = self.get_ref_from_bedpe(chr1, end1, chr2, end2)

        return self.get_alt_with_breakpoints(
            chr1,
            end1,
            strand1,
            ref1,
            chr2,
            end2,
            strand2,
            ref2,
        )
