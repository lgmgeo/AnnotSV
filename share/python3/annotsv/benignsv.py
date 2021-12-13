from pathlib import Path

from annotsv.constants import annotation_dir
from annotsv.enums import GenomeBuild


def checkBenignFiles(organism, genome_build: GenomeBuild):
    benign_dir = annotation_dir / f"Annotations_{organism}/SVincludedInFt/BenignSV/{genome_build}"
    benign_files = {"k": Path()}
    ...


def checkClinVar_benignFile(genome_build: GenomeBuild):
    ...


def checkClinGenHITS_benignFile(genome_build: GenomeBuild):
    ...


def checkDGV_benignFile(genome_build: GenomeBuild):
    ...


def checkGnomAD_benignFile(genome_build: GenomeBuild):
    ...


def checkDDD_benignFile(genome_build: GenomeBuild):
    ...


def check1000g_benignFile(genome_build: GenomeBuild):
    ...


def checkIMH_benignFile(genome_build: GenomeBuild):
    ...


def benignSVannotation(SVchrom, SVstart, SVend):
    ...
