from __future__ import annotations

from annotsv.context import Context
from annotsv.enums import GenomeBuild

PATHOGENIC_SVTYPES = ["Loss", "Gain", "Ins", "Inv"]


def check_pathogenic_files(app: Context):
    for build in [GenomeBuild.GRCh37, GenomeBuild.GRCh37]:
        check_clinvar_file(app, build)
        check_clingen_hits_file(app, build)
        check_dbvar_file(app, build)
        check_snv_indel_file(app, build)


def check_clinvar_file(app: Context, build: GenomeBuild):
    downloaded_files = list(app.config.pathogenic_sv_root.glob("{build}/clinvar*.vcf.gz"))
    label = "ClinVar"

    if downloaded_files:
        raise NotImplementedError()
    else:
        app.log.debug(f"No new pathogenic {label} annotation to format")


def check_clingen_hits_file(app: Context, build: GenomeBuild):
    gene_file = app.config.pathogenic_sv_root / f"{build}/ClinGen_gene_curation_list_{build}.tsv"
    region_file = (
        app.config.pathogenic_sv_root / f"{build}/ClinGen_region_curation_list_{build}.tsv"
    )
    label = "ClinGen HITS"

    if gene_file.exists() or region_file.exists():
        raise NotImplementedError()
    else:
        app.log.debug(f"No new pathogenic {label} annotation to format")


def check_dbvar_file(app: Context, build: GenomeBuild):
    del_file = (
        app.config.pathogenic_sv_root / f"{build}/{build}.nr_deletions.pathogenic.tsv.gz"
    )  # Loss
    dup_file = (
        app.config.pathogenic_sv_root / f"{build}/{build}.nr_duplications.pathogenic.tsv.gz"
    )  # Gain
    ins_file = (
        app.config.pathogenic_sv_root / f"{build}/{build}.nr_insertions.pathogenic.tsv.gz"
    )  # Ins
    label = "dbVar"

    if del_file.exists() or dup_file.exists() or ins_file.exists():
        raise NotImplementedError()
    else:
        app.log.debug(f"No new pathogenic {label} annotation to format")


def check_omim_file(app: Context, build: GenomeBuild):
    downloaded_files = list(app.config.pathogenic_sv_root.glob(f"{build}/*_morbid.tsv.gz"))
    label = "OMIM"

    if downloaded_files:
        raise NotImplementedError()
    else:
        app.log.debug(f"No new pathogenic {label} annotation to format")


def check_snv_indel_file(app: Context, build: GenomeBuild):
    downloaded_files = list(app.config.pathogenic_snv_root.glob(f"{build}/clinvar*.vcf.gz"))
    label = "SNV indel"

    if downloaded_files:
        raise NotImplementedError()
    else:
        app.log.debug(f"No new pathogenic {label} annotation to format")


def _pathogenic_sv_files(app: Context, build: GenomeBuild, tmp: bool = False):
    "returns Path objects for Loss, Gain, Ins, Inv files"
    files = []
    if tmp:
        ftype = "tmp"
    else:
        ftype = "sorted"
    for sv_type in PATHOGENIC_SVTYPES:
        files.append(
            app.config.pathogenic_sv_root / f"{build}/pathogenic_{sv_type}_SV_{build}.{ftype}.bed"
        )
    return files


def pathogenic_sv_annotation(app: Context, sv_chrom: str, sv_start: int, sv_end: int):
    ...


def pathogenic_snv_indel_annotation(app: Context, sv_chrom: str, sv_start: int, sv_end: int):
    ...
