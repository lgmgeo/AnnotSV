from __future__ import annotations

from typing import TYPE_CHECKING
from annotsv.schemas import AnnotationValidator, ResolvedFiles
from annotsv.enums import GenomeBuild

if TYPE_CHECKING:
    from annotsv.context import Context

BENIGN_SVTYPES = ["Loss", "Gain", "Ins", "Inv"]


class BenignValidator(AnnotationValidator):
    def __init__(self, app: Context):
        # , *, label: str, downloaded: ResolvedFiles, formatted: ResolvedFiles, extra_downloaded: List[ResolvedFiles] = None, extra_formatted: List[ResolvedFiles] = None):
        # super().__init__(app, label, downloaded, formatted, extra_downloaded=extra_downloaded, extra_formatted=extra_formatted)
        ...


# Creation / update (if some data source files are presents) of:
# - $benignLossFile_Sorted
# - $benignGainFile_Sorted
# - $benignInsFile_Sorted
# - $benignInvFile_Sorted
def check_benign_files(app: Context):
    for build in [GenomeBuild.GRCh37, GenomeBuild.GRCh38]:
        check_gnomad_file(app, build)
        check_dgv_file(app, build)
        check_ddd_file(app, build)
        check_1000g_file(app, build)
        check_clingen_hits_file(app, build)
        check_clinvar_file(app, build)
        check_imh_file(app, build)
        check_cmri_file(app, build)

        sorted_files = _benign_files(app, build)
        tmp_files = _benign_files(app, build, tmp=True)
        for idx, sv_type in enumerate(BENIGN_SVTYPES):
            if tmp_files[idx].exists():
                raise NotImplementedError()


def check_1000g_file(app: Context, build: GenomeBuild):
    downloaded_files = list(app.config.benign_root.glob(f"{build}/ALL.wgs.mergedSV*.vcf.gz"))

    if downloaded_files:
        raise NotImplementedError()
    else:
        app.log.debug(f"No new benign {build} 1000g annotation to format")


def check_clingen_hits_file(app: Context, build: GenomeBuild):
    downloaded_gene_file = (
        app.config.benign_root / f"{build}/ClinGen_gene_curation_list_{build}.tsv"
    )
    downloaded_region_file = (
        app.config.benign_root / f"{build}/ClinGen_region_curation_list_{build}.tsv"
    )

    if downloaded_gene_file.exists() or downloaded_region_file.exists():
        raise NotImplementedError()
    else:
        app.log.debug(f"No new benign {build} clingen annotation to format")


def check_clinvar_file(app: Context, build: GenomeBuild):
    files_downloaded = list(app.config.benign_root.glob(f"{build}/clinvar*.vcf.gz"))

    if files_downloaded:
        raise NotImplementedError()
    else:
        app.log.debug(f"No new benign {build} clinvar annotation to format")


def check_cmri_file(app: Context, build: GenomeBuild):
    downloaded_file = app.config.benign_root / f"{build}/pb_joint_merged.sv.vcf"

    if downloaded_file.exists():
        raise NotImplementedError()
    else:
        app.log.debug(f"No new benign {build} CMRI annotation to format")


def check_ddd_file(app: Context, build: GenomeBuild):
    downloaded_file = app.config.benign_root / f"{build}/population_cnv_{build.name.lower()}.txt.gz"

    if downloaded_file.exists():
        raise NotImplementedError()
    else:
        app.log.debug(f"No new benign {build} DDD annotation to format")


def check_dgv_file(app: Context, build: GenomeBuild):
    downloaded_files = list(app.config.benign_root.glob(f"{build}/GRCh3*_hg*_variants_*.txt"))

    if downloaded_files:
        raise NotImplementedError()
    else:
        app.log.debug(f"No new benign {build} DGV annotation to format")


def check_gnomad_file(app: Context, build: GenomeBuild):
    # GRCh37 only (not yet available in GRCh38, 2021/10/05)
    downloaded_file = app.config.benign_root / f"{build}/gnomad_v2.1_sv.sites.bed.gz"

    if downloaded_file.exists():
        raise NotImplementedError()
    else:
        app.log.debug(f"No new benign {build} gnomAD annotation to format")


def check_imh_file(app: Context, build: GenomeBuild):
    downloaded_files = list(app.config.benign_root.glob(f"{build}/*callset.public.bedpe.gz"))

    if downloaded_files:
        raise NotImplementedError()
    else:
        app.log.debug(f"No new benign {build} IMH annotation to format")


def _benign_files(app: Context, build: GenomeBuild, tmp: bool = False):
    "returns Path objects for Loss, Gain, Ins, Inv files"
    files = []
    if tmp:
        ftype = "tmp"
    else:
        ftype = "sorted"
    for sv_type in BENIGN_SVTYPES:
        files.append(app.config.benign_root / f"{build}/{sv_type}_SV_{build}.{ftype}.bed")
    return files


def benignSVannotation(app: Context, SVchrom, SVstart, SVend):
    ...
