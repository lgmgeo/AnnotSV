from pathlib import Path

from annotsv.context import Context
from annotsv.enums import GenomeBuild


def check_benign_files(app: Context):
    for build in GenomeBuild:
        benign_dir = app.config.benign_dir / build.name
        benign_loss_file = benign_dir / f"benign_Loss_SV_{build}.sorted.bed"
        benign_gain_file = benign_dir / f"benign_Gain_SV_{build}.sorted.bed"
        benign_ins_file = benign_dir / f"benign_Ins_SV_{build}.sorted.bed"
        benign_inv_file = benign_dir / f"benign_Inv_SV_{build}.sorted.bed"

        check_gnomad_file(app)
        check_dgv_file(app)
        check_ddd_file(app)
        check_1000g_file(app)
        check_clingen_hits_file(app)
        check_clinvar_file(app)
        check_imh_file(app)
        check_cmri_file(app)


def check_clinvar_file(app: Context):
    ...


def check_clingen_hits_file(app: Context):
    ...


def check_cmri_file(app: Context):
    ...


def check_dgv_file(app: Context):
    ...


def check_gnomad_file(app: Context):
    ...


def check_ddd_file(app: Context):
    ...


def check_1000g_file(app: Context):
    ...


def check_imh_file(app: Context):
    ...


def _tmp_benign_files(app: Context, build: GenomeBuild):
    benign_dir = app.config.benign_dir / build.name
    for ftype in ["Loss", "Gain", "Ins", "inv"]:
        ...


def benignSVannotation(SVchrom, SVstart, SVend):
    ...
