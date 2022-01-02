from __future__ import annotations

from typing import TYPE_CHECKING, List
from annotsv.schemas import AnnotationValidator, ResolvedFiles
from annotsv.enums import GenomeBuild

if TYPE_CHECKING:
    from annotsv.context import Context

BENIGN_SVTYPES = ["Loss", "Gain", "Ins", "Inv"]


class BenignValidator(AnnotationValidator):
    def __init__(
        self,
        app: Context,
        *,
        label: str,
        downloaded_pattern: str,
        extra_patterns: List[str] = None,
    ):
        self._app = app
        self.label = f"benign {label}"
        self.formatted = []
        self.downloaded = []
        plist = [downloaded_pattern]
        if extra_patterns:
            plist.extend(extra_patterns)
        for build in [GenomeBuild.GRCh37, GenomeBuild.GRCh38]:
            for pattern in plist:
                # some patterns also have build in filename, so do a str.format just in case
                self.downloaded.append(
                    ResolvedFiles(
                        self._app.config.benign_root / build.name,
                        pattern.format(build=build.name, build_lower=build.name.lower()),
                    )
                )

    def check(self):
        if self.downloaded_exist():
            self.update()
        else:
            self._app.log.debug(f"No new {self.label} annotation to format")
        return True


class ClingenValidator(BenignValidator):
    def __init__(self, app: Context):
        super().__init__(
            app,
            label="ClinGen HITS",
            downloaded_pattern="ClinGen_gene_curation_list_{build}.tsv",
            extra_patterns=["ClinGen_region_curation_list_{build}.tsv"],
        )

    def update(self):
        raise NotImplementedError()


class ClinvarValidator(BenignValidator):
    def __init__(self, app: Context):
        super().__init__(app, label="ClinVar", downloaded_pattern="clinvar*.vcf.gz")

    def update(self):
        raise NotImplementedError()


class CMRIValidator(BenignValidator):
    def __init__(self, app: Context):
        super().__init__(app, label="CMRI", downloaded_pattern="pb_joint_merged.sv.vcf")

    def update(self):
        raise NotImplementedError()


class DDDValidator(BenignValidator):
    def __init__(self, app: Context):
        super().__init__(app, label="DDD", downloaded_pattern="population_cnv_{build_lower}.txt.gz")

    def update(self):
        raise NotImplementedError()


class DGVValidator(BenignValidator):
    def __init__(self, app: Context):
        super().__init__(app, label="DGV", downloaded_pattern="GRCh3*_hg*_variants_*.txt")

    def update(self):
        raise NotImplementedError()


class GnomadValidator(BenignValidator):
    # GRCh37 only (not yet available in GRCh38, 2021/10/05)
    def __init__(self, app: Context):
        super().__init__(app, label="gnomAD", downloaded_pattern="gnomad_v2.1_sv.sites.bed.gz")

    def update(self):
        raise NotImplementedError()


class IMHValidator(BenignValidator):
    def __init__(self, app: Context):
        super().__init__(app, label="IMH", downloaded_pattern="*callset.public.bedpe.gz")

    def update(self):
        raise NotImplementedError()


class ThousandGenomesValidator(BenignValidator):
    def __init__(self, app: Context):
        super().__init__(app, label="1000g", downloaded_pattern="ALL.wgs.mergedSV*.vcf.gz")

    def update(self):
        raise NotImplementedError()


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
