from __future__ import annotations
from os import dup
from pathlib import Path
from typing import TYPE_CHECKING, List
from annotsv.enums import GenomeBuild
from annotsv.schemas import AnnotationValidator, ResolvedFiles

if TYPE_CHECKING:
    from annotsv.context import Context

PATHOGENIC_SVTYPES = ["Loss", "Gain", "Ins", "Inv"]


class PathogenicValidator(AnnotationValidator):
    def __init__(
        self,
        app: Context,
        *,
        label: str,
        downloaded_pattern: str,
        extra_patterns: List[str] = None,
        root_dir: Path = None,
    ):
        self._app = app
        self.label = f"pathogenic {label}"
        self.formatted = []
        self.downloaded = []
        if root_dir is None:
            root_dir = self._app.config.pathogenic_sv_root
        plist = [downloaded_pattern]
        if extra_patterns:
            plist.extend(extra_patterns)
        for build in [GenomeBuild.GRCh37, GenomeBuild.GRCh38]:
            for pattern in plist:
                # some patterns also have build in filename, so do a str.format just in case
                self.downloaded.append(
                    ResolvedFiles(
                        root_dir / build.name,
                        pattern.format(build=build.name, build_lower=build.name.lower()),
                    )
                )

    def check(self):
        if self.downloaded_exist():
            self.update()
        else:
            self._app.log.debug(f"No new {self.label} annotation to format")
        return True


class ClingenValidator(PathogenicValidator):
    def __init__(self, app: Context):
        super().__init__(
            app,
            label="ClinGen HITS",
            downloaded_pattern="ClinGen_gene_curation_list_{build}.tsv",
            extra_patterns=["ClinGen_region_curation_list_{build}.tsv"],
        )

    def update(self):
        raise NotImplementedError()


class ClinvarValidator(PathogenicValidator):
    def __init__(self, app: Context):
        super().__init__(app, label="ClinVar", downloaded_pattern="clinvar*.vcf.gz")

    def update(self):
        raise NotImplementedError()


class DBVarValidator(PathogenicValidator):
    def __init__(self, app: Context):
        del_pattern = "{build}.nr_deletions.pathogenic.tsv.gz"
        dup_pattern = "{build}.nr_duplications.pathogenic.tsv.gz"
        ins_pattern = "{build}.nr_insertions.pathogenic.tsv.gz"
        super().__init__(
            app,
            label="dbVar",
            downloaded_pattern=del_pattern,
            extra_patterns=[dup_pattern, ins_pattern],
        )

    def update(self):
        raise NotImplementedError()


class OMIMValidator(PathogenicValidator):
    def __init__(self, app: Context):
        super().__init__(app, label="OMIM", downloaded_pattern="*_morbid.tsv.gz")

    def update(self):
        raise NotImplementedError()


class SnvIndelValidator(PathogenicValidator):
    def __init__(self, app: Context):
        super().__init__(
            app,
            label="SNV indel",
            downloaded_pattern="clinvar*.vcf.gz",
            root_dir=app.config.pathogenic_snv_root,
        )

    def update(self):
        raise NotImplementedError()


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
