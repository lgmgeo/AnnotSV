from __future__ import annotations

from typing import TYPE_CHECKING, List

from annotsv.enums import GenomeBuild, Organisms, TranscriptSource
from annotsv.schemas import AnnotationValidator, ResolvedFiles

if TYPE_CHECKING:
    from annotsv.context import Context


class PromoterValidator(AnnotationValidator):
    def __init__(self, app: Context):
        self.label = "promoter"
        self.downloaded = []
        self.formatted = [
            ResolvedFiles(
                app.config.reg_elements_dir,
                f"promoter_{app.config.promoter_size}bp_{app.config.tx}_{app.config.genome_build}.sorted.bed",
            )
        ]

    def check(self):
        if not self.formatted_exist():
            self.update()
        return True

    def update(self):
        raise NotImplementedError()


##  EnhancerAtlas
#################
## - Check if some "*_EP.txt" files have been downloaded
##
## - Check and create if necessary:
##   - EA_RefSeq_GRCh37.sorted.bed
##   - EA_ENSEMBL_GRCh37.sorted.bed
##
## - The GRCh38 version should be manually created by lift over with the UCSC web server, sorted, and
##   move in the “$ANNOTSV/share/AnnotSV/Annotations_Human/FtIncludedInSV/RegulatoryElements/GRCh38” directory.
class EnhancerAtlasValidator(AnnotationValidator):
    def __init__(self, app: Context):
        super().__init__(
            app,
            label="EnhancerAtlas",
            downloaded=ResolvedFiles(app.config.reg_elements_dir, "*_EP.txt"),
            formatted=ResolvedFiles(
                app.config.reg_elements_dir,
                f"EA_{app.config.tx}_{app.config.genome_build}.sorted.bed",
            ),
        )

    def update(self):
        raise NotImplementedError()


class GenehancerValidator(AnnotationValidator):
    def __init__(self, app: Context):
        rf_dir = app.config.reg_elements_dir
        rf_refseq = ResolvedFiles(rf_dir, f"GH_RefSeq_{app.config.genome_build}.sorted.bed")
        rf_ensembl = ResolvedFiles(rf_dir, f"GH_ENSEMBL_{app.config.genome_build}.sorted.bed")
        rf_elements = ResolvedFiles(rf_dir, "GeneHancer_elements.txt")
        rf_associations = ResolvedFiles(rf_dir, "GeneHancer_gene_associations_scores.txt")
        rf_hg19 = ResolvedFiles(rf_dir, "GeneHancer_hg19.txt")
        super().__init__(
            app,
            label="GeneHancer",
            downloaded=rf_elements,
            extra_downloaded=[rf_associations, rf_hg19],
            formatted=rf_refseq,
            extra_formatted=[rf_ensembl],
        )

    def update(self):
        raise NotImplementedError()

    def check(self):
        success = True
        if self._app.config.organism is not Organisms.Human:
            self._app.log.debug(
                f"No {self.label} annotation for {self._app.config.organism}, ignoring"
            )
            success = False
        else:
            success = super().check()
            if not success:
                self._app.log.warning(
                    f"No {self.label} annotations available. Please, see in the README file how to add these annotations. Users need to contact the GeneCards team."
                )

        return success


## Human:
#########
## - Check if the following miRTargetLink files has been downloaded:
##   - Validated_miRNA-gene_pairs_hsa_miRBase_v22.1_GRCh38_location_augmented.tsv
##
## - Check and create if necessary:
##   - miRTargetLink_RefSeq_GRCh38.sorted.bed
##   - miRTargetLink_ENSEMBL_GRCh38.sorted.bed
##
## - GRCh37 miRTargetLink are not provided (only GRCh38)
##   - miRTargetLink_ENSEMBL_GRCh37.sorted.bed => created with a UCSC liftover
##   - miRTargetLink_RefSeq_GRCh37.sorted.bed  => created with a UCSC liftover
#
## Mouse:
#########
## - Check if the following miRTargetLink files has been downloaded:
##   - Validated_miRNA-gene_pairs_mmu_miRBase_v22.1_GRCm38_location_augmented.tsv
##
## - Check and create if necessary:
##   - miRTargetLink_RefSeq_mm10.sorted.bed
##
## - mm9 miRTargetLink is not provided (only mm10)
##   - miRTargetLink_RefSeq_mm9.sorted.bed  => created with a UCSC liftover
class MirTargetLinkValidator(AnnotationValidator):
    def __init__(self, app: Context):
        # GRCh38 or GRCm38
        build_str = f"GRC{app.config.organism[0].lower()}38"
        super().__init__(
            app,
            label="miRTargetLink",
            downloaded=ResolvedFiles(
                app.config.reg_elements_dir,
                f"*_miRNA-gene_pairs_hsa_miRBase_*_GRC{build_str}38_location_augmented.tsv",
            ),
            formatted=ResolvedFiles(
                app.config.reg_elements_dir,
                f"miRTargetLink_{app.config.tx}_{app.config.genome_build}.sorted.bed",
            ),
        )

    def check(self):
        success = True
        if self._app.config.genome_build not in [GenomeBuild.GRCh38, GenomeBuild.mm10]:
            self._app.log.debug(
                f"{self.label} annotation not available for {self._app.config.genome_build}, ignoring"
            )
            success = False
        elif (
            self._app.config.genome_build is GenomeBuild.mm10
            and self._app.config.tx is not TranscriptSource.RefSeq
        ):
            self._app.log.debug(
                f"{self.label} annotation not available for {self._app.config.genome_build}/{self._app.config.tx}, ignoring"
            )
            success = False
        else:
            success = super().check()

        return success

    def update(self):
        raise NotImplementedError()


def is_refseq_gene_name(app: Context, gene_name: str):
    if not app.refseq_genes:
        refseq_file = app.config.genes_dir / "genes.RefSeq.sorted.bed"
        with refseq_file.open("rt") as fh:
            for line in fh:
                if line.strip():
                    app.refseq_genes.add(line.split("\t")[4])

    return gene_name in app.refseq_genes


def is_ensembl_gene_name(app: Context, gene_name: str):
    if not app.ensembl_genes:
        ensembl_file = app.config.genes_dir / "genes.ENSEMBL.sorted.bed"
        with ensembl_file.open("rt") as fh:
            for line in fh:
                if line.strip():
                    app.ensembl_genes.add(line.split("\t")[4])

    return gene_name in app.ensembl_genes


def regulatory_elements_annotation(app: Context, overlapped_genes: List[str]):
    ...
