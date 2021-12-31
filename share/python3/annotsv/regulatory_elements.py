from __future__ import annotations
from typing import List

from annotsv.context import Context
from annotsv.enums import Organisms


def is_refseq_gene_name(app: Context, gene_name: str):
    if not app.refseq_genes:
        refseq_file = app.config.genes_dir / "genes.RefSeq.sorted.bed"
        with refseq_file.open("rt") as fh:
            for line in fh:
                if line.strip():
                    app.refseq_genes.add(line.split("\t")[4])

    return gene_name in app.refseq_genes


def is_ensembel_gene_name(app: Context, gene_name: str):
    if not app.ensembl_genes:
        ensembl_file = app.config.genes_dir / "genes.ENSEMBL.sorted.bed"
        with ensembl_file.open("rt") as fh:
            for line in fh:
                if line.strip():
                    app.ensembl_genes.add(line.split("\t")[4])

    return gene_name in app.ensembl_genes


## - Check and create if necessary the "promoter_XXXbp_*_GRCh*.sorted.bed" file.
def check_promoter_file(app: Context):
    formatted_file = (
        app.config.reg_elements_dir
        / f"promoter_{app.config.promoter_size}bp_{app.config.tx}_{app.config.genome_build}.sorted.bed"
    )

    if not formatted_file.exists():
        raise NotImplementedError()

    app.promoter_ann = True


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
def check_ea_files(app: Context):
    # GRCh38 version should be manually created by lift over
    formatted_refseq = (
        app.config.reg_elements_dir / f"EA_RefSeq_{app.config.genome_build}.sorted.bed"
    )
    formatted_ensembl = (
        app.config.reg_elements_dir / f"EA_ENSEMBL_{app.config.genome_build}.sorted.bed"
    )
    ea_file = (
        app.config.reg_elements_dir / f"EA_{app.config.tx}_{app.config.genome_build}.sorted.bed"
    )
    downloaded_files = list(app.config.reg_elements_dir.glob("*_EP.txt"))
    label = "EnhancerAtlas"

    if ea_file.exists():
        app.log.debug(f"Enabling {label} annotation")
        app.ea_ann = True

    if formatted_refseq.exists() and formatted_ensembl.exists():
        app.log.debug(f"Using existing {label} annotation files")
    elif downloaded_files:
        raise NotImplementedError()
    else:
        app.log.debug(f"No {label} annotation")
        assert app.ea_ann is False, f"{label} annotation is enabled, but no annotation files found"


## - Check if the following GH files has been downloaded:
##   - GeneHancer_elements.txt
##   - GeneHancer_gene_associations_scores.txt
##   - GeneHancer_hg19.txt
##
## - Check and create if necessary:
##   - GH_RefSeq_GRCh37.sorted.bed
##   - GH_RefSeq_GRCh38.sorted.bed
##   - GH_ENSEMBL_GRCh37.sorted.bed
##   - GH_ENSEMBL_GRCh38.sorted.bed
def check_gh_files(app: Context):
    refseq_file = app.config.reg_elements_dir / f"GH_RefSeq_{app.config.genome_build}.sorted.bed"
    ensembl_file = app.config.reg_elements_dir / f"GH_ENSEMBL_{app.config.genome_build}.sorted.bed"
    elements_file = app.config.reg_elements_dir / "GeneHancer_elements.txt"
    associations_file = app.config.reg_elements_dir / "GeneHancer_gene_associations_scores.txt"
    hg19_file = app.config.reg_elements_dir / "GeneHancer_hg19.txt"
    label = "GeneHancer"

    if app.config.organism is not Organisms.Human:
        app.log.debug(f"No {label} annotation for {app.config.organism}, ignoring")
    elif refseq_file.exists() and ensembl_file.exists():
        app.log.debug(f"Enabling {label} annotation")
        app.gh_ann = True
    elif any(not f.exists() for f in [elements_file, associations_file, hg19_file]):
        app.log.warning(
            f"No {label} annotations available. Please, see in the README file how to add these annotations. Users need to contact the GeneCards team."
        )
    else:
        raise NotImplementedError()


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
def check_mir_target_link_files(app: Context):
    refseq_file = (
        app.config.reg_elements_dir / f"miRTargetLink_RefSeq_{app.config.genome_build}.sorted.bed"
    )
    ensembl_file = (
        app.config.reg_elements_dir / f"miRTargetLink_ENSEMBL_{app.config.genome_build}.sorted.bed"
    )
    label = "miRTargetLink"

    if (
        app.config.organism is Organisms.Human and refseq_file.exists() and ensembl_file.exists()
    ) or (app.config.organism is Organisms.Mouse and refseq_file.exists()):
        app.log.debug(f"Enabling {app.config.genome_build} {label} annotation")
        app.mirna_ann = True
    else:
        raise NotImplementedError()


def regulatory_elements_annotation(app: Context, overlapped_genes: List[str]):
    ...
