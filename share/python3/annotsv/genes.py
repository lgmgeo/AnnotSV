from __future__ import annotations
import gzip
from annotsv.context import Context
from annotsv.constants import VALID_CHROM


def check_genes_refseq_file(app: Context):
    """checks if the Genes file was downloaded, and sorts/formats it. no-op if using AnnotSV
    packaged annotation files"""
    gene_file_downloaded = app.config.genes_dir / "refGene.txt.gz"
    gene_file_formatted = app.config.genes_dir / "genes.RefSeq.sorted.bed"

    if not gene_file_downloaded.exists() and not gene_file_formatted.exists():
        app.abort(
            f"{gene_file_downloaded} and {gene_file_formatted} do not exist. Please check your install.",
        )

    if not gene_file_formatted.exists():
        ## Delete promoters files (need to be updated after the creation of new genes file)
        for prom_file in app.config.promoter_dir.glob(
            f"promoter_*bp_RefSeq_{app.config.genome_build}.sorted.bed"
        ):
            prom_file.unlink()

        ## - Create the "genes.RefSeq.sorted.bed"
        app.log.info(f"creating {gene_file_formatted} - done only once during the first annotation")

        # Removing non-standard contigs (other than the standard 1-22,X,Y,MT) and sorting the file in karyotypic order
        ## Save the line of the genesTXT by chromosome (L_lines($chrom))

        lines_by_chr = {}
        with gzip.open(gene_file_downloaded, "rt") as gf_in:
            for line in gf_in:
                cols = line.rstrip("\n").split("\t")
                cols[1] = cols[1].replace("chr", "")
                if cols[1] in VALID_CHROM:
                    lines_by_chr[cols[1]] = cols

        # Chromosome nomenclature used: 1-22,X,Y,M,T (without "chr")
        ## INPUT:    #bin name chrom strand txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds score name2 cdsStartStat cdsEndStat exonFrames
        ## OUTPUT:   chrom txStart txEnd strand name2 name cdsStart cdsEnd exonStarts exonEnds
        # WARNING : NR_* are for non coding RNA. However, cdsStart=cdsEnd, => CDSlength=1
        output_cols = [4, 5, 3, 12, 1, 6, 7, 9, 10]
        # tmp_format = gene_file_formatted.with_suffix(".tmp")
        with gene_file_formatted.open("rt") as gf_out:
            for chrom in VALID_CHROM:
                if chrom not in lines_by_chr:
                    continue
                # list.sort is in-place
                lines_by_chr[chrom].sort(key=lambda x: int(x[5]))
                lines_by_chr[chrom].sort(key=lambda x: int(x[4]))
                for line in lines_by_chr[chrom]:
                    l_str = "\t".join(line[i] for i in output_cols)
                    gf_out.write(f"{l_str}\n")

        gene_file_downloaded.unlink()

    app.genes_file = gene_file_formatted


def check_genes_ensembl_file(app: Context):
    "Checks if the formatted ENSEMBL genes file is present"
    gene_file_formatted = app.config.genes_dir / "genes.ENSEMBL.sorted.bed"

    if not gene_file_formatted.exists():
        app.abort(f"{gene_file_formatted} does not exist. Please check your install.")

    app.genes_file = gene_file_formatted


## Annotate the SV bedFile with the genes file.
## Keep only 1 transcript annotation by gene:
##   - the one selected by the user with the "-txFile" option
##   - the one with the most of "bp from CDS" (=overlappedCDSlength)
##   - if x transcript with same "bp from CDS", the one with the most of "bp from UTR, exon, intron" (=overlappedTxLength)
##
## Creation of FullAndSplitBedFile ($g_AnnotSV(outputDir)/$g_AnnotSV(outputFile).tmp)
## -> formatted and sorted
def genes_annotation(app: Context):
    ...
