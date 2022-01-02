from __future__ import annotations

from typing import TYPE_CHECKING
from annotsv.schemas import AnnotationValidator, ResolvedFiles

if TYPE_CHECKING:
    from annotsv.context import Context

# ExAC downloaded file: "fordist_cleaned_nonpsych_z_pli_rec_null_data.txt"
# Header:
# transcript  gene  chr  n_exons tx_start  tx_end  bp   p_syn  p_mis   p_lof   n_syn   n_mis   n_lof   adj_exp_syn  adj_exp_mis   adj_exp_lof   syn_z   mis_z  lof_z   pLI   pRecessive  pNull

## - Check if the ExAC file has been downloaded:
#    - fordist_cleaned_nonpsych_z_pli_rec_null_data.txt
#
## - Check and create if necessary the following file:
#    - 'date'_ExAC-Zscore.annotations.tsv.gz
def check_gene_intolerance_file(app: Context):
    downloaded_file = app.config.exac_dir / "fordist_cleaned_nonpsych_z_pli_rec_null_data.txt"
    formatted_files = list(app.config.exac_dir.glob("*_GeneIntolerance-Zscore.annotations.tsv.gz"))

    if not downloaded_file.exists() and not formatted_files:
        app.log.debug(f"No gene intolerance annotation")
    elif len(formatted_files) > 1:
        app.keep_last_file("Gene Intolerance", formatted_files)
    elif not formatted_files:
        raise NotImplementedError()
    else:
        app.log.debug("No new gene intolerance annotation to format")


class GeneIntoleranceValidator(AnnotationValidator):
    def __init__(self, app: Context):
        super().__init__(
            app,
            label="Gene Inolerance",
            downloaded=ResolvedFiles(
                app.config.exac_dir, "fordist_cleaned_nonpsych_z_pli_rec_null_data.txt"
            ),
            formatted=ResolvedFiles(
                app.config.exac_dir, "*_GeneIntolerance-Zscore.annotations.tsv.gz"
            ),
        )

    def update(self):
        raise NotImplementedError()


# gene: Gencode ID
# chr: chromosome
# start: genomic coordinate gene start
# end: genomic coordinate gene end
# gene_symbol: HUGO gene symbol
# mean_rd: Average read depth of gene across all individuals
# gc_content: Proportion of gene that is GC
# complexity: sequence complexity
# cds_len: Number of targeted coding bases
# gene_length: Total gene length
# num_targ: Number of targets of the gene included in CNV calling
# segdups: number of pairs of segmental duplications gene is within
# dip: Number of confident diploid calls
# del: Number of confident deletion calls
# dup: Number of confident duplication calls
# del.sing: Number of confident deletion calls spanning only a single gene
# dup.sing: Number of confident duplication calls spanning only a single gene
# del.sing.score: Winsorised single-gene deletion intolerance z-score
# dup.sing.score: Winsorised single-gene duplication intolerance z-score
# del.score: Winsorised deletion intolerance z-score
# dup.score: Winsorised duplication intolerance z-score
# cnv.score: Winsorised cnv intolerance z-score
# flag: Gene is in a known region of recurrent CNVs mediated by tandem segmental duplications and intolerance scores are more likely to be biased or noisy.
#

## - Check if the ExAC file has been downloaded:
#    - exac-final-cnv.gene.scores071316
#
## - Check and create if necessary the following file:
#    - 'date'_ExAC.CNV-Zscore.annotations.tsv.gz
def check_cnv_intolerance_file(app: Context):
    downloaded_file = app.config.exac_dir / "exac-final-cnv.gene.scores071316"
    formatted_files = list(app.config.exac_dir.glob("*_ExAC.CNV-Zscore.annotations.tsv.gz"))

    if not downloaded_file.exists() and not formatted_files:
        app.log.debug("No cnv intolerance annotation")
    elif len(formatted_files) > 1:
        app.keep_last_file("CNV Intolerance", formatted_files)
    elif not formatted_files:
        raise NotImplementedError()
    else:
        app.log.debug("No new cnv intolerance annotation to format")


class CNVIntoleranceValidator(AnnotationValidator):
    def __init__(self, app: Context):
        super().__init__(
            app,
            label="CNV Intolerance",
            downloaded=ResolvedFiles(app.config.exac_dir, "exac-final-cnv.gene.scores071316"),
            formatted=ResolvedFiles(app.config.exac_dir, "*_ExAC.CNV-Zscore.annotations.tsv.gz"),
        )

    def update(self):
        raise NotImplementedError()
