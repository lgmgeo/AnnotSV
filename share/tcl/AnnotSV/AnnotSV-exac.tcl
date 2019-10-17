############################################################################################################
# AnnotSV 2.2.4                                                                                            #
#                                                                                                          #
# AnnotSV: An integrated tool for Structural Variations annotation and ranking                             #
#                                                                                                          #
# Copyright (C) 2017-2019 Veronique Geoffroy (veronique.geoffroy@inserm.fr)                                #
#                                                                                                          #
# This is part of AnnotSV source code.                                                                     #
#                                                                                                          #
# This program is free software; you can redistribute it and/or                                            #
# modify it under the terms of the GNU General Public License                                              #
# as published by the Free Software Foundation; either version 3                                           #
# of the License, or (at your option) any later version.                                                   #
#                                                                                                          #
# This program is distributed in the hope that it will be useful,                                          #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                           #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                            #
# GNU General Public License for more details.                                                             #
#                                                                                                          #
# You should have received a copy of the GNU General Public License                                        #
# along with this program; If not, see <http://www.gnu.org/licenses/>.                                     #
############################################################################################################


# ExAC downloaded file: "fordist_cleaned_nonpsych_z_pli_rec_null_data.txt"
# Header:
# transcript  gene  chr  n_exons tx_start  tx_end  bp   p_syn  p_mis   p_lof   n_syn   n_mis   n_lof   adj_exp_syn  adj_exp_mis   adj_exp_lof   syn_z   mis_z  lof_z   pLI   pRecessive  pNull

## - Check if the ExAC file has been downloaded:
#    - fordist_cleaned_nonpsych_z_pli_rec_null_data.txt
#
## - Check and create if necessary the following file:
#    - 'date'_ExAC.pLI-Zscore.annotations.tsv.gz
proc checkGeneIntoleranceFile {} {

    global g_AnnotSV

    ## Check if the GeneIntolerance file has been downloaded then formatted
    ######################################################################
    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Genes-based"

    set GeneIntoleranceFileDownloaded [glob -nocomplain "$extannDir/ExAC/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt"]
    set GeneIntoleranceFileFormattedGzip [glob -nocomplain "$extannDir/ExAC/*_GeneIntolerance.pLI-Zscore.annotations.tsv.gz"]

    if {$GeneIntoleranceFileDownloaded eq "" && $GeneIntoleranceFileFormattedGzip eq ""} {
	# No "Gene Intolerance" annotation
	return
    }

    if {[llength $GeneIntoleranceFileFormattedGzip]>1} {
	puts "Several Gene Intolerance files exist:"
	puts "$GeneIntoleranceFileFormattedGzip"
	puts "Keep only one: [lindex $GeneIntoleranceFileFormattedGzip end]\n"
	foreach gi [lrange $GeneIntoleranceFileFormattedGzip 0 end-1] {
	    file rename -force $gi $gi.notused
	}
	return
    }

    if {$GeneIntoleranceFileFormattedGzip eq ""} {
	## Create : 'date'_GeneIntolerance.pLI-Zscore.annotations.tsv.gz   ; # Header: chr start end syn_z mis_z pLI
	set GeneIntoleranceFileFormatted "$extannDir/ExAC/[clock format [clock seconds] -format "%Y%m%d"]_GeneIntolerance.pLI-Zscore.annotations.tsv"

	puts "...GeneIntolerance configuration for pLI and Z scores annotation from ExAC ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"

	puts "\t...creation of $GeneIntoleranceFileFormatted.gz"
	puts "\t   (done only once during the first GeneIntolerance annotation)\n"

	set TexteToWrite {genes\tsynZ_ExAC\tmisZ_ExAC\tpLI_ExAC}
	foreach L [LinesFromFile $GeneIntoleranceFileDownloaded] {
	    if {[regexp  "^transcript" $L]} {
		set i_gene    [lsearch -regexp $L "^gene"]; if {$i_gene == -1} {puts "Bad syntax into $GeneIntoleranceFileDownloaded.\ngene field not found - Exit with error"; exit 2}
		set i_synZ    [lsearch -regexp $L "^syn_z"]; if {$i_synZ == -1} {puts "Bad syntax into $GeneIntoleranceFileDownloaded.\nsyn_z field not found - Exit with error"; exit 2}
		set i_misZ    [lsearch -regexp $L "^mis_z"]; if {$i_misZ == -1} {puts "Bad syntax into $GeneIntoleranceFileDownloaded.\nmis_z field not found - Exit with error"; exit 2}
		set i_pLI     [lsearch -regexp $L "^pLI"];   if {$i_pLI == -1} {puts "Bad syntax into $GeneIntoleranceFileDownloaded.\npLI field not found - Exit with error"; exit 2}
		continue
	    }
	    set Ls [split $L "\t"]

	    set gene [lindex $Ls $i_gene]
	    set synZ [lindex $Ls $i_synZ]
	    set misZ [lindex $Ls $i_misZ]
	    set pLI  [lindex $Ls $i_pLI]

	    lappend TexteToWrite "$gene\t$synZ\t$misZ\t$pLI"
	}
	WriteTextInFile [join $TexteToWrite "\n"] $GeneIntoleranceFileFormatted
	if {[catch {exec gzip $GeneIntoleranceFileFormatted} Message]} {
	    puts "-- checkGeneIntoleranceFile --"
	    puts "gzip $GeneIntoleranceFileFormatted"
	    puts "$Message\n"
	}

	file delete -force $GeneIntoleranceFileDownloaded
    }
}


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
# ExAC downloaded file: "fordist_cleaned_nonpsych_z_pli_rec_null_data.txt"
# Header:
# transcript  gene  chr  n_exons tx_start  tx_end  bp   p_syn  p_mis   p_lof   n_syn   n_mis   n_lof   adj_exp_syn  adj_exp_mis   adj_exp_lof   syn_z   mis_z  lof_z   pLI   pRecessive  pNull

## - Check if the ExAC file has been downloaded:
#    - exac-final-cnv.gene.scores071316
#
## - Check and create if necessary the following file:
#    - 'date'_ExAC.CNV-Zscore.annotations.tsv.gz
proc checkCNVintoleranceFile {} {

    global g_AnnotSV

    ## Check if the CNVintoleranceFile has been downloaded then formatted
    ######################################################################
    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Genes-based"

    set CNVintoleranceFileDownloaded [glob -nocomplain "$extannDir/ExAC/exac-final-cnv.gene.scores071316"]
    set CNVintoleranceFileFormattedGzip [glob -nocomplain "$extannDir/ExAC/*_ExAC.CNV-Zscore.annotations.tsv.gz"]

    if {$CNVintoleranceFileDownloaded eq "" && $CNVintoleranceFileFormattedGzip eq ""} {
	# No "CNV Intolerance" annotation
	return
    }

    if {[llength $CNVintoleranceFileFormattedGzip]>1} {
	puts "Several CNV-Intolerant-Genes files exist:"
	puts "$CNVintoleranceFileFormattedGzip"
	puts "Keep only one: [lindex $CNVintoleranceFileFormattedGzip end]\n"
	foreach ci [lrange $CNVintoleranceFileFormattedGzip 0 end-1] {
	    file rename -force $ci $ci.notused
	}
	return
    }

    if {$CNVintoleranceFileFormattedGzip eq ""} {
	## Create : 'date'_ExAC.CNV-Zscore.annotations.tsv.gz ; # Header: chr start end syn_z mis_z pLI
	set CNVintoleranceFileFormatted "$extannDir/ExAC/[clock format [clock seconds] -format "%Y%m%d"]_ExAC.CNV-Zscore.annotations.tsv"

	puts "...ExAC CNV Intolerant Genes configuration ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"

	puts "\t...creation of $CNVintoleranceFileFormatted.gz"
	puts "\t   (done only once during the first CNV Intolerant Genes  annotation)\n"

	set TexteToWrite {genes\tdelZ_ExAC\tdupZ_ExAC\tcnvZ_ExAC}
	foreach L [LinesFromFile $CNVintoleranceFileDownloaded] {

	    set Ls [split $L " "]

	    if {[regexp  "^gene" $Ls]} {
		set i_gene    [lsearch -regexp $Ls "^gene_symbol"]; if {$i_gene == -1} {puts "Bad syntax into $CNVintoleranceFileDownloaded.\ngene_symbol field not found - Exit with error"; exit 2}
		set i_cnvZ    [lsearch -regexp $Ls "^cnv.score"];   if {$i_cnvZ == -1} {puts "Bad syntax into $CNVintoleranceFileDownloaded.\ncnv.score field not found - Exit with error"; exit 2}
		set i_delZ    [lsearch -regexp $Ls "^del.score"];   if {$i_delZ == -1} {puts "Bad syntax into $CNVintoleranceFileDownloaded.\ndel.score field not found - Exit with error"; exit 2}
		set i_dupZ    [lsearch -regexp $Ls "^dup.score"];   if {$i_dupZ == -1} {puts "Bad syntax into $CNVintoleranceFileDownloaded.\ndup.score field not found - Exit with error"; exit 2}
		continue
	    }

	    set gene [lindex $Ls $i_gene]
	    set delZ [lindex $Ls $i_delZ]
	    set dupZ [lindex $Ls $i_dupZ]
	    set cnvZ [lindex $Ls $i_cnvZ]

	    if {![info exists del($gene)]} {
		set del($gene) $delZ
		set dup($gene) $dupZ
		set cnv($gene) $cnvZ
	    } else {
		## The file has several entries for 63 genes (several "Gencode ID" for 1 "Gene symbo")
		## We keep the highest value
		if {$delZ > $del($gene)} {set del($gene) $delZ}
		if {$dupZ > $dup($gene)} {set dup($gene) $dupZ}
		if {$cnvZ > $cnv($gene)} {set cnv($gene) $cnvZ}
	    }

	}

	foreach gene [array names del] {
	    lappend TexteToWrite "$gene\t$del($gene)\t$dup($gene)\t$cnv($gene)"
	}
	WriteTextInFile [join $TexteToWrite "\n"] $CNVintoleranceFileFormatted

	if {[catch {exec gzip $CNVintoleranceFileFormatted} Message]} {
	    puts "-- checkCNVintoleranceFile --"
	    puts "gzip $CNVintoleranceFileFormatted"
	    puts "$Message\n"
	}

	file delete -force $CNVintoleranceFileDownloaded
    }
}
