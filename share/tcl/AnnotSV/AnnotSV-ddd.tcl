############################################################################################################
# AnnotSV 2.3.3                                                                                            #
#                                                                                                          #
# AnnotSV: An integrated tool for Structural Variations annotation and ranking                             #
#                                                                                                          #
# Copyright (C) 2017-2020 Veronique Geoffroy (veronique.geoffroy@inserm.fr)                                #
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




### The Deciphering Developmental Disorders (DDD) Study: all be made publicly available via DECIPHER
#####################################
# The Development Disorder Genotype - Phenotype Database = DDG2P
#
# DDD downloaded file: "DDG2P.csv.gz"
# Header:
# "gene symbol","gene mim","disease name","disease mim","DDD category","allelic requirement","mutation consequence",phenotypes,"organ specificity list",pmids,panel,"prev symbols","hgnc id"

## - Check if the following DDD file has been downloaded:
#    - DDG2P.csv.gz
#
## - Check and create if necessary the following file:
#    - 'date'_DDG2P.sorted.tsv.gz
proc checkDDDgeneFile {} {

    global g_AnnotSV

    ## Check if the "DDD gene" file has been downloaded then formatted
    #################################################################
    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Genes-based"
    set DDDfileDownloaded [glob -nocomplain "$extannDir/DDD/DDG2P.csv.gz"]
    set DDDfileFormattedGzip [glob -nocomplain "$extannDir/DDD/*_DDG2P.tsv.gz"]

    if {$DDDfileDownloaded eq "" && $DDDfileFormattedGzip eq ""} {
	# No "DDD gene" annotation (with extAnn procedure)
	return
    }

    if {[llength $DDDfileFormattedGzip]>1} {
	puts "Several DDG2P files exist:"
	puts "$DDDfileFormattedGzip"
	puts "Keep only one: [lindex $DDDfileFormattedGzip end]\n"
	foreach ddd [lrange $DDDfileFormattedGzip 0 end-1] {
	    file rename -force $ddd $ddd.notused
	}
	return
    }

    if {$DDDfileFormattedGzip eq ""} {
	# The downloaded file exist but not the formatted.
	## Create : 'date'_DDG2P.tsv.gz
	updateDDDgeneFile
    }
}

proc updateDDDgeneFile {} {

    global g_AnnotSV

    if {[catch {package require csv} Message]} {
	puts "Tcl package \"csv\" is required for DECIPHER gene annotations."
	puts "$Message"
	puts "-> no DECIPHER gene annotations.\n"
	# No "DDD gene" annotation (with extAnn procedure)
	return
    }

    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Genes-based"
    set DDDfileDownloaded [glob -nocomplain "$extannDir/DDD/DDG2P.csv.gz"]

    ## Create : 'date'_DDG2P.tsv.gz
    set DDDfileFormatted "$extannDir/DDD/[clock format [clock seconds] -format "%Y%m%d"]_DDG2P.tsv"

    puts "...DDD genes configuration"

    puts "\t...creation of $DDDfileFormatted.gz"
    puts "\t   (done only once during the first DDD annotation)\n"

    set TexteToWrite {genes\tDDD_status\tDDD_mode\tDDD_consequence\tDDD_disease\tDDD_pmids}
    foreach L [LinesFromGZFile $DDDfileDownloaded] {
	if {[regexp  "^\"gene symbol" $L]} {
	    set Ls [::csv::split $L]
	    set i_gene      [lsearch -regexp $Ls "gene symbol"];          if {$i_gene == -1} {puts "Bad syntax into $DDDfileDownloaded.\ngene symbol field not found - Exit with error"; exit 2}
	    set i_disease   [lsearch -regexp $Ls "disease name"];         if {$i_disease == -1} {puts "Bad syntax into $DDDfileDownloaded.\ndisease name field not found - Exit with error"; exit 2}
	    set i_category  [lsearch -regexp $Ls "DDD category"];         if {$i_category == -1} {puts "Bad syntax into $DDDfileDownloaded.\nDDD category field not found - Exit with error"; exit 2}
	    set i_allelic   [lsearch -regexp $Ls "allelic requirement"];  if {$i_allelic == -1} {puts "Bad syntax into $DDDfileDownloaded.\nallelic requirement field not found - Exit with error"; exit 2}
	    set i_mutation  [lsearch -regexp $Ls "mutation consequence"]; if {$i_mutation == -1} {puts "Bad syntax into $DDDfileDownloaded.\nmutation consequence field not found - Exit with error"; exit 2}
	    set i_pmids     [lsearch -regexp $Ls "pmids"];                if {$i_pmids == -1} {puts "Bad syntax into $DDDfileDownloaded.\npmids field not found - Exit with error"; exit 2}
	    continue
	}
	set Ls [::csv::split $L]

	set gene [lindex $Ls $i_gene]
	lappend L_genes "$gene"
	set disease [lindex $Ls $i_disease]
	set category [lindex $Ls $i_category]
	set allelic [lindex $Ls $i_allelic]
	set mutation [lindex $Ls $i_mutation]
	set pmids [lindex $Ls $i_pmids]

	lappend L_category($gene) "$category"
	lappend L_allelic($gene) "$allelic"
	lappend L_mutation($gene) "$mutation"
	lappend L_disease($gene) "$disease"
	lappend L_pmids($gene) "$pmids"
    }

    # Write outputfile (will be used as external annotation)
    set L_genes [lsort -unique $L_genes]
    foreach gene $L_genes {
	set a [join $L_category($gene) "/"]; if {[regexp "^/+$" $a]} {set a ""}
	set b [join $L_allelic($gene) "/"] ; if {[regexp "^/+$" $b]} {set b ""}
	set c [join $L_mutation($gene) "/"]; if {[regexp "^/+$" $c]} {set c ""}
	set d [join $L_disease($gene) "/"] ; if {[regexp "^/+$" $d]} {set d ""}
	set e [join $L_pmids($gene) "/"]   ; if {[regexp "^/+$" $e]} {set e ""}
	lappend TexteToWrite "$gene\t$a\t$b\t$c\t$d\t$e"
    }

    WriteTextInFile [join $TexteToWrite "\n"] $DDDfileFormatted
    if {[catch {exec gzip $DDDfileFormatted} Message]} {
	puts "-- updateDDDgeneFile --"
	puts "gzip $DDDfileFormatted"
	puts "$Message\n"
    }

    file delete -force $DDDfileDownloaded
}


# DDD downloaded file: "population_cnv.txt.gz"
##############################################
# Header:
# #population_cnv_id      chr     start   end     deletion_observations   deletion_frequency      deletion_standard_error duplication_observations        duplication_frequency   duplication_standard_error      observations     frequency       standard_error  type    sample_size     study

## - Check if the following DDD files have been downloaded:
#    - population_cnv.txt.gz
#
## - Check and create if necessary the following files:
#    - 'date'_DDD_population_cnv.sorted.bed
proc checkDDDfrequencyFile {} {

    global g_AnnotSV

    ## Check if the DDD frequency file has been downloaded then formatted
    ####################################################################
    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/"

    set DDDfrequencyFileDownloaded [glob -nocomplain "$extannDir/SVincludedInFt/DDD/$g_AnnotSV(genomeBuild)/population_cnv.txt.gz"]
    set DDDfrequencyFileFormatted [glob -nocomplain "$extannDir/SVincludedInFt/DDD/$g_AnnotSV(genomeBuild)/*_DDD_population_cnv.sorted.bed"]

    if {$DDDfrequencyFileDownloaded eq "" && $DDDfrequencyFileFormatted eq ""} {
	# No "DDD freq" annotation
	set g_AnnotSV(DDDfreqAnn) 0
	return
    } else {
	set g_AnnotSV(DDDfreqAnn) 1
	# Check if the user asked for these annotations in the configfile
	set test 0
	foreach col "DDD_SV DDD_DUP_n_samples_with_SV DDD_DUP_Frequency DDD_DEL_n_samples_with_SV DDD_DEL_Frequency" {
	    if {[lsearch -exact "$g_AnnotSV(outputColHeader)" $col] ne -1} {set test 1;break}
	}
	if {$test eq 0} {set g_AnnotSV(DDDfreqAnn) 0; return}
    }

   if {[llength $DDDfrequencyFileFormatted]>1} {
	puts "Several population_cnv files exist:"
	puts "$DDDfrequencyFileFormatted"
	puts "Keep only one: [lindex $DDDfrequencyFileFormatted end]\n"
       foreach ddd [lrange $DDDfrequencyFileFormatted 0 end-1] {
	    file rename -force $ddd $ddd.notused
	}
	return
    }

    if {$DDDfrequencyFileFormatted eq ""} {
	# The downloaded file exist but not the formatted:
	##   - 'date'_DDD_population_cnv.bed ; # Header: chr start end DDD_DUP_n_samples_with_SV DDD_DUP_Frequency DDD_DEL_n_samples_with_SV DDD_DEL_Frequency
	set DDDfrequencyFileFormatted "$extannDir/SVincludedInFt/DDD/$g_AnnotSV(genomeBuild)/[clock format [clock seconds] -format "%Y%m%d"]_DDD_population_cnv.sorted.bed"

	puts "...DDD frequency configuration"

	puts "\t...creation of $DDDfrequencyFileFormatted"
	puts "\t   (done only once during the first DDD annotation)\n"

	set TexteToWrite ""
	foreach L [LinesFromGZFile $DDDfrequencyFileDownloaded] {
	    if {[regexp  "^#population_cnv_id" $L]} {
		set Ls [split $L "\t"]
		set i_chr      [lsearch -regexp $Ls "chr"];     if {$i_chr == -1} {puts "Bad syntax into $DDDfrequencyFileDownloaded.\nchr symbol field not found - Exit with error"; exit 2}
		set i_start    [lsearch -regexp $Ls "start"];   if {$i_start == -1} {puts "Bad syntax into $DDDfrequencyFileDownloaded.\nstart symbol field not found - Exit with error"; exit 2}
		set i_end      [lsearch -regexp $Ls "end"];     if {$i_end == -1} {puts "Bad syntax into $DDDfrequencyFileDownloaded.\nend symbol field not found - Exit with error"; exit 2}
		set i_delobs   [lsearch -regexp $Ls "deletion_observations"];   if {$i_delobs == -1} {puts "Bad syntax into $DDDfrequencyFileDownloaded.\ndeletion_observations symbol field not found - Exit with error"; exit 2}
		set i_delfreq  [lsearch -regexp $Ls "deletion_frequency"];      if {$i_delfreq == -1} {puts "Bad syntax into $DDDfrequencyFileDownloaded.\ndeletion_frequency symbol field not found - Exit with error"; exit 2}
		set i_dupobs  [lsearch -regexp $Ls "duplication_observations"]; if {$i_dupobs == -1} {puts "Bad syntax into $DDDfrequencyFileDownloaded.\nduplication_observations symbol field not found - Exit with error"; exit 2}
		set i_dupfreq     [lsearch -regexp $Ls "duplication_frequency"];if {$i_dupfreq == -1} {puts "Bad syntax into $DDDfrequencyFileDownloaded.\nduplication_frequency symbol field not found - Exit with error"; exit 2}
		set i_study     [lsearch -regexp $Ls "study"];     if {$i_study == -1} {puts "Bad syntax into $DDDfrequencyFileDownloaded.\nstudy symbol field not found - Exit with error"; exit 2}
		continue
	    }
	    set Ls [split $L "\t"]

	    set chr [lindex $Ls $i_chr]
	    set start [lindex $Ls $i_start]
	    set end [lindex $Ls $i_end]
	    set dupobs [lindex $Ls $i_dupobs]
	    set dupfreq [lindex $Ls $i_dupfreq]
	    set delobs [lindex $Ls $i_delobs]
	    set delfreq [lindex $Ls $i_delfreq]
	    set study [lindex $Ls $i_study]

	    if {$study eq "DDD"} {
		lappend TexteToWrite "$chr\t$start\t$end\t$dupobs\t$dupfreq\t$delobs\t$delfreq"
	    }
	}

	# Write outputfile
	WriteTextInFile [join $TexteToWrite "\n"] $DDDfrequencyFileFormatted.tmp

	# Sorting of the bedfile:
	# Intersection with very large files can cause trouble with excessive memory usage.
	# A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm.
	set sortTmpFile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_sort.tmp.bash"
	ReplaceTextInFile "#!/bin/bash" $sortTmpFile
	WriteTextInFile "# The locale specified by the environment can affects the traditional sort order. We need to use native byte values." $sortTmpFile
	WriteTextInFile "export LC_ALL=C" $sortTmpFile
	WriteTextInFile "sort -k1,1 -k2,2n $DDDfrequencyFileFormatted.tmp > $DDDfrequencyFileFormatted" $sortTmpFile
	file attributes $sortTmpFile -permissions 0755
	if {[catch {eval exec bash $sortTmpFile} Message]} {
	    puts "-- updateDDDgeneFile --"
	    puts "sort -k1,1 -k2,2n $DDDfileFormatted.tmp > $DDDfileFormatted"
	    puts "$Message"
	    puts "Exit with error"
	    exit 2
	}
	file delete -force $sortTmpFile 
	# Delete the downloaded and tmp files
	file delete -force $DDDfrequencyFileDownloaded
	file delete -force $DDDfrequencyFileFormatted.tmp
    }
}


proc DDDfrequencyAnnotation {SVchrom SVstart SVend L_i} {

    global g_AnnotSV
    global dddText


    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)"
    set DDDfrequencyFileFormatted [glob -nocomplain "$extannDir/SVincludedInFt/DDD/$g_AnnotSV(genomeBuild)/*_DDD_population_cnv.sorted.bed"]

    if {![info exists dddText(DONE)]} {

	# headerOutput "DDD_SV DDD_DUP_n_samples_with_SV DDD_DUP_Frequency DDD_DEL_n_samples_with_SV DDD_DEL_Frequency"
	set L_dddText(Empty) "{} {} {-1} {} {-1}"
	foreach i $L_i {
	    lappend dddText(Empty) "[lindex $L_dddText(Empty) $i]"
	}
	set dddText(Empty) "[join $dddText(Empty) "\t"]"

	# Intersect
	regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".intersect.ddd" tmpFile
	set tmpFile "$g_AnnotSV(outputDir)/[file tail $tmpFile]"

	file delete -force $tmpFile
	if {[catch {exec $g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(fullAndSplitBedFile) -b $DDDfrequencyFileFormatted -wa -wb > $tmpFile} Message]} {
	    puts "-- DDDfrequencyAnnotation --"
	    puts "$g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(fullAndSplitBedFile) -b $DDDfrequencyFileFormatted -wa -wb > $tmpFile"
	    puts "$Message"
	    puts "Exit with error"
	    exit 2
	}

	# Parse
	set f [open $tmpFile]
	while {![eof $f]} {
	    set L [gets $f]
	    if {$L eq ""} {continue}
	    set Ls [split $L "\t"]
	    set SVtoAnn_chrom [lindex $Ls 0]
	    set SVtoAnn_start [lindex $Ls 1]
	    set SVtoAnn_end   [lindex $Ls 2]

	    set DDDSV_start [lindex $Ls end-5]
	    set DDDSV_end   [lindex $Ls end-4]

	    set DUP_n    [lindex $Ls end-3]
	    set DUP_freq [lindex $Ls end-2]
	    set DEL_n    [lindex $Ls end-1]
	    set DEL_freq [lindex $Ls end]

	    # Select:
	    # - DDD SV that share > XX % length with the SV to annotate
	    set DDD_length [expr {$DDDSV_end-$DDDSV_start}]
	    set SVtoAnn_length [expr {$SVtoAnn_end - $SVtoAnn_start}]
	    # The DDD SV is an insertion or a breakpoint
	    if {$DDD_length<1} {
		set DDD_length 1
	    }
	    # The SV to annotate is an insertion or a breakpoint
	    if {$SVtoAnn_length<1} {
		set SVtoAnn_length 1
	    }

	    if {$SVtoAnn_start < $DDDSV_start} {
		set overlap_start $DDDSV_start
	    } else {
		set overlap_start $SVtoAnn_start
	    }
	    if {$SVtoAnn_end < $DDDSV_end} {
		set overlap_end $SVtoAnn_end
	    } else {
		set overlap_end $DDDSV_end
	    }
	    set overlap_length [expr {$overlap_end - $overlap_start}]

	    # Keeping only DDD respecting the overlaps (reciprocal or not reciprocal)
	    if {[expr {$overlap_length*100.0/$SVtoAnn_length}] < $g_AnnotSV(overlap)} {continue}
	    if {$g_AnnotSV(reciprocal) eq "yes"} {
		if {[expr {$overlap_length*100.0/$DDD_length}] < $g_AnnotSV(overlap)} {continue}
	    }

	    lappend L_SVddd($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) "$SVtoAnn_chrom:${DDDSV_start}-${DDDSV_end}"
	    lappend L_DUPn($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) "$DUP_n"
	    lappend L_DUPfreq($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) "$DUP_freq"
	    lappend L_DELn($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) "$DEL_n"
	    lappend L_DELfreq($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) "$DEL_freq"

	}


	# Loading DDD final annotation for each SV
	foreach SVtoAnn [array names L_SVddd] {
	    set L_dddText($SVtoAnn) ""
	    # Change metrics from "." to ","
	    if {[set g_AnnotSV(metrics)] eq "fr"} {
		regsub -all {\.} $L_DUPfreq($SVtoAnn) "," L_DUPfreq($SVtoAnn)
		regsub -all {\.} $L_DELfreq($SVtoAnn) "," L_DELfreq($SVtoAnn)
	    }
	    # For DUPn, DUPfreq, DELn and DELfreq: keep only the maximum value
	    lappend L_dddText($SVtoAnn) "[join $L_SVddd($SVtoAnn) ";"]"
	    set max -1
	    foreach a $L_DUPn($SVtoAnn) {if {$a > $max} {set max $a}}
	    lappend L_dddText($SVtoAnn) "$max"
	    set max -1
	    foreach a $L_DUPfreq($SVtoAnn) {if {$a > $max} {set max $a}}
	    lappend L_dddText($SVtoAnn) "$max"
	    set max -1
	    foreach a $L_DELn($SVtoAnn) {if {$a > $max} {set max $a}}
	    lappend L_dddText($SVtoAnn) "$max"
	    set max -1
	    foreach a $L_DELfreq($SVtoAnn) {if {$a > $max} {set max $a}}
	    lappend L_dddText($SVtoAnn) "$max"

	    # Keep only the user requested columns (defined in the configfile)
	    set dddText($SVtoAnn) ""
	    foreach i $L_i {
		lappend dddText($SVtoAnn) "[lindex $L_dddText($SVtoAnn) $i]"
	    }
	    set dddText($SVtoAnn) [join $dddText($SVtoAnn) "\t"]
	}
	catch {unset L_dddText}
	set dddText(DONE) 1
	file delete -force $tmpFile
    }

    if {[info exist dddText($SVchrom,$SVstart,$SVend)]} {
	return $dddText($SVchrom,$SVstart,$SVend)
    } else {
	return $dddText(Empty)
    }
}
