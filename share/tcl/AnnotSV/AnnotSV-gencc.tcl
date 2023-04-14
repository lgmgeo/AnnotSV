############################################################################################################
# AnnotSV 3.3.5                                                                                            #
#                                                                                                          #
# AnnotSV: An integrated tool for Structural Variations annotation and ranking                             #
#                                                                                                          #
# Copyright (C) 2017-2023 Veronique Geoffroy (veronique.geoffroy@inserm.fr)                                #
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




### GenCC
#####################################
# GenCC downloaded file: https://search.thegencc.org/download/action/submissions-export-tsv
#
# Header:
# "uuid","gene_curie","gene_symbol","disease_curie","disease_title","disease_original_curie","disease_original_title","classification_curie","classification_title","moi_curie","moi_title","submitter_curie","submitter_title","submitted_as_hgnc_id","submitted_as_hgnc_symbol","submitted_as_disease_id","submitted_as_disease_name","submitted_as_moi_id","submitted_as_moi_name","submitted_as_submitter_id","submitted_as_submitter_name","submitted_as_classification_id","submitted_as_classification_name","submitted_as_date","submitted_as_public_report_url","submitted_as_notes","submitted_as_pmids","submitted_as_assertion_criteria_url","submitted_as_submission_id","submitted_run_date"

## - Check if the following GenCC file has been downloaded:
#    - submissions-export-tsv
#
## - Check and create if necessary the following file:
#    - 'date'_GenCC.sorted.tsv.gz
proc checkGenCCgeneFile {} {

    global g_AnnotSV

    ## Check if the "GenCC gene" file has been downloaded then formatted
    #################################################################
    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Gene-based"
    set GenCCfileDownloaded [glob -nocomplain "$extannDir/GenCC/submissions-export-tsv"]
    set GenCCfileFormattedGzip [glob -nocomplain "$extannDir/GenCC/*_GenCC.tsv" "$extannDir/GenCC/*_GenCC.tsv.gz"]

    if {$GenCCfileDownloaded eq "" && $GenCCfileFormattedGzip eq ""} {
	# No "GenCC gene" annotation (with extAnn procedure)
	return
    }

    if {[llength $GenCCfileFormattedGzip]>1} {
	puts "Several GenCC files exist:"
	puts "$GenCCfileFormattedGzip"
	puts "Keep only one: [lindex $GenCCfileFormattedGzip end]\n"
	foreach gencc [lrange $GenCCfileFormattedGzip 0 end-1] {
	    file rename -force $gencc $gencc.notused
	}
	return
    }

    if {$GenCCfileFormattedGzip eq ""} {
	# The downloaded file exist but not the formatted.
	## Create : 'date'_GenCC.tsv.gz
	updateGenCCgeneFile
    }
}

proc updateGenCCgeneFile {} {

    global g_AnnotSV

    if {[catch {package require csv} Message]} {
	puts "Tcl package \"csv\" is required for GenCC gene annotations."
	puts "$Message"
	puts "-> no GenCC gene annotations.\n"
	# No "GenCC gene" annotation (with extAnn procedure)
	return
    }

    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Gene-based"
    set GenCCfileDownloaded [glob -nocomplain "$extannDir/GenCC/submissions-export-tsv"]

    ## Create : 'date'_GenCC.tsv.gz
    set GenCCfileFormatted "$extannDir/GenCC/[clock format [clock seconds] -format "%Y%m%d"]_GenCC.tsv"

    puts "\t...GenCC genes configuration"
    puts "\t\t...creation of $GenCCfileFormatted.gz"
    puts "\t\t   (done only once during the first GenCC annotation)"

    set TexteToWrite {genes\tGenCC_disease\tGenCC_moi\tGenCC_classification\tGenCC_pmid}
    foreach L [LinesFromFile $GenCCfileDownloaded] {
	if {[regexp  "^\"uuid" $L]} {
	    set ntext "field not found - Exit with error"
	    set Ls [::csv::split $L "\t"]
	    set i_gene    [lsearch -regexp $Ls "gene_symbol"];         if {$i_gene == -1} {puts "Bad syntax into $GenCCfileDownloaded.\ngene_symbol $ntext"; exit 2}
	    set i_disease [lsearch -regexp $Ls "disease_title"];       if {$i_disease == -1} {puts "Bad syntax into $GenCCfileDownloaded.\ndisease_title $ntext"; exit 2}
	    set i_moi     [lsearch -regexp $Ls "moi_title"];           if {$i_moi == -1} {puts "Bad syntax into $GenCCfileDownloaded.\nmoi_title $ntext"; exit 2}
	    set i_classif [lsearch -regexp $Ls "classification_title"]; if {$i_classif == -1} {puts "Bad syntax into $GenCCfileDownloaded.\nclassification_title $ntext"; exit 2}
	    set i_pmid    [lsearch -regexp $Ls "submitted_as_pmids"];  if {$i_pmid == -1} {puts "Bad syntax into $GenCCfileDownloaded.\nsubmitted_as_pmids $ntext"; exit 2}
	    continue
	}
	set Ls [::csv::split $L "\t"]

	set gene [lindex $Ls $i_gene]
	lappend L_genes "$gene"
	set disease [lindex $Ls $i_disease]
	set moi [lindex $Ls $i_moi]
	set classif [lindex $Ls $i_classif]
	set pmid [lindex $Ls $i_pmid]

	if {$disease ne ""} {
	    lappend L_disease($gene) "$disease"
	}
	if {[regexp -nocase "Autosomal dominant" $moi]} {
	    if {[regexp -nocase "with maternal imprinting" $moi]} {
		set moi "ADm"
	    } elseif {[regexp -nocase "with paternal imprinting" $moi]} {
		set moi "ADm"
	    } else {
		set moi "AD"
	    }
	} elseif {[regexp -nocase "Autosomal recessive" $moi]} {
	    set moi "AR"
	} elseif {[regexp -nocase "Digenic" $moi]} {
	    set moi "2G"
	} elseif {[regexp -nocase "Mitochondrial" $moi]} {
	    set moi "MT"
	} elseif {[regexp -nocase "Semidominant" $moi]} {
	    set moi "sD"
	} elseif {[regexp -nocase "Somatic mosaicism" $moi]} {
	    set moi "SOM"
	} elseif {[regexp -nocase "X-linked" $moi]} {
	    if {[regexp -nocase "dominant" $moi]} {
		set moi "XLD"
	    } elseif {[regexp -nocase "recessive" $moi]} {
		set moi "XLR"
	    } else {
		set moi "XL"
	    }
	} elseif {[regexp -nocase "Y-linked" $moi]} {
	    if {[regexp -nocase "dominant" $moi]} {
		set moi "YLD"
	    } elseif {[regexp -nocase "recessive" $moi]} {
		set moi "YLR"
	    } else {
		set moi "YL"
	    }
	} else {
	    set moi ""
	}
	
	if {$moi ne ""} {
	    lappend L_moi($gene) "$moi"
	    #lappend all_moi $moi ;# To use during the annotation update (and see if new moi is used by GenCC)
	    if {[regexp -nocase "incomplete penetrance|Reduced penetrance" $L]} {
		lappend L_moi($gene) "IPVE"
	    }
	}
	if {$classif ne ""} {
	    lappend L_classif($gene) "$classif"
	}
	regsub -all "\"" $pmid "" pmid
	if {$pmid ne ""} {
	    regsub -all " ," $pmid ";" pmid
	    regsub -all "," $pmid ";" pmid
	    regsub -all " " $pmid ";" pmid
	    foreach p [split "$pmid" ";"] {
		if {$p eq 0} {continue}
		if {$p eq ""} {continue}
		lappend L_pmid($gene) $p
	    }
	}
    }
    #puts [join [lsort -unique $all_moi] "\n"]; exit
    
    # Write outputfile (will be used as external annotation)
    set L_genes [lsort -unique $L_genes]
	puts "\t"
    foreach gene $L_genes {
	if {![isARefSeqGeneName $gene] && ![isAnENSEMBLgeneName $gene]} {continue}
	if {[info exists L_disease($gene)]}   {set a [join [lsort -unique $L_disease($gene)] ";"]}  else {set a ""}
	if {[info exists L_moi($gene)]}       {set b [join [lsort -unique $L_moi($gene)] ";"]}      else {set b ""}
	if {[info exists L_classif($gene) ]}  {set c [join [lsort -unique $L_classif($gene) ] ";"]} else {set c ""}
	if {[info exists L_pmid($gene)]}      {set d [join [lsort -unique $L_pmid($gene)] ";"]}     else {set d ""} 
	lappend TexteToWrite "$gene\t$a\t$b\t$c\t$d"
    }

    WriteTextInFile [join $TexteToWrite "\n"] $GenCCfileFormatted
    if {[catch {exec gzip $GenCCfileFormatted} Message]} {
	puts "-- updateGenCCgeneFile --"
	puts "gzip $GenCCfileFormatted"
	puts "$Message\n"
    }

    file delete -force $GenCCfileDownloaded
}
