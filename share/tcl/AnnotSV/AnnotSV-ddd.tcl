############################################################################################################
# AnnotSV 2.5.2                                                                                            #
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
    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Gene-based"
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

    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Gene-based"
    set DDDfileDownloaded [glob -nocomplain "$extannDir/DDD/DDG2P.csv.gz"]

    ## Create : 'date'_DDG2P.tsv.gz
    set DDDfileFormatted "$extannDir/DDD/[clock format [clock seconds] -format "%Y%m%d"]_DDG2P.tsv"

    puts "...DDD genes configuration"

    puts "\t...creation of $DDDfileFormatted.gz"
    puts "\t   (done only once during the first DDD annotation)\n"

    set TexteToWrite {genes\tDDD_status\tDDD_mode\tDDD_consequence\tDDD_disease\tDDD_pmid}
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

	if {$category ne ""} {
	    lappend L_category($gene) "$category"
	}
	if {$allelic ne ""} {
	    lappend L_allelic($gene) "$allelic"
	}
	if {$mutation ne ""} {
	    lappend L_mutation($gene) "$mutation"
	}
	if {$disease ne ""} {
	    lappend L_disease($gene) "$disease"
	}
	if {$pmids ne ""} {
	    lappend L_pmids($gene) "$pmids"
	}
    }

    # Write outputfile (will be used as external annotation)
    set L_genes [lsort -unique $L_genes]
    foreach gene $L_genes {
	if {[info exists L_category($gene)]} {set a [join [lsort -unique $L_category($gene)] ";"]} else {set a ""}
	if {[info exists L_allelic($gene) ]} {set b [join [lsort -unique $L_allelic($gene) ] ";"]} else {set b ""}
	if {[info exists L_mutation($gene)]} {set c [join [lsort -unique $L_mutation($gene)] ";"]} else {set c ""}
	if {[info exists L_disease($gene) ]} {set d [join [lsort -unique $L_disease($gene) ] ";"]} else {set d ""}
	if {[info exists L_pmids($gene)   ]} {set e [join [lsort -unique $L_pmids($gene)   ] ";"]} else {set e ""} 
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
