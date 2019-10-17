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

## - Check if the "genemap2.txt" file exists.
## - Check and create if necessary the 'date'_OMIMannotations.tsv file.
proc checkOMIMfile {} {

    global g_AnnotSV


    set omimDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Genes-based/OMIM" 
 
    ## Check if the OMIM file has been downloaded the formatted
    ##########################################################
    set OMIMfileDownloaded [glob -nocomplain "$omimDir/genemap2.txt"]
    set OMIMfileFormattedGzip [glob -nocomplain "$omimDir/*_OMIMannotations.tsv.gz"]

    if {$OMIMfileDownloaded eq "" && $OMIMfileFormattedGzip eq ""} {
	# No OMIM annotation
	return
    } 

    if {[llength $OMIMfileFormattedGzip]>1} {
	puts "Several OMIM files exist:"
	puts "$OMIMfileFormattedGzip"
	puts "Keep only one: [lindex $OMIMfileFormattedGzip end]\n"
	foreach omim [lrange $OMIMfileFormattedGzip 0 end-1] {
	    file rename -force $omim $omim.notused
	}
	return
    } 

    if {$OMIMfileFormattedGzip eq ""} {     
	## - Create the 'date'_OMIMannotations.tsv file.
	##   Header: genes, Mim Number, Phenotypes, Inheritance

	set OMIMfileFormatted "$omimDir/[clock format [clock seconds] -format "%Y%m%d"]_OMIMannotations.tsv"
	puts "...creation of $OMIMfileFormatted.gz ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])\n"
	ReplaceTextInFile "genes\tMim Number\tPhenotypes\tInheritance" $OMIMfileFormatted
	
	# Parsing of $OMIMfileDownloaded
	foreach L [LinesFromFile $OMIMfileDownloaded] {
	    set Ls [split $L "\t"]
	    
	    if {[regexp "^# Chromosome\tGenomic" $L]} {
		set i_gene1    [lsearch -exact $Ls "Gene Symbols"];    if {$i_gene1 == -1} {puts "Bad header line syntax. Gene Symbols column not found - Exit with error"; exit 2}
		set i_gene2    [lsearch -exact $Ls "Approved Symbol"]; if {$i_gene2 == -1} {puts "Bad header line syntax. Approved Symbol column not found - Exit with error"; exit 2}
		set i_mim      [lsearch -exact $Ls "Mim Number"];      if {$i_mim == -1} {puts "Bad header line syntax. Mim Number column not found - Exit with error"; exit 2}
		set i_pheno    [lsearch -exact $Ls "Phenotypes"];      if {$i_pheno == -1} {puts "Bad header line syntax. Phenotypes column not found - Exit with error"; exit 2}
		continue
	    }
	    if {[regexp "^#" $L]} {continue}
	    
	    set gene1 [lindex $Ls $i_gene1]
	    set gene2 [lindex $Ls $i_gene2]
	    set mim [lindex $Ls $i_mim]
	    set phenoTmp [lindex $Ls $i_pheno]
	    regsub -all -nocase "Autosomal dominant" $phenoTmp "AD" phenoTmp
	    regsub -all -nocase "Autosomal recessive" $phenoTmp "AR" phenoTmp
	    regsub -all -nocase "X-linked dominant" $phenoTmp "XLD" phenoTmp
	    regsub -all -nocase "X-linked recessive" $phenoTmp "XLR" phenoTmp
	    regsub -all -nocase "Y-linked dominant" $phenoTmp "YLD" phenoTmp
	    regsub -all -nocase "Y-linked recessive" $phenoTmp "YLR" phenoTmp
	    regsub -all -nocase "X-linked" $phenoTmp "XL" phenoTmp
	    regsub -all -nocase "Y-linked" $phenoTmp "YL" phenoTmp
	    
	    set pheno {}
	    set inher {}
	    foreach p [split $phenoTmp ";"] {
		if {[regexp "(.*\\\(\[0-9\]+\\\)),? *(.*)" $p match 1pheno 1inher]} {
		    lappend pheno $1pheno
		    lappend inher $1inher
		}
	    }
	    set test [lsort -unique $pheno]
	    if {[llength $test] eq 1} {
		set pheno [lindex $pheno 0]
	    } else {
		set pheno [join $pheno "/"]
	    }
	    #if {[regexp "^/+$" $pheno]} {set pheno ""}
	    if {$pheno eq ""} {continue}
	    
	    set test [lsort -unique $inher]
	    if {[llength $test] eq 1} {
		set inher [lindex $inher 0]
	    } else {
		set inher [join $inher "/"]
	    }
	    #if {[regexp "^/+$" $inher]} {set inher ""}
	    
	    set allGenes [split $gene1 ","]
	    lappend allGenes $gene2
	    set allGenes [lsort -unique $allGenes]
	    foreach g $allGenes {
		regsub -all " " $g "" g
		if {$g eq ""} {continue}
		lappend L_genes $g
		
		lappend Mim($g) $mim
		lappend Pheno($g) $pheno
		lappend Inheritance($g) $inher
	    }
	}
	
	# creation of $OMIMfileFormatted.gz
	set L_genes [lsort -unique $L_genes]
	# Header :
	# Gene    Mim Number      Phenotypes      Inheritance
	foreach g $L_genes {
	    set Pheno($g) "[join $Pheno($g) ";"]"
	    
	    set Inheritance($g) "[join $Inheritance($g) ";"]"
	    if {[regexp "^;+$" $Inheritance($g)]} {set Inheritance($g) ""}
	    WriteTextInFile "$g\t[join $Mim($g) ";"]\t$Pheno($g)\t$Inheritance($g)" $OMIMfileFormatted
	}
	if {[catch {exec gzip $OMIMfileFormatted} Message]} {
	    puts "-- checkOMIMfile --"
	    puts "gzip $OMIMfileFormatted"
	    puts "$Message\n"
	}

	file delete -force $OMIMfileDownloaded
    }
}

## - Check if the "morbidmap.txt" file exists.
## - Check and create if necessary the 'date'_morbidgenes.tsv.gz file.
## - Check and create if necessary the 'date'_morbidgenescandidates.tsv.gz file.
proc checkMorbidGenesfile {} {

    global g_AnnotSV


    set omimDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Genes-based/OMIM" 
 
    ## Check if the MorbidGenes file has been downloaded the formatted
    ##################################################################
    set MorbidGenesFileDownloaded [glob -nocomplain "$omimDir/morbidmap.txt"]
    set MorbidGenesFileFormattedGzip [glob -nocomplain "$omimDir/*_morbidGenes.tsv.gz"]
    set MorbidGenesCandidatesFileFormattedGzip [glob -nocomplain "$omimDir/*_morbidGenescandidates.tsv.gz"]

    if {$MorbidGenesFileDownloaded eq "" && $MorbidGenesFileFormattedGzip eq ""} {
	# No MorbidGenes annotation
	return
    } 

    if {[llength $MorbidGenesFileFormattedGzip]>1} {
	puts "Several Morbid Genes files exist:"
	puts "$MorbidGenesFileFormattedGzip"
	puts "Keep only one: [lindex $MorbidGenesFileFormattedGzip end]\n"
	foreach morbid [lrange $MorbidGenesFileFormattedGzip 0 end-1] {
	    file rename -force $morbid $morbid.notused
	}
	return
    } 
    if {[llength $MorbidGenesCandidatesFileFormattedGzip]>1} {
	puts "Several Morbid Genes candidates files exist:"
	puts "$MorbidGenesCandidatesFileFormattedGzip"
	puts "Keep only one: [lindex $MorbidGenesCandidatesFileFormattedGzip end]\n"
	foreach morbid [lrange $MorbidGenesCandidatesFileFormattedGzip 0 end-1] {
	    file rename -force $morbid $morbid.notused
	}
	return
    } 

    if {$MorbidGenesFileFormattedGzip eq ""} {     
	## - Create the 'date'_morbidGenes.tsv file.
	##   Header: genes, morbidGenes
	## - Create the 'date'_morbidgenescandidates.tsv.gz file.
	##   Header: genes, morbidGenesCandidates

	set MorbidGenesFileFormatted "$omimDir/[clock format [clock seconds] -format "%Y%m%d"]_morbidGenes.tsv"
	set MorbidGenesCandidatesFileFormatted "$omimDir/[clock format [clock seconds] -format "%Y%m%d"]_morbidGenesCandidates.tsv"
	puts "...creation of $MorbidGenesFileFormatted.gz and $MorbidGenesCandidatesFileFormatted ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])\n"
	ReplaceTextInFile "genes\tmorbidGenes" $MorbidGenesFileFormatted
	ReplaceTextInFile "genes\tmorbidGenesCandidates" $MorbidGenesCandidatesFileFormatted
	
	# Parsing of $MorbidGenesFileDownloaded
	set L_morbidgenes {}
	set L_morbidgenescandidates {}
	foreach L [LinesFromFile $MorbidGenesFileDownloaded] {
	    set Ls [split $L "\t"]

	    # Header from the downloaded file:
	    # Phenotype     Gene Symbols    MIM Number      Cyto Location
	    if {[regexp "^# Phenotype\t" $L]} {
		set i_gene    [lsearch $Ls "Gene Symbols"];    if {$i_gene == -1} {puts "Bad header line syntax. Gene Symbols column not found - Exit with error"; exit 2}
		set i_pheno   [lsearch $Ls "# Phenotype"];     if {$i_pheno == -1} {puts "Bad header line syntax. Phenotype column not found - Exit with error"; exit 2}
		continue
	    }
	    if {[regexp "^#" $L]} {continue}

	    set pheno [lindex $Ls $i_pheno]

	    ## "[ ]", indicate "nondiseases", mainly genetic variations that lead to apparently abnormal laboratory test values (e.g., dysalbuminemic euthyroidal hyperthyroxinemia).
	    if {[regexp "^\\\[.+\\\]" $pheno]} {continue}  

	    ## (1) The disorder is placed on the map based on its association with a gene, but the underlying defect is not known.
	    ## (2) The disorder has been placed on the map by linkage or other statistical method; no mutation has been found.
	    if {[regexp "\\\(1\\\)|\\\(2\\\)" $pheno]} {continue}

	    set gene [lindex $Ls $i_gene]
	    set allGenes [split $gene ","]
	    set allGenes [lsort -unique $allGenes]
	    foreach g $allGenes {
		regsub -all " " $g "" g
		if {$g eq ""} {continue}
		## "{ }", indicate mutations that contribute to susceptibility to multifactorial disorders (e.g., diabetes, asthma) or to susceptibility to infection (e.g., malaria).
		## "?", before the phenotype name indicates that the relationship between the phenotype and gene is provisional. 	   
		if {[regexp "^{.+}|^\\\?" $pheno]} {
		    lappend L_morbidgenescandidates $g
		} else {
		    lappend L_morbidgenes $g
		}
	    }
	}
	
	# creation of $MorbidGenesFileFormatted.gz
	set L_morbidgenes [lsort -unique $L_morbidgenes]
	## creted header: genes, morbidGenes
	foreach g $L_morbidgenes {
	    WriteTextInFile "$g\tyes" $MorbidGenesFileFormatted
	}
	if {[catch {exec gzip $MorbidGenesFileFormatted} Message]} {
	    puts "-- checkMorbidGenesfile --"
	    puts "gzip $MorbidGenesFileFormatted"
	    puts "$Message\n"
	}

	# creation of $MorbidGenesCandidatesFileFormatted.gz
	set L_morbidgenescandidates [lsort -unique $L_morbidgenescandidates]
	## creted header: genes, morbidGenesCandidates
	foreach g $L_morbidgenescandidates {
	    WriteTextInFile "$g\tyes" $MorbidGenesCandidatesFileFormatted
	}
	if {[catch {exec gzip $MorbidGenesCandidatesFileFormatted} Message]} {
	    puts "-- checkMorbidGenesfile --"
	    puts "gzip $MorbidGenesCandidatesFileFormatted"
	    puts "$Message\n"
	}

	file delete -force $MorbidGenesFileDownloaded
    }
}
