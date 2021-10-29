############################################################################################################
# AnnotSV 3.1                                                                                              #
#                                                                                                          #
# AnnotSV: An integrated tool for Structural Variations annotation and ranking                             #
#                                                                                                          #
# Copyright (C) 2017-2021 Veronique Geoffroy (veronique.geoffroy@inserm.fr)                                #
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


    set omimDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Gene-based/OMIM" 
 
    ## Check if the OMIM file has been downloaded the formatted
    ##########################################################
    set OMIMfileDownloaded [glob -nocomplain "$omimDir/genemap2.txt"]
    set OMIMfile1FormattedGzip [glob -nocomplain "$omimDir/*_OMIM-1-annotations.tsv.gz"]
    set OMIMfile2FormattedGzip [glob -nocomplain "$omimDir/*_OMIM-2-annotations.tsv.gz"]

    if {$OMIMfileDownloaded eq ""} {
	# No OMIM annotation to do
	return
    } 

    if {[llength $OMIMfile1FormattedGzip]>1} {
	puts "Several OMIM files exist:"
	puts "$OMIMfile1FormattedGzip"
	puts "Keep only one: [lindex $OMIMfile1FormattedGzip end]\n"
	foreach omim [lrange $OMIMfile1FormattedGzip 0 end-1] {
	    file rename -force $omim $omim.notused
	}
	return
    } 
    if {[llength $OMIMfile2FormattedGzip]>1} {
	puts "Several OMIM files exist:"
	puts "$OMIMfile2FormattedGzip"
	puts "Keep only one: [lindex $OMIMfile2FormattedGzip end]\n"
	foreach omim [lrange $OMIMfile2FormattedGzip 0 end-1] {
	    file rename -force $omim $omim.notused
	}
	return
    } 

    if {$OMIMfile1FormattedGzip eq "" || $OMIMfile2FormattedGzip eq ""} {     
	## - Create the 'date'_OMIM-1-annotations.tsv and 'date'_OMIM-2-annotations.tsv files.
	##   Header1: genes, OMIM_ID
	##   Header2: genes, OMIM_phenotype, OMIM_inheritance

	set OMIMfile1Formatted "$omimDir/[clock format [clock seconds] -format "%Y%m%d"]_OMIM-1-annotations.tsv"
	set OMIMfile2Formatted "$omimDir/[clock format [clock seconds] -format "%Y%m%d"]_OMIM-2-annotations.tsv"
	puts "\t...creation of [file tail $OMIMfile1Formatted.gz] and [file tail $OMIMfile2Formatted] ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	puts "\t   (in [file dirname $OMIMfile1Formatted.gz])"
	ReplaceTextInFile "genes\tOMIM_ID" $OMIMfile1Formatted
	ReplaceTextInFile "genes\tOMIM_phenotype\tOMIM_inheritance" $OMIMfile2Formatted
	
	# Parsing of $OMIMfileDownloaded
	foreach L [LinesFromFile $OMIMfileDownloaded] {
	    set Ls [split $L "\t"]
	    
	    if {[regexp "^# Chromosome\tGenomic" $L]} {
		set i_gene1    [lsearch -exact $Ls "Gene Symbols"];         if {$i_gene1 == -1} {puts "Bad header line syntax. Gene Symbols column not found - Exit with error"; exit 2}
		set i_gene2    [lsearch -exact $Ls "Approved Gene Symbol"]; if {$i_gene2 == -1} {puts "Bad header line syntax. Approved Symbol column not found - Exit with error"; exit 2}
		set i_mim      [lsearch -exact $Ls "MIM Number"];           if {$i_mim == -1} {puts "Bad header line syntax. Mim Number column not found - Exit with error"; exit 2}
		set i_pheno    [lsearch -exact $Ls "Phenotypes"];           if {$i_pheno == -1} {puts "Bad header line syntax. Phenotypes column not found - Exit with error"; exit 2}
		continue
	    }
	    if {[regexp "^#" $L]} {continue}
	    
	    set gene1 [lindex $Ls $i_gene1]
	    set gene2 [lindex $Ls $i_gene2]
	    set allGenes [split $gene1 ","]
	    lappend allGenes $gene2
	    set allGenes [lsort -unique $allGenes]

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

	    foreach p [split $phenoTmp ";"] {
		if {[regexp " *(.*\\\(\[0-9\]+\\\)),? *(.*)" $p match 1pheno 1inher]} {
		    if {$1inher ne ""} {
			regsub -all " *, +" $1inher "," 1inher
			foreach g $allGenes {
			    regsub -all " " $g "" g
			    if {$g eq ""} {continue}	
			    lappend Pheno($g) "$1pheno $1inher"
			    lappend Inheritance($g) {*}[split $1inher ";"]
			}
		    } else {
			foreach g $allGenes {
			    regsub -all " " $g "" g
			    if {$g eq ""} {continue}	
			    lappend Pheno($g) "$1pheno"
			}
		    }
		    foreach g $allGenes {
			regsub -all " " $g "" g
			if {$g eq ""} {continue}
			lappend L_genes $g
			lappend Mim($g) $mim
		    }
		}
	    }
	}
	
	# creation of $OMIMfileFormatted.gz
	set L_genes [lsort -unique $L_genes]
	# Header :
	# Gene    OMIM_ID    
	# Gene    OMIM_phenotype      OMIM_inheritance
	foreach g $L_genes {
	    set Pheno($g) "[join [lsort -unique $Pheno($g)] ";"]"
	    if {![info exists Inheritance($g)]} {set Inheritance($g) ""}
	    set Inheritance($g) "[join [lsort -unique $Inheritance($g)] ";"]"
	    if {[regexp "^;+$" $Inheritance($g)]} {set Inheritance($g) ""}
	    WriteTextInFile "$g\t[join [lsort -unique $Mim($g)] ";"]" $OMIMfile1Formatted
	    WriteTextInFile "$g\t$Pheno($g)\t$Inheritance($g)" $OMIMfile2Formatted
	}
	if {[catch {exec gzip $OMIMfile1Formatted} Message]} {
	    puts "-- checkOMIMfile --"
	    puts "gzip $OMIMfile1Formatted"
	    puts "$Message\n"
	}
	if {[catch {exec gzip $OMIMfile2Formatted} Message]} {
	    puts "-- checkOMIMfile --"
	    puts "gzip $OMIMfile2Formatted"
	    puts "$Message\n"
	}

	file delete -force $OMIMfileDownloaded
    }
}

## - Check if the "morbidmap.txt" file exists.
## - Check and create if necessary the 'date'_morbid.tsv.gz file.
## - Check and create if necessary the 'date'_morbidCandidate.tsv.gz file.
proc checkMorbidfile {} {

    global g_AnnotSV


    set omimDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Gene-based/OMIM" 
 
    ## Check if the Morbid file has been downloaded the formatted
    ##################################################################
    set MorbidFileDownloaded [glob -nocomplain "$omimDir/morbidmap.txt"]
    set MorbidFileFormattedGzip [glob -nocomplain "$omimDir/*_morbid.tsv.gz"]
    set MorbidCandidateFileFormattedGzip [glob -nocomplain "$omimDir/*_morbidCandidate.tsv.gz"]

    if {$MorbidFileDownloaded eq "" && $MorbidFileFormattedGzip eq ""} {
	# No Morbid annotation
	return
    } 

    if {[llength $MorbidFileFormattedGzip]>1} {
	puts "Several Morbid Genes files exist:"
	puts "$MorbidFileFormattedGzip"
	puts "Keep only one: [lindex $MorbidFileFormattedGzip end]\n"
	foreach morbid [lrange $MorbidFileFormattedGzip 0 end-1] {
	    file rename -force $morbid $morbid.notused
	}
	return
    } 
    if {[llength $MorbidCandidateFileFormattedGzip]>1} {
	puts "Several Morbid Genes candidate files exist:"
	puts "$MorbidCandidateFileFormattedGzip"
	puts "Keep only one: [lindex $MorbidCandidateFileFormattedGzip end]\n"
	foreach morbid [lrange $MorbidCandidateFileFormattedGzip 0 end-1] {
	    file rename -force $morbid $morbid.notused
	}
	return
    } 

    if {$MorbidFileFormattedGzip eq ""} {     
	## - Create the 'date'_morbid.tsv file.
	##   Header: genes, OMIM_morbid
	## - Create the 'date'_morbidCandidate.tsv.gz file.
	##   Header: genes, OMIM_morbid_candidate

	set MorbidFileFormatted "$omimDir/[clock format [clock seconds] -format "%Y%m%d"]_morbid.tsv"
	set MorbidCandidateFileFormatted "$omimDir/[clock format [clock seconds] -format "%Y%m%d"]_morbidCandidate.tsv"
	puts "\t...creation of [file tail $MorbidFileFormatted.gz] and [file tail $MorbidCandidateFileFormatted] ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	puts "\t   (in [file dirname $MorbidFileFormatted.gz])"
	ReplaceTextInFile "genes\tOMIM_morbid" $MorbidFileFormatted
	ReplaceTextInFile "genes\tOMIM_morbid_candidate" $MorbidCandidateFileFormatted
	
	# Parsing of $MorbidFileDownloaded
	set L_morbid {}
	set L_morbidCandidate {}
	foreach L [LinesFromFile $MorbidFileDownloaded] {
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

	    ## The following lines are automatically either (3) or (4):
	    ## (3) indicates that the molecular basis of the disorder is known; a mutation has been found in the gene.
	    ## (4) indicates that a contiguous gene deletion or duplication syndrome, multiple genes are deleted or duplicated causing the phenotype.
	    set gene [lindex $Ls $i_gene]
	    set allGenes [split $gene ","]
	    set allGenes [lsort -unique $allGenes]
	    foreach g $allGenes {
		regsub -all " " $g "" g
		if {$g eq ""} {continue}
		## "{ }", indicate mutations that contribute to susceptibility to multifactorial disorders (e.g., diabetes, asthma) or to susceptibility to infection (e.g., malaria).
		## "?", before the phenotype name indicates that the relationship between the phenotype and gene is provisional. 	   
		if {[regexp "^{.+}|^\\\?" $pheno]} {
		    lappend L_morbidCandidate $g
		} else {
		    lappend L_morbid $g
		}
	    }
	}
	
	# creation of $MorbidFileFormatted.gz
	set L_morbid [lsort -unique $L_morbid]
	## creted header: genes, OMIM_morbid
	foreach g $L_morbid {
	    WriteTextInFile "$g\tyes" $MorbidFileFormatted
	}
	if {[catch {exec gzip $MorbidFileFormatted} Message]} {
	    puts "-- checkMorbidfile --"
	    puts "gzip $MorbidFileFormatted"
	    puts "$Message\n"
	}

	# creation of $MorbidCandidateFileFormatted.gz
	set L_morbidCandidate [lsort -unique $L_morbidCandidate]
	## creted header: genes, OMIM_morbid_candidate
	foreach g $L_morbidCandidate {
	    WriteTextInFile "$g\tyes" $MorbidCandidateFileFormatted
	}
	if {[catch {exec gzip $MorbidCandidateFileFormatted} Message]} {
	    puts "-- checkMorbidfile --"
	    puts "gzip $MorbidCandidateFileFormatted"
	    puts "$Message\n"
	}

	file delete -force $MorbidFileDownloaded
    }
}


proc fromOMIMtoPhenotype {OMIM} {
    global g_AnnotSV
    global g_phen
    global g_gene
    
    if { ! [info exists g_phen(DONE)]} {
	set OMIMfile1 [glob -nocomplain "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Gene-based/OMIM/*_OMIM-1-annotations.tsv.gz"]
	# 20201107_OMIM-1-annotations.tsv.gz 
	# Header = genes   OMIM_ID
	foreach L [LinesFromGZFile $OMIMfile1] {
	    set Ls [split $L "\t"]
	    set gene [lindex $Ls 0]
	    set L_omim [split [lindex $Ls 1] ";"]
	    foreach omim $L_omim {
		lappend g_gene($omim) $gene
	    }
	}
	
	set OMIMfile2 [glob -nocomplain "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Gene-based/OMIM/*_OMIM-2-annotations.tsv.gz"]
	# 20201107_OMIM-2-annotations.tsv.gz 
	# Header = genes   OMIM_phenotype      OMIM_inheritance	
	foreach L [LinesFromGZFile $OMIMfile2] {
	    set Ls [split $L "\t"]
	    set gene [lindex $Ls 0]
	    set phenotype [lindex $Ls 1]
	    regsub -all "\\\{|\\\}|\\\[|\\\]" $phenotype "" phenotype
	    regsub -all "/" $phenotype ";" phenotype
	    set g_phen($gene) $phenotype
	}
	
	set g_phen(DONE) 1 
    }
    
    set phenotype ""
    if {[info exists g_gene($OMIM)]} {
	foreach gene $g_gene($OMIM) {
	    if {[info exists g_phen($gene)]} {
		lappend phenotype $g_phen($gene)
	    }
	}
	set phenotype [join $phenotype ";"]
    } 

    return $phenotype
}

proc fromOMIMgeneToPhenotype {OMIMgene} {
    global g_AnnotSV
    global g_phen
    global g_gene
    
    if { ! [info exists g_phen(DONE)]} {
	fromOMIMtoPhenotype 617631 ;# use to initialize "g_phen" and "g_gene"
    }

    if {[info exists g_phen($OMIMgene)]} {
	return $g_phen($OMIMgene)
    } else {
	return ""
    }
}


proc isMorbid {gene} {
    global g_AnnotSV
    global g_morbid
    
    if {![info exists g_morbid(DONE)]} {
	set g_morbid(DONE) 1 
	set morbidFile [glob -nocomplain "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Gene-based/OMIM/*_morbid.tsv.gz"]
	if {$morbidFile ne ""} { ;# Doesn't exists for Mouse
	    set f [open "| gzip -cd $morbidFile"]
	    while {! [eof $f]} {
		set L [gets $f]
		set g_morbid([lindex $L 0]) 1
	    }
	    close $f
	}
    }
    
    if {[info exists g_morbid($gene)]} {
	return 1
    } else {
	return 0
    } 
}
