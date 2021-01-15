############################################################################################################
# AnnotSV 3.0.2                                                                                            #
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

# Is it a RefSeq gene name?
proc isARefSeqGeneName {geneName} {

    global g_AnnotSV
    global g_exist 
    
    if {![info exists g_exist(RefSeq)]} {
	set g_exist(RefSeq) "1"

	set genesDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Genes/$g_AnnotSV(genomeBuild)"
	foreach L [LinesFromFile "$genesDir/genes.RefSeq.sorted.bed"] {
	    set Ls [split $L "\t"]
	    set g_exist(RefSeq,[lindex $Ls 4]) 1
	}
    }
    
    if {[info exists g_exist(RefSeq,$geneName)]} {set test 1} else {set test 0}

    return $test
}

# Is it an ENSEMBL gene name?
proc isAnENSEMBLgeneName {geneName} {

    global g_AnnotSV
    global g_exist 
    
    if {![info exists g_exist(ENSEMBL)]} {
	set g_exist(ENSEMBL) "1"

	set genesDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Genes/$g_AnnotSV(genomeBuild)"
	foreach L [LinesFromFile "$genesDir/genes.ENSEMBL.sorted.bed"] {
	    set Ls [split $L "\t"]
	    set g_exist(ENSEMBL,[lindex $Ls 4]) 1
	}
    }
    
    if {[info exists g_exist(ENSEMBL,$geneName)]} {set test 1} else {set test 0}

    return $test
}



## - Check and create if necessary the "promoter_XXXbp_*_GRCh*.sorted.bed" file.
proc checkPromoterFile {} {

    global g_AnnotSV

    ## Check if the promoters file has been formatted
    #################################################
    set g_AnnotSV(promAnn) 0
    set genesDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Genes/$g_AnnotSV(genomeBuild)"
    set regElementsDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/FtIncludedInSV/RegulatoryElements/$g_AnnotSV(genomeBuild)"

    if {![file exists $regElementsDir]} {file mkdir $regElementsDir}

    set genesFileFormatted "$genesDir/genes.$g_AnnotSV(tx).sorted.bed"
    set promoterFormatted "$regElementsDir/promoter_$g_AnnotSV(promoterSize)bp_${g_AnnotSV(tx)}_$g_AnnotSV(genomeBuild).sorted.bed"

    if {![file exists $promoterFormatted]} {
	
	# Creation of the promoterFormatted 
	puts "\t...$g_AnnotSV(tx) promoters configuration ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	puts "\t   (-promoterSize = $g_AnnotSV(promoterSize))\n"
	foreach L [LinesFromFile $genesFileFormatted] {
	    # Line example:
	    # 1       17368   17436   -       MIR6859-2       NR_107062       17436   17436   17368,  17436,
	    set Ls [split $L "\t"]	    
	    set chrom [lindex $Ls 0]
	    set strand [lindex $L 3]
	    if {$strand eq "+"} {
		set promRight [lindex $Ls 1]
		set promLeft [expr $promRight-$g_AnnotSV(promoterSize)]
		if {$promLeft < 0} {set promLeft 0}
	    } elseif {$strand eq "-"} {
		set promLeft [lindex $Ls 2]
		set promRight [expr $promLeft+$g_AnnotSV(promoterSize)]
	    } else {
		puts "strand eq \"$strand\" -- $L"
	    }
	    lappend genesProm($chrom\t$promLeft\t$promRight) [lindex $Ls 4]
	    lappend L_prom($chrom) "$promLeft\t$promRight"
	}
	set L_allProm {}
	foreach chrom {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M MT} {
	    if {![info exists L_prom($chrom)]} {continue}
	    set L_prom($chrom) [lsort -command AscendingSortOnElement0 [lsort -command AscendingSortOnElement1 [lsort -unique $L_prom($chrom)]]]
	    foreach prom $L_prom($chrom) {
		lappend L_allProm "$chrom\t$prom\t${g_AnnotSV(tx)}_promoter\t[join [lsort -unique $genesProm($chrom\t$prom)] ";"]"
	    }
	}
	WriteTextInFile [join $L_allProm "\n"] $promoterFormatted.tmp
	# Sorting of the bedfile:
	# Intersection with very large files can cause trouble with excessive memory usage.
	# A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm. 
	set sortTmpFile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_sort.tmp.bash"
	ReplaceTextInFile "#!/bin/bash" $sortTmpFile
	WriteTextInFile "# The locale specified by the environment can affects the traditional sort order. We need to use native byte values." $sortTmpFile
	WriteTextInFile "export LC_ALL=C" $sortTmpFile
	WriteTextInFile "sort -k1,1 -k2,2n $promoterFormatted.tmp > $promoterFormatted" $sortTmpFile
	file attributes $sortTmpFile -permissions 0755
	if {[catch {eval exec bash $sortTmpFile} Message]} {
	    puts "-- checkPromoterFile --"
	    puts "sort -k1,1 -k2,2n $promoterFormatted.tmp > $promoterFormatted"
	    puts "$Message"
	    puts "Exit with error"
	    exit 2
	}
	file delete -force $sortTmpFile 
	file delete -force $promoterFormatted.tmp
    }
    set g_AnnotSV(promAnn) 1
}


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
proc checkEAfiles {} {

    global g_AnnotSV


    ## Check if EA files have been formatted
    ########################################
    set g_AnnotSV(EAann) 0
    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)"
    set regElementsDir "$extannDir/FtIncludedInSV/RegulatoryElements"
    set EARefSeqFileFormattedGRCh37 "$regElementsDir/$g_AnnotSV(genomeBuild)/EA_RefSeq_$g_AnnotSV(genomeBuild).sorted.bed" ;# GRCh38 version should be manually created by lift over
    set EAENSEMBLfileFormattedGRCh37 "$regElementsDir/$g_AnnotSV(genomeBuild)/EA_ENSEMBL_$g_AnnotSV(genomeBuild).sorted.bed" ;# GRCh38 version should be manually created by lift over

    set necessaryEAfile "$regElementsDir/$g_AnnotSV(genomeBuild)/EA_${g_AnnotSV(tx)}_$g_AnnotSV(genomeBuild).sorted.bed"
    if {[file exists $necessaryEAfile]} {set g_AnnotSV(EAann) 1}

    if {[file exists $EARefSeqFileFormattedGRCh37] && [file exists $EAENSEMBLfileFormattedGRCh37]} {
	return
    }

    ###############################################
    ## No formatted EA files available at this step
    ###############################################

    ## Check if some EA files have been downloaded
    ##############################################
    set L_EAfiles [glob -nocomplain "$regElementsDir/GRCh37/*_EP.txt"] ;# only distributed in GRCh37
    if {$L_EAfiles eq ""} {
	return
    }
    
    ## Downloaded EA files are available
    ####################################

    # Line example: 
    # chr1:1252950-1254670_ENSG00000160087$UBE2J2$chr1$1209265$-      1.205829
    #
    # Format:
    # For human, the data format in the file listed as "A:B-C_D$E$F$G$H	I":
    # A - Chromosome of enhancer.
    # B - The starting position of enhancer.
    # C - The ending position of enhancer.
    # D - Gene ensembl ID.
    # E - Gene ensembl Name.
    # F - chromosome of gene.
    # G - the position of the gene transcription start site.
    # H - The strand of DNA the gene located.
    # I - The predition score of the enhancer-gene interaction.

    puts "\t...EnhancerAtlas configuration ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    
    ## Creation of EA_*_GRCh37.sorted.bed
    ##########################################
    ##   - # Header: chr EAstart EAend EAtype EAgene

    puts "\t\t...creation of [file tail $EARefSeqFileFormattedGRCh37] and [file tail $EAENSEMBLfileFormattedGRCh37]"
    puts "\t\t   (done only once)"
    puts "\t\t   For GRCh38 use, please lift over this GRCh37 file to GRCh38"

    file delete -force $EARefSeqFileFormattedGRCh37
    file delete -force $EAENSEMBLfileFormattedGRCh37
    
    # Writings:
    file delete -force $EARefSeqFileFormattedGRCh37.tmp
    file delete -force $EAENSEMBLfileFormattedGRCh37.tmp
    foreach EAfile $L_EAfiles {
	foreach L [LinesFromFile $EAfile] {
	    if {[regexp "chr(\[0-9XYMT\]+):(\[0-9\]+)-(\[0-9\]+)_\[^\\\$\]+\\\$(\[^\\\$\]+)" $L match chrom start end gene]} {
		# RefSeq
		if {[isARefSeqGeneName $gene]} {
		    lappend L_coordRefSeq($chrom) "$chrom\t$start\t$end"
		    lappend L_genesRefSeq($chrom\t$start\t$end) "$gene"
		}
		# ENSEMBL
		if {[isAnENSEMBLgeneName $gene]} {
		    lappend L_coordENSEMBL($chrom) "$chrom\t$start\t$end"
		    lappend L_genesENSEMBL($chrom\t$start\t$end) "$gene"
		}		
	    } 
	}
    }		   
    foreach chrom {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M MT} {
	# RefSeq
	if {[info exists L_coordRefSeq($chrom)]} {
	    set L_coordRefSeq($chrom) [lsort -command AscendingSortOnElement1 [lsort -command AscendingSortOnElement2 [lsort -unique $L_coordRefSeq($chrom)]]]
	    set L_linesToWriteGRCh37 {}
	    foreach coord $L_coordRefSeq($chrom) {
		set L_genesRefSeq($coord) [lsort -unique $L_genesRefSeq($coord)]
		lappend L_linesToWriteGRCh37 "$coord\tEA_enhancer\t[join $L_genesRefSeq($coord) ";"]"
	    }
	    WriteTextInFile "[join $L_linesToWriteGRCh37 "\n"]" $EARefSeqFileFormattedGRCh37.tmp
	}
	# ENSEMBL
	if {[info exists L_coordENSEMBL($chrom)]} {
	    set L_coordENSEMBL($chrom) [lsort -command AscendingSortOnElement1 [lsort -command AscendingSortOnElement2 [lsort -unique $L_coordENSEMBL($chrom)]]]
	    set L_linesToWriteGRCh37 {}
	    foreach coord $L_coordENSEMBL($chrom) {
		set L_genesENSEMBL($coord) [lsort -unique $L_genesENSEMBL($coord)]
		lappend L_linesToWriteGRCh37 "$coord\tEA_enhancer\t[join $L_genesENSEMBL($coord) ";"]"
	    }
	    WriteTextInFile "[join $L_linesToWriteGRCh37 "\n"]" $EAENSEMBLfileFormattedGRCh37.tmp
	}
    }
        
    # Sorting of the bedfile:
    # Intersection with very large files can cause trouble with excessive memory usage.
    # A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm. 
    # RefSeq
    set sortTmpFile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_sort.tmp.bash"
    ReplaceTextInFile "#!/bin/bash" $sortTmpFile
    WriteTextInFile "# The locale specified by the environment can affects the traditional sort order. We need to use native byte values." $sortTmpFile
    WriteTextInFile "export LC_ALL=C" $sortTmpFile
    WriteTextInFile "sort -k1,1 -k2,2n $EARefSeqFileFormattedGRCh37.tmp >> $EARefSeqFileFormattedGRCh37" $sortTmpFile
    file attributes $sortTmpFile -permissions 0755
    if {[catch {eval exec bash $sortTmpFile} Message]} {
	puts "-- checkEAfiles --"
	puts "sort -k1,1 -k2,2n $EARefSeqFileFormattedGRCh37.tmp >> $EARefSeqFileFormattedGRCh37"
	puts "$Message"
	puts "Exit with error"
	exit 2
    }
    file delete -force $sortTmpFile 
    file delete -force $EARefSeqFileFormattedGRCh37.tmp
    # ENSEMBL
    set sortTmpFile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_sort.tmp.bash"
    ReplaceTextInFile "#!/bin/bash" $sortTmpFile
    WriteTextInFile "# The locale specified by the environment can affects the traditional sort order. We need to use native byte values." $sortTmpFile
    WriteTextInFile "export LC_ALL=C" $sortTmpFile
    WriteTextInFile "sort -k1,1 -k2,2n $EAENSEMBLfileFormattedGRCh37.tmp >> $EAENSEMBLfileFormattedGRCh37" $sortTmpFile
    file attributes $sortTmpFile -permissions 0755
    if {[catch {eval exec bash $sortTmpFile} Message]} {
	puts "-- checkEAfiles --"
	puts "sort -k1,1 -k2,2n $EAENSEMBLfileFormattedGRCh37.tmp >> $EAENSEMBLfileFormattedGRCh37"
	puts "$Message"
	puts "Exit with error"
	exit 2
    }
    file delete -force $sortTmpFile 
    file delete -force $EAENSEMBLfileFormattedGRCh37.tmp

    # Delete the downloaded files
    if {[file exists $EARefSeqFileFormattedGRCh37] && [file exists $EAENSEMBLfileFormattedGRCh37]} {
	foreach EAfile $L_EAfiles {
	    file delete -force $EAfile
	}
    }
    
    set necessaryEAfile "$regElementsDir/$g_AnnotSV(genomeBuild)/EA_${g_AnnotSV(tx)}_$g_AnnotSV(genomeBuild).sorted.bed"
    if {[file exists $necessaryEAfile]} {set g_AnnotSV(EAann) 1}
}




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
proc checkGHfiles {} {

    global g_AnnotSV


    ## Check if GH files has been formatted
    #######################################
    set g_AnnotSV(GHAnn) 0
    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)"
    set regElementsDir "$extannDir/FtIncludedInSV/RegulatoryElements"
    set GHRefSeqFileFormatted  "$regElementsDir/$g_AnnotSV(genomeBuild)/GH_RefSeq_$g_AnnotSV(genomeBuild).sorted.bed"
    set GHENSEMBLfileFormatted "$regElementsDir/$g_AnnotSV(genomeBuild)/GH_ENSEMBL_$g_AnnotSV(genomeBuild).sorted.bed"

    if {[file exists $GHRefSeqFileFormatted] && [file exists $GHENSEMBLfileFormatted]} {
	set g_AnnotSV(GHAnn) 1
	return
    }

    ###############################################
    ## No formatted GH files available at this step
    ###############################################

    ## Check if GH file has been downloaded and unzipped
    ####################################################
    # Checks if the unzipped GH files exist:
    set GHelementsF "$regElementsDir/GeneHancer_elements.txt" ;# GRCh38
    set GHassociationsF "$regElementsDir/GeneHancer_gene_associations_scores.txt" ;# genes
    set GHhg19F "$regElementsDir/GeneHancer_hg19.txt" ;# GRCh37
    foreach GHfile "$GHelementsF $GHassociationsF $GHhg19F" {
	if {![file exists $GHfile]} {
	    if {$g_AnnotSV(organism) eq "Human"} {
		puts "\nWARNING: No GeneHancer annotations available."
		puts "(Please, see in the README file how to add these annotations. Users need to contact the GeneCards team.)\n"
	    }
	    return
	}
    }
    
    ## Downloaded GH files have been unzipped and are available
    ###########################################################
    puts "\t...GeneHancer configuration ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    
    ## Creation of the 4 GH formatted files (RefSeq and ENSEMBL, GRCh37 and GRCh38)
    ##   - # Header: chr GHstart GHend GHtype GHgene
    set GHRefSeqFileFormatted37  "$regElementsDir/GRCh37/GH_RefSeq_GRCh37.sorted.bed"
    set GHRefSeqFileFormatted38  "$regElementsDir/GRCh38/GH_RefSeq_GRCh38.sorted.bed"
    set GHENSEMBLfileFormatted37 "$regElementsDir/GRCh37/GH_ENSEMBL_GRCh37.sorted.bed"
    set GHENSEMBLfileFormatted38 "$regElementsDir/GRCh38/GH_ENSEMBL_GRCh38.sorted.bed"

    puts "\t\t...creation of the 4 GH_*_GRCh3*.sorted.bed files"
    puts "\t\t   (done only once)"

    file delete -force $GHRefSeqFileFormatted37
    file delete -force $GHRefSeqFileFormatted38
    file delete -force $GHENSEMBLfileFormatted37
    file delete -force $GHENSEMBLfileFormatted38
    
    # $GHelementsF (GRCh38)
    # Header: chr element_start element_end GHid is_elite regulatory_element_type
    foreach L [LinesFromFile $GHelementsF] {
	set Ls [split $L "\t"]
	if {[regexp "^chr" $L]} {
	    set i_GHstart [lsearch -exact $Ls "element_start"]
	    set i_GHend [lsearch -exact $Ls "element_end"]
	    set i_GHid [lsearch -exact $Ls "GHid"]
	    set i_GHtype [lsearch -exact $Ls "regulatory_element_type"]
	    continue
	}
	regsub "chr" [lindex $Ls 0] "" chrom
	set GHid [lindex $Ls $i_GHid]
	lappend L_GHid($chrom) "$GHid"

	set GHstart [lindex $Ls $i_GHstart]
	set GHend [lindex $Ls $i_GHend]
	set coordGRCh38($GHid) "$chrom\t$GHstart\t$GHend"
	set type($GHid) "GH_[lindex $Ls $i_GHtype]"
	set L_GHgene(RefSeq,$GHid) "" ;# initialisation before the next foreach
	set L_GHgene(ENSEMBL,$GHid) "" ;# initialisation before the next foreach
	
    }

    # $GHassociationsF => genes
    # Header: GHid symbol combined_score is_elite
    foreach L [LinesFromFile $GHassociationsF] {
	set Ls [split $L "\t"]
	if {[regexp "^GHid" $L]} {
	    set i_GHid [lsearch -exact $Ls "GHid"]
	    set i_regulatedGene [lsearch -exact $Ls "symbol"]
	    continue
	}
	set GHid [lindex $Ls $i_GHid]
	set gene [lindex $Ls $i_regulatedGene]
	if {[isARefSeqGeneName $gene]} {
	    lappend L_GHgene(RefSeq,$GHid) $gene
	}
	if {[isAnENSEMBLgeneName $gene]} {
	    lappend L_GHgene(ENSEMBL,$GHid) $gene
	}
    }

    # $GHhg19F
    # Header: chromosome start end GHid
    foreach L [LinesFromFile $GHhg19F] {
	set Ls [split $L "\t"]
	if {[regexp "^chromosome" $L]} {
	    set i_GHstart [lsearch -exact $Ls "start"]
	    set i_GHend [lsearch -exact $Ls "end"]
	    set i_GHid [lsearch -exact $Ls "GHid"]
	    continue
	}
	set GHid [lindex $Ls $i_GHid]
	regsub "chr" [lindex $Ls 0] "" chrom
	set GHstart [lindex $Ls $i_GHstart]
	set GHend [lindex $Ls $i_GHend]
	set coordGRCh37($GHid) "$chrom\t$GHstart\t$GHend"
    }

    # Writings:
    file delete -force $GHRefSeqFileFormatted37.tmp
    file delete -force $GHRefSeqFileFormatted38.tmp
    file delete -force $GHENSEMBLfileFormatted37.tmp
    file delete -force $GHENSEMBLfileFormatted38.tmp
    foreach chrom {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M MT} {
	if {![info exists L_GHid($chrom)]} {continue}
	set L_linesToWriteRefSeqGRCh37 {}
	set L_linesToWriteRefSeqGRCh38 {}
	set L_linesToWriteENSEMBLGRCh37 {}
	set L_linesToWriteENSEMBLGRCh38 {}
	foreach GHid $L_GHid($chrom) {
	    # remove some enhancers (those for which no external identifiers are available, e.g. "hsa-miR-1273e-082")
	    ## Header: chr GHstart GHend GHtype GHgene
	    if {$L_GHgene(RefSeq,$GHid) ne ""} {
		catch {lappend L_linesToWriteRefSeqGRCh37 "$coordGRCh37($GHid)\t$type($GHid)\t[join $L_GHgene(RefSeq,$GHid) ";"]"}
		catch {lappend L_linesToWriteRefSeqGRCh38 "$coordGRCh38($GHid)\t$type($GHid)\t[join $L_GHgene(RefSeq,$GHid) ";"]"}
	    }
	    if {$L_GHgene(ENSEMBL,$GHid) ne ""} {
		catch {lappend L_linesToWriteENSEMBLGRCh37 "$coordGRCh37($GHid)\t$type($GHid)\t[join $L_GHgene(ENSEMBL,$GHid) ";"]"}
		catch {lappend L_linesToWriteENSEMBLGRCh38 "$coordGRCh38($GHid)\t$type($GHid)\t[join $L_GHgene(ENSEMBL,$GHid) ";"]"}
	    }
	}
	set L_linesToWriteRefSeqGRCh37 [lsort -command AscendingSortOnElement1 [lsort -command AscendingSortOnElement2 $L_linesToWriteRefSeqGRCh37]]
	set L_linesToWriteRefSeqGRCh38 [lsort -command AscendingSortOnElement1 [lsort -command AscendingSortOnElement2 $L_linesToWriteRefSeqGRCh38]]
	set L_linesToWriteENSEMBLGRCh37 [lsort -command AscendingSortOnElement1 [lsort -command AscendingSortOnElement2 $L_linesToWriteENSEMBLGRCh37]]
	set L_linesToWriteENSEMBLGRCh38 [lsort -command AscendingSortOnElement1 [lsort -command AscendingSortOnElement2 $L_linesToWriteENSEMBLGRCh38]]
	WriteTextInFile "[join $L_linesToWriteRefSeqGRCh37 "\n"]" $GHRefSeqFileFormatted37.tmp
	WriteTextInFile "[join $L_linesToWriteRefSeqGRCh38 "\n"]" $GHRefSeqFileFormatted38.tmp
	WriteTextInFile "[join $L_linesToWriteENSEMBLGRCh37 "\n"]" $GHENSEMBLfileFormatted37.tmp
	WriteTextInFile "[join $L_linesToWriteENSEMBLGRCh38 "\n"]" $GHENSEMBLfileFormatted38.tmp
    }
    
    # Sorting of the bedfile:
    # Intersection with very large files can cause trouble with excessive memory usage.
    # A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm. 
    foreach GHfile {GHRefSeqFileFormatted37 GHRefSeqFileFormatted38 GHENSEMBLfileFormatted37 GHENSEMBLfileFormatted38} {
	set sortTmpFile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_sort.tmp.bash"
	ReplaceTextInFile "#!/bin/bash" $sortTmpFile
	WriteTextInFile "# The locale specified by the environment can affects the traditional sort order. We need to use native byte values." $sortTmpFile
	WriteTextInFile "export LC_ALL=C" $sortTmpFile
	WriteTextInFile "sort -k1,1 -k2,2n [set $GHfile].tmp >> [set $GHfile]" $sortTmpFile
	file attributes $sortTmpFile -permissions 0755
	if {[catch {eval exec bash $sortTmpFile} Message]} {
	    puts "-- checkGHfiles --"
	    puts "sort -k1,1 -k2,2n [set $GHfile].tmp >> [set $GHfile]"
	    puts "$Message"
	    puts "Exit with error"
	    exit 2
	}
	file delete -force $sortTmpFile 
	file delete -force [set $GHfile].tmp
    }

    # Delete the downloaded files
    if {[file exists $GHRefSeqFileFormatted37] && [file exists $GHRefSeqFileFormatted38] && [file exists $GHENSEMBLfileFormatted37] && [file exists $GHENSEMBLfileFormatted38]} {
	foreach GHfile "$GHelementsF $GHassociationsF $GHhg19F $regElementsDir/GeneHancer_tissues.txt $regElementsDir/ReadMe.txt" {
	    file delete -force $GHfile
	}
    } 
    set g_AnnotSV(GHAnn) 1
}





## Annotate the SV bedFile with the regulatory elements files
## Creation of the g_re global variable
proc regulatoryElementsAnnotation {L_allGenesOverlapped} {

    global g_AnnotSV
    global g_re
    global g_HITS
    
    # Intersect of the input SV bedfile with the regulatory elements files
    ######################################################################
    ## Sorted SV bed file: $g_AnnotSV(bedFile)
    regsub ".annotated.tsv$" $g_AnnotSV(outputDir)/$g_AnnotSV(outputFile) ".SV_RE_intersect.tmp" SV_RE_intersectBEDfile
    file delete -force $SV_RE_intersectBEDfile

    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)"
    set L_REfiles {}
    lappend L_REfiles "$extannDir/FtIncludedInSV/RegulatoryElements/$g_AnnotSV(genomeBuild)/EA_${g_AnnotSV(tx)}_$g_AnnotSV(genomeBuild).sorted.bed"     ;#EnhancerAtlas
    lappend L_REfiles "$extannDir/FtIncludedInSV/RegulatoryElements/$g_AnnotSV(genomeBuild)/GH_${g_AnnotSV(tx)}_$g_AnnotSV(genomeBuild).sorted.bed"     ;#GeneHancer
    lappend L_REfiles "$extannDir/FtIncludedInSV/RegulatoryElements/$g_AnnotSV(genomeBuild)/promoter_${g_AnnotSV(promoterSize)}bp_${g_AnnotSV(tx)}_$g_AnnotSV(genomeBuild).sorted.bed" ;# promoter
    foreach reFile $L_REfiles {
	if {![file exists $reFile]} {continue}
	if {[catch {exec $g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(bedFile) -b $reFile -wa -wb >> $SV_RE_intersectBEDfile} Message]} {
	    if {[catch {exec $g_AnnotSV(bedtools) intersect -a $g_AnnotSV(bedFile) -b $reFile -wa -wb >> $SV_RE_intersectBEDfile} Message]} {
		puts "-- regulatoryElementsAnnotation --"
		puts "$g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(bedFile) -b $reFile -wa -wb >> $SV_RE_intersectBEDfile"
		puts "$Message"
		puts "Exit with error"
		exit 2
	    }
	}
    }

    # Definition of the g_re($SVfromBED) variable (list all the regulated genes impacted by 1 SV)
    set L_allRegulatedGenes {} ;# used for exomiser
    if {![file exists $SV_RE_intersectBEDfile ] || [file size $SV_RE_intersectBEDfile] eq 0} {
	# No intersection between SV and regulatory elements
	## Delete temporary file
	file delete -force $SV_RE_intersectBEDfile
    } else {
	set f [open "$SV_RE_intersectBEDfile"]
	while {! [eof $f]} {
	    set L [gets $f]
	    set Ls [split $L "\t"]
	    # SVfromBED ("chrom, start, end")
	    set SVfromBED "[join [lrange $Ls 0 2] "\t"]"
	    lappend g_re($SVfromBED) {*}[split [lindex $Ls end] ";"] ;# regulated genes
	    lappend L_allRegulatedGenes {*}[split [lindex $Ls end] ";"]
	}
	close $f
	foreach sv [array names g_re] {
	    set g_re($sv) [lsort -unique $g_re($sv)]
	}
	## Delete temporary file
	if {$g_AnnotSV(REreport) eq "no"} {
	    file delete -force $SV_RE_intersectBEDfile
	}
    }

    ## Display
    puts "...searching for SV overlaps with a gene or a regulatory elements"
    puts "\t...[llength $L_allGenesOverlapped] genes overlapped with an SV"
    set L_allRegulatedGenes [lsort -unique $L_allRegulatedGenes]
    puts "\t...[llength $L_allRegulatedGenes] genes regulated by a regulatory element which is overlapped with an SV\n"
    
    ## Preparation of the phenotype-driven analysis (Exomiser)
    set L_allGenes $L_allGenesOverlapped
    lappend L_allGenes {*}$L_allRegulatedGenes
    set L_allGenes [lsort -unique $L_allGenes]
    if {$g_AnnotSV(hpo) ne "" && $L_allGenes ne ""} {
	runExomiser "$L_allGenes" "$g_AnnotSV(hpo)" 
    }
    ## HI/TS information for these regulated genes
    set clingenDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Gene-based/ClinGen"
    set ClinGenFileFormattedGzip [lindex [glob -nocomplain "$clingenDir/*_ClinGenAnnotations.tsv.gz"] end]
    foreach g $L_allRegulatedGenes {
	set t($g) 1
    }
    if {$ClinGenFileFormattedGzip ne ""} {
	set f [open "| gzip -cd $ClinGenFileFormattedGzip"]
	while {! [eof $f]} {
	    set L [gets $f]
	    set Ls [split $L "\t"]
	    set gene [lindex $Ls 0]
	    if {[info exists t($gene)]} {
		set L_ann ""
		set HI [lindex $Ls 1]; if {$HI ne "Not yet evaluated"} {lappend L_ann "HI=$HI"}
		set TS [lindex $Ls 2]; if {$TS ne "Not yet evaluated"} {lappend L_ann "TS=$TS"}
		if {$L_ann ne ""} {
		    set g_HITS($gene) "[join $L_ann "/"]" ;# used as sub-annotation (between parentheses) => slashes separated
		}
	    }
	}
    }
    
    return
	
    
##     
##    # Complete the "fullAndSplitBedFile" ($g_AnnotSV(outputDir)/$g_AnnotSV(outputFile).tmp)
##    #######################################################################################
##    # Used for the insertion of the "full/split" information
##    set L_Bed [LinesFromFile $g_AnnotSV(bedFile)]
##    set L "[FirstLineFromFile $SV_RE_intersectBEDfile]"
##    set Ls [split $L "\t"]
##
##    set previousSV ""
##    set L_LinesToWrite {}
##    foreach L [LinesFromFile $fullAndSplitBedFile] {
##	set Ls [split $L "\t"]
##
##	# Annotation_mode (full or split)
##	set AnnotationMode [lindex $Ls end]
##	
##	# To write
##	if {$AnnotationMode eq "full"} {
##	    if {$g_AnnotSV(svtBEDcol) ne "-1"} {
##		set SVtype "\t[lindex $Ls "$g_AnnotSV(svtBEDcol)"]"
##	    } else {
##		set SVtype ""
##	    }
##	    set currentSV "[join [lrange $Ls 0 2] "\t"]$SVtype"
##
##	    if {[info exists g_re($previousSV)]} {
##		foreach l $g_re($previousSV) {
##		    lappend L_LinesToWrite "$previousFullLine\t[join $l "\t"]\tsplit-re"
##		}
##	    }
##	    set previousSV $currentSV
##	    set previousFullLine $L
##	}
##	lappend L_LinesToWrite $L
##    }
##
##    ReplaceTextInFile "[join $L_LinesToWrite "\n"]" $fullAndSplitBedFile
##    ## Delete temporary file
##    file delete -force $SV_RE_intersectBEDfile
##

}
