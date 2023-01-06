############################################################################################################
# AnnotSV 3.2.2                                                                                            #
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
	set ENSEMBLfile "$genesDir/genes.ENSEMBL.sorted.bed"
	if {[file exists $ENSEMBLfile]} {
	    foreach L [LinesFromFile $ENSEMBLfile] {
		set Ls [split $L "\t"]
		set g_exist(ENSEMBL,[lindex $Ls 4]) 1
	    }
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
	regsub -all "/" [lindex $Ls $i_GHtype] "_" type($GHid)
	regsub "promoter_enhancer" $type($GHid) "enhancer" type($GHid) ;# to simplify the annotation in the AnnotSV output \
	    # (we don't want: "RE=GH_enhancer+GH_promoter+GH_promoter_enhancer" in the RE_gene output column)
	set type($GHid) "GH_[string tolower $type($GHid)]"
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
proc checkMiRTargetLinkFiles {} {

    global g_AnnotSV

    set g_AnnotSV(miRNAann) 0

    ## Check if miRTargetLink files has been formatted
    ##################################################
    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)"
    set regElementsDir "$extannDir/FtIncludedInSV/RegulatoryElements"
    set miRNARefSeqFileFormatted  "$regElementsDir/$g_AnnotSV(genomeBuild)/miRTargetLink_RefSeq_$g_AnnotSV(genomeBuild).sorted.bed"
    set miRNAENSEMBLfileFormatted "$regElementsDir/$g_AnnotSV(genomeBuild)/miRTargetLink_ENSEMBL_$g_AnnotSV(genomeBuild).sorted.bed"
    if {$g_AnnotSV(organism) eq "Human"} {
	if {[file exists $miRNARefSeqFileFormatted] && [file exists $miRNAENSEMBLfileFormatted]} {
	    set g_AnnotSV(miRNAann) 1
	    return
	}
    } elseif {$g_AnnotSV(organism) eq "Mouse"} {
	if {[file exists $miRNARefSeqFileFormatted]} {
	    set g_AnnotSV(miRNAann) 1
	    return
	}
    }

    ##########################################################
    ## No formatted miRTargetLink files available at this step
    ##########################################################

    ## miRTargetLink is only available in GRCh38 and mm10
    #####################################################
    if {$g_AnnotSV(genomeBuild) ne "GRCh38" && $g_AnnotSV(genomeBuild) ne "mm10"} {
	return
    }
    
    ## Check if miRTargetLink file has been downloaded
    ##################################################
    # Checks if the miRTargetLink file exist:
    # (we only keep the validated miRNA file. But the code is OK to also keep the predicted miRNA file if wanted)
    if {$g_AnnotSV(genomeBuild) eq "GRCh38"} {
	set miRNAfilesDownloaded [glob -nocomplain "$regElementsDir/*_miRNA-gene_pairs_hsa_miRBase_*_GRCh38_location_augmented.tsv"];# Validated (+ Predicted) miRNA files
    } elseif {$g_AnnotSV(genomeBuild) eq "mm10"} {
	set miRNAfilesDownloaded [glob -nocomplain "$regElementsDir/*_miRNA-gene_pairs_mmu_miRBase_*_GRCm38_location_augmented.tsv"];# Validated (+ Predicted) miRNA files
    }
    if {$miRNAfilesDownloaded eq ""} {
	return
    }

    ## Downloaded miRTargetLink file is available
    #############################################
    puts "\t...miRTargetLink configuration ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    
    ## Creation of the 2 miRNA formatted files: RefSeq (GRCh38 + mm10)and ENSEMBL (GRCh38)
    ##   - # Header: chr miRNAstart miRNAend miRNAtype miRNAgene
    set miRNARefSeqFileFormatted38or10  "$regElementsDir/$g_AnnotSV(genomeBuild)/miRTargetLink_RefSeq_$g_AnnotSV(genomeBuild).sorted.bed"
    set miRNAENSEMBLfileFormatted38     "$regElementsDir/$g_AnnotSV(genomeBuild)/miRTargetLink_ENSEMBL_$g_AnnotSV(genomeBuild).sorted.bed"

    puts "\t\t...creation of the miRNA_*_$g_AnnotSV(genomeBuild).sorted.bed files"
    puts "\t\t   (done only once)"

    file delete -force $miRNARefSeqFileFormatted38or10
    file delete -force $miRNAENSEMBLfileFormatted38
    
    # Validated_miRNA-gene_pairs_hsa_miRBase_v22.1_GRCh38_location_augmented.tsv:
    #    seqnames  start          end             strand  ID              Alias           Name             species   target_gene     confidence
    #    chr11     122146523      122146544       -       MIMAT0010195    MIMAT0010195    hsa-let-7a-2-3p  hsa       CADM1           Functional MTI (Weak)
    #    chr11     122146523      122146544       -       MIMAT0010195    MIMAT0010195    hsa-let-7a-2-3p  hsa       ARID1A          Functional MTI (Weak)
    # Validated_miRNA-gene_pairs_mmu_miRBase_v22.1_GRCm38_location_augmented.tsv
    #    seqnames  start          end             strand  ID              Alias           Name            species    target_gene     confidence
    #    chr13     48538239       48538260        -       MIMAT0000521    MIMAT0000521    mmu-let-7a-5p   mmu        Hoxa9           Functional MTI
    #    chr13     48538239       48538260        -       MIMAT0000521    MIMAT0000521    mmu-let-7a-5p   mmu        Hmga2           Functional MTI
    foreach miRNAfile $miRNAfilesDownloaded {   ;# Validated and predicted
	foreach L [LinesFromFile $miRNAfile] {
	    set Ls [split $L "\t"]
	    if {[regexp "^seqnames" $L]} {
		set i_miRNAchrom [lsearch -exact $Ls "seqnames"]
		set i_miRNAstart [lsearch -exact $Ls "start"]
		set i_miRNAend   [lsearch -exact $Ls "end"]
		set i_miRNAgene  [lsearch -exact $Ls "target_gene"]
		continue
	    }
	    regsub "chr" [lindex $Ls $i_miRNAchrom] "" chrom
	    set miRNAstart [lindex $Ls $i_miRNAstart]
	    set miRNAend [lindex $Ls $i_miRNAend]
	    set coord "$chrom\t$miRNAstart\t$miRNAend"
	    lappend L_Genes($coord) [lindex $Ls $i_miRNAgene]
	    lappend L_coord($chrom) "$coord"
	}
    }
    # Writings:
    file delete -force $miRNARefSeqFileFormatted38or10.tmp
    file delete -force $miRNAENSEMBLfileFormatted38.tmp
    foreach chrom {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M MT} {
	if {![info exists L_coord($chrom)]} {continue}
	set L_linesToWriteRefSeqGRCh38or10 {}
	set L_linesToWriteENSEMBLGRCh38 {}
	foreach coord [lsort -unique $L_coord($chrom)] {
	    set allRefSeqGenes ""
	    set allEnsemblGenes ""
	    foreach gene $L_Genes($coord) {
		if {[isARefSeqGeneName $gene]} {
		    lappend allRefSeqGenes $gene
		}
		if {[isAnENSEMBLgeneName $gene]} {
		    lappend allEnsemblGenes $gene
		}
	    }
	    if {$allRefSeqGenes ne ""} {lappend L_linesToWriteRefSeqGRCh38or10 "$coord\tmTL_miRNA\t[join [lsort -unique $allRefSeqGenes] ";"]"}
	    if {$allEnsemblGenes ne ""} {lappend L_linesToWriteENSEMBLGRCh38 "$coord\tmTL_miRNA\t[join [lsort -unique $allEnsemblGenes] ";"]"}	    
	}
	set L_linesToWriteRefSeqGRCh38or10 [lsort -command AscendingSortOnElement1 [lsort -command AscendingSortOnElement2 $L_linesToWriteRefSeqGRCh38or10]]
	set L_linesToWriteENSEMBLGRCh38 [lsort -command AscendingSortOnElement1 [lsort -command AscendingSortOnElement2 $L_linesToWriteENSEMBLGRCh38]]
	WriteTextInFile "[join $L_linesToWriteRefSeqGRCh38or10 "\n"]" $miRNARefSeqFileFormatted38or10.tmp
	if {$g_AnnotSV(genomeBuild) eq "GRCh38"} {
	    WriteTextInFile "[join $L_linesToWriteENSEMBLGRCh38 "\n"]" $miRNAENSEMBLfileFormatted38.tmp
	}
    }
    
    # Sorting of the bedfile:
    # Intersection with very large files can cause trouble with excessive memory usage.
    # A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm. 
    foreach miRNAfile {miRNARefSeqFileFormatted38or10 miRNAENSEMBLfileFormatted38} {
	if {![file exists [set $miRNAfile].tmp]} {continue} ;# use for Mouse annotation, ENSEMBL annotation doesn't exist
	set sortTmpFile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_sort.tmp.bash"
	ReplaceTextInFile "#!/bin/bash" $sortTmpFile
	WriteTextInFile "# The locale specified by the environment can affects the traditional sort order. We need to use native byte values." $sortTmpFile
	WriteTextInFile "export LC_ALL=C" $sortTmpFile
	WriteTextInFile "sort -k1,1 -k2,2n [set $miRNAfile].tmp >> [set $miRNAfile]" $sortTmpFile
	file attributes $sortTmpFile -permissions 0755
	if {[catch {eval exec bash $sortTmpFile} Message]} {
	    puts "-- checkMiRTargetLinkFiles --"
	    puts "sort -k1,1 -k2,2n [set $miRNAfile].tmp >> [set $miRNAfile]"
	    puts "$Message"
	    puts "Exit with error"
	    exit 2
	}
	file delete -force $sortTmpFile 
	file delete -force [set $miRNAfile].tmp
    }

    # Delete the downloaded files
    foreach miRNAfile $miRNAfilesDownloaded {
	file delete -force $miRNAfile
    }
    
    set g_AnnotSV(miRNAann) 1
}




## Annotate the SV bedFile with the regulatory elements files
## Creation of the g_re global variable
proc regulatoryElementsAnnotation {L_allGenesOverlapped} {

    global g_AnnotSV
    global g_re
    global g_genesReg
    global g_HI
    global g_TS
    global g_reDB

    # Intersect of the input SV bedfile with the regulatory elements files
    ######################################################################
    ## Sorted SV bed file: $g_AnnotSV(bedFile)
    regsub "(.annotated)?.tsv$" $g_AnnotSV(outputDir)/$g_AnnotSV(outputFile) ".SV_RE_intersect.tmp" SV_RE_intersectBEDfile
    file delete -force $SV_RE_intersectBEDfile

    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)"
    set regElementsDir "$extannDir/FtIncludedInSV/RegulatoryElements/$g_AnnotSV(genomeBuild)"
    set L_REfiles {}
    lappend L_REfiles "$regElementsDir/EA_${g_AnnotSV(tx)}_$g_AnnotSV(genomeBuild).sorted.bed"     ;#EnhancerAtlas
    lappend L_REfiles "$regElementsDir/GH_${g_AnnotSV(tx)}_$g_AnnotSV(genomeBuild).sorted.bed"     ;#GeneHancer
    lappend L_REfiles "$regElementsDir/promoter_${g_AnnotSV(promoterSize)}bp_${g_AnnotSV(tx)}_$g_AnnotSV(genomeBuild).sorted.bed" ;# promoter
    lappend L_REfiles "$regElementsDir/miRTargetLink_${g_AnnotSV(tx)}_$g_AnnotSV(genomeBuild).sorted.bed"     ;#miRTargetLink
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

    # Definition of the g_genesReg($SVfromBED) variable (list all the regulated genes impacted by 1 SV)
    ###################################################################################################
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
	    # SVfromBED ("chrom\tstart\tend")
	    set SVfromBED "[join [lrange $Ls 0 2] "\t"]"
	    set L_genesFromTheLine [split [lindex $Ls end] ";"]
	    lappend g_genesReg($SVfromBED) {*}$L_genesFromTheLine ;# regulated genes
	    lappend L_allRegulatedGenes {*}$L_genesFromTheLine
	    set db [lindex $Ls end-1]
	    foreach g $L_genesFromTheLine {
		lappend g_reDB($SVfromBED,$g) $db
	    }
	}
	close $f
	foreach SV [array names g_genesReg] {
	    set g_genesReg($SV) [lsort -unique $g_genesReg($SV)]
	}
	foreach key [array names g_reDB] {
	    set g_reDB($key) [lsort -unique $g_reDB($key)]
	}
	
	## Delete temporary file
	if {$g_AnnotSV(REreport)} {
	    regsub "(.annotated)?.tsv$" $g_AnnotSV(outputDir)/$g_AnnotSV(outputFile) ".SV_RE_intersect.report" permanentSV_RE_intersectBEDfile
	    file delete -force $permanentSV_RE_intersectBEDfile
	    file rename $SV_RE_intersectBEDfile $permanentSV_RE_intersectBEDfile
	} else {
	    file delete -force $SV_RE_intersectBEDfile
	}
    }

    ## Display
    ##########
    puts "...searching for SV overlaps with a gene or a regulatory elements"
    puts "\t...[llength $L_allGenesOverlapped] genes overlapped with an SV"
    set L_allRegulatedGenes [lsort -unique $L_allRegulatedGenes]
    puts "\t...[llength $L_allRegulatedGenes] genes regulated by a regulatory element which is overlapped with an SV\n"
    
    ## Preparation of the phenotype-driven analysis (Exomiser)
    ## (to be able to access to the exomiser score of a gene with [ExomiserAnnotation $gName "score"])
    ##################################################################################################
    set L_allGenes $L_allGenesOverlapped
    lappend L_allGenes {*}$L_allRegulatedGenes
    set L_allGenes [lsort -unique $L_allGenes]
    if {$g_AnnotSV(hpo) ne "" && $L_allGenes ne ""} {
	runExomiser "$L_allGenes" "$g_AnnotSV(hpo)" 
    }
    
    ## HI/TS information for these regulated genes
    ## -> definition of g_HI($gene) and g_TS($gene)
    ###############################################
    set clingenDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Gene-based/ClinGen"
    set ClinGenFileFormatted [lindex [glob -nocomplain "$clingenDir/*_ClinGenAnnotations.tsv*"] end]
    foreach g $L_allRegulatedGenes {
	set t($g) 1
    }
    if {$ClinGenFileFormatted ne ""} {
	if {[regexp ".gz$" $ClinGenFileFormatted]} {
	    set f [open "| gzip -cd $ClinGenFileFormatted"]	 
	} else {		
	    set f [open "$ClinGenFileFormatted"]
	}	 
	while {! [eof $f]} {
	    set L [gets $f]
	    set Ls [split $L "\t"]
	    set gene [lindex $Ls 0]
	    if {[info exists t($gene)]} {
		set HI [lindex $Ls 1]; if {$HI ne "Not yet evaluated"} {set g_HI($gene) "$HI"}
		set TS [lindex $Ls 2]; if {$TS ne "Not yet evaluated"} {set g_TS($gene) "$TS"}
	    }
	}
    }

    ## Resume RE information for each SV
    ## -> definition of g_re($SVfromBED)
    ####################################
    foreach SV [array names g_genesReg] {
	set g_re($SV) ""
	foreach gName "$g_genesReg($SV)" {
	    if {$g_AnnotSV(REselect2)} {
		# AnnotSV restrict the report of the regulated genes to the ones not present in "Gene_name".
		if {[lsearch -exact $L_allGenesOverlapped $gName] ne "-1"} {continue}
	    }
	    set HI ""
	    catch {set HI "$g_HI($gName)"}
	    set TS ""
	    catch {set TS "$g_TS($gName)"}
	    if {$g_AnnotSV(hpo) ne ""} {
		set exomiserScore "[ExomiserAnnotation $gName "score"]"
	    } else {set exomiserScore ""}
	    
	    set lAnn {}
	    if {$g_AnnotSV(REselect1)} {
		# By default, only the genes entering in one of the following categories are reported:
		#  - OMIM morbid genes
		#  - HI genes (ClinGen HI = 3)
		#  - TS genes (ClinGen TS = 3)
		#  - Phenotype matched genes (Exomiser gene score > 0.7)
		#  - User candidate genes 
		if {$HI eq "3"} {lappend lAnn "HI=$HI"}
		if {$TS eq "3"} {lappend lAnn "TS=$TS"}
		if {$exomiserScore ne "" && $exomiserScore > 0.7} {lappend lAnn "EX=$exomiserScore"}
		if {[isMorbid $gName]} {lappend lAnn "morbid"}
		if {[isCandidate $gName]} {lappend lAnn "candidate"}
		if {$lAnn eq ""} {
		    continue
		} else {
		    lappend lAnn "RE=[join $g_reDB($SV,$gName) "+"]"
		}
	    } else {
		if {$HI ne ""} {lappend lAnn "HI=$HI"}
		if {$TS ne ""} {lappend lAnn "TS=$TS"}
		if {$exomiserScore ne "" && $exomiserScore ne "0.0000" && $exomiserScore ne "-1.0"} {lappend lAnn "EX=$exomiserScore"}
		if {[isMorbid $gName]} {lappend lAnn "morbid"}
		if {[isCandidate $gName]} {lappend lAnn "candidate"}
		lappend lAnn "RE=[join $g_reDB($SV,$gName) "+"]"
	    }
	    if {$lAnn ne ""} {
		lappend g_re($SV) "$gName ([join $lAnn "/"])"
	    } else {
		lappend g_re($SV) "$gName"
	    }
	}
	set g_re($SV) [join $g_re($SV) ";"]
    }

    return
}
