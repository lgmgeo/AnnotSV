############################################################################################################
# AnnotSV 3.4.2                                                                                            #
#                                                                                                          #
# AnnotSV: An integrated tool for Structural Variations annotation and ranking                             #
#                                                                                                          #
# Copyright (C) 2017-2024 Veronique Geoffroy (veronique.geoffroy@inserm.fr)                                #
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
    puts "\t\t   For GRCh38 use, please lift over these GRCh37 files to GRCh38"
    
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
    set g_AnnotSV(GHann) 0
    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)"
    set regElementsDir "$extannDir/FtIncludedInSV/RegulatoryElements"
    set GHRefSeqFileFormatted  "$regElementsDir/$g_AnnotSV(genomeBuild)/GH_RefSeq_$g_AnnotSV(genomeBuild).sorted.bed"
    set GHENSEMBLfileFormatted "$regElementsDir/$g_AnnotSV(genomeBuild)/GH_ENSEMBL_$g_AnnotSV(genomeBuild).sorted.bed"
    
    if {[file exists $GHRefSeqFileFormatted] && [file exists $GHENSEMBLfileFormatted]} {
        set g_AnnotSV(GHann) 1
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
    set g_AnnotSV(GHann) 1
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


##  ABC
#################
## - Check if "AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz" has been downloaded
##
## - Check and create if necessary:
##   - ABC_RefSeq_GRCh37.sorted.bed
##   - ABC_ENSEMBL_GRCh37.sorted.bed
##
## - The GRCh38 version should be manually created by lift over with the UCSC web server, sorted, and
##   move in the “$ANNOTSV/share/AnnotSV/Annotations_Human/FtIncludedInSV/RegulatoryElements/GRCh38” directory.
proc checkABCfiles {} {
    
    global g_AnnotSV
    
    
    ## Check if ABC files have been formatted
    ########################################
    set g_AnnotSV(ABCann) 0
    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)"
    set regElementsDir "$extannDir/FtIncludedInSV/RegulatoryElements"
    set necessaryABCfile "$regElementsDir/$g_AnnotSV(genomeBuild)/ABC_${g_AnnotSV(tx)}_$g_AnnotSV(genomeBuild).sorted.bed"
    if {[file exists $necessaryABCfile]} {set g_AnnotSV(ABCann) 1}
    
    set ABCRefSeqFileFormattedGRCh37 "$regElementsDir/GRCh37/ABC_RefSeq_GRCh37.sorted.bed" ;# GRCh38 version should be manually created by lift over
    set ABCENSEMBLfileFormattedGRCh37 "$regElementsDir/GRCh37/ABC_ENSEMBL_GRCh37.sorted.bed" ;# GRCh38 version should be manually created by lift over
    if {[file exists $ABCRefSeqFileFormattedGRCh37] && [file exists $ABCENSEMBLfileFormattedGRCh37]} {
        return
    }
    
    ###############################################
    ## No formatted ABC file available at this step
    ###############################################
    
    ## Check if the ABC file has been downloaded
    ############################################
    set DownloadedABCfile "$regElementsDir/GRCh37/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz" ;# only distributed in GRCh37
    if {![file exists $DownloadedABCfile]} {
        return
    }
    
    ## Downloaded ABC file is available
    ###################################
    # Line example:
    # chr1    710010  710210  genic|chr1:709860-710360        genic   0.420794        LOC100288069    714068  NA      0.815301        True    3958    False   0.125372        0.125372        0.123436        0.123436      0.001216 0.124652        0.052453        0.025955        0.051941        0.021762        brite_adipose-Loft2014
    #
    # Format:
    # 1. chr: enhancer chromosome (eg chr1)
    # 2. start: enhancer start position (eg 710010)
    # 3. end: enhancer end position (eg 710210)
    # 4. name: enhancer name made of its class and its initial (before size normalization) coordinates (eg genic|chr1:709860-710360)
    # 5. class: enhancer class, i.e. position with respect to genes (eg genic)
    # 6. activity_base: geometric mean DHS (or ATAC) and H3K27ac. This is the Activity used in the ABC Score (eg 0.420794)
    # 7. TargetGene: enhancer target gene (eg LOC100288069)
    # 8. TargetGeneTSS: enhancer target gene TSS position (eg 714068)
    # 9. TargetGeneExpression: enhancer target gene expression (eg NA if no gene expression but the activity of the gene promoter was used)
    # 10. TargetGenePromoterActivityQuantile: quantile of Activity at Promoter of target gene (eg 0.815301)
    # 11. TargetGeneIsExpressed: whether enhancer target gene is expressed (eg True)
    # 12. distance: distance in bp between enhancer and TSS of target gene (eg 3958)
    # 13. isSelfPromoter: boolean denoting whether element is the promoter of the TargetGene (eg False)
    # 14. hic_contact: K-R normalized Hi-C contacts between element and TargetGene TSS (eg 0.125372)
    # 15. powerlaw_contact: contact between enhancer and gene predicted by power law in the cell type (eg 0.125372)
    # 16. powerlaw_contact_reference: contact between enhancer and gene predicted by power law in K562 (eg 0.123436)
    # 17. hic_contact_pl_scaled: hic_contact scaled by the difference in powerlaw fits between target cell type and reference cell type (eg 0.123436)
    # 18. hic_pseudocount: pseudocount added to HiC Contact (eg 0.001216)
    # 19. hic_contact_pl_scaled_adj: powerlaw scaled KR Normalized HiC plus pseudocount. This is the Contact used in the ABC Score (eg 0.124652)
    # 20. ABC.Score.Numerator: the numator of the ABC Score. Activity x Contact without dividing by the scores for other elements near the Target Gene (eg 0.052453)
    # 21. ABC.Score: ABC Score for the enhancer/gene relationship in this cell type (eg 0.025955)
    # 22. powerlaw.Score.Numerator: power law score numerator (eg 0.051941)
    # 23. powerlaw.Score: power law score (eg 0.021762)
    # 24. CellType: cell type where the enhancer/gene relationship is seen (eg brite_adipose-Loft2014)
    
    puts "\t...ABC configuration ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    
    ## Creation of ABC_*_GRCh37.sorted.bed
    ##########################################
    ##   - # Header: chr ABCstart ABCend ABCtype ABCgenes
    
    puts "\t\t...creation of [file tail $ABCRefSeqFileFormattedGRCh37] and [file tail $ABCENSEMBLfileFormattedGRCh37]"
    puts "\t\t   (done only once)"
    puts "\t\t   For GRCh38 use, please lift over these GRCh37 files to GRCh38"
    
    file delete -force $ABCRefSeqFileFormattedGRCh37
    file delete -force $ABCENSEMBLfileFormattedGRCh37
    
    # Writings:
    file delete -force $ABCRefSeqFileFormattedGRCh37.tmp
    file delete -force $ABCENSEMBLfileFormattedGRCh37.tmp
    
    set f [open "| gzip -cd $DownloadedABCfile"]
    while {! [eof $f]} {
        set L [gets $f]
        set Ls [split $L "\t"]
        if {[regexp "^chr\t" $L]} {
            set i_start [lsearch -exact $Ls "start"]
            set i_end [lsearch -exact $Ls "end"]
            set i_gene [lsearch -exact $Ls "TargetGene"]
            continue
        }
        regsub "chr" [lindex $Ls 0] "" chrom
        set start [lindex $Ls $i_start]
        set end [lindex $Ls $i_end]
        set gene [lindex $Ls $i_gene]
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
    foreach chrom {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M MT} {
        # RefSeq
        if {[info exists L_coordRefSeq($chrom)]} {
            set L_coordRefSeq($chrom) [lsort -command AscendingSortOnElement1 [lsort -command AscendingSortOnElement2 [lsort -unique $L_coordRefSeq($chrom)]]]
            set L_linesToWriteGRCh37 {}
            foreach coord $L_coordRefSeq($chrom) {
                set L_genesRefSeq($coord) [lsort -unique $L_genesRefSeq($coord)]
                lappend L_linesToWriteGRCh37 "$coord\tABC_enhancer\t[join $L_genesRefSeq($coord) ";"]"
            }
            WriteTextInFile "[join $L_linesToWriteGRCh37 "\n"]" $ABCRefSeqFileFormattedGRCh37.tmp
        }
        # ENSEMBL
        if {[info exists L_coordENSEMBL($chrom)]} {
            set L_coordENSEMBL($chrom) [lsort -command AscendingSortOnElement1 [lsort -command AscendingSortOnElement2 [lsort -unique $L_coordENSEMBL($chrom)]]]
            set L_linesToWriteGRCh37 {}
            foreach coord $L_coordENSEMBL($chrom) {
                set L_genesENSEMBL($coord) [lsort -unique $L_genesENSEMBL($coord)]
                lappend L_linesToWriteGRCh37 "$coord\tABC_enhancer\t[join $L_genesENSEMBL($coord) ";"]"
            }
            WriteTextInFile "[join $L_linesToWriteGRCh37 "\n"]" $ABCENSEMBLfileFormattedGRCh37.tmp
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
    WriteTextInFile "sort -k1,1 -k2,2n $ABCRefSeqFileFormattedGRCh37.tmp >> $ABCRefSeqFileFormattedGRCh37" $sortTmpFile
    file attributes $sortTmpFile -permissions 0755
    if {[catch {eval exec bash $sortTmpFile} Message]} {
        puts "-- checkABCfiles --"
        puts "sort -k1,1 -k2,2n $ABCRefSeqFileFormattedGRCh37.tmp >> $ABCRefSeqFileFormattedGRCh37"
        puts "$Message"
        puts "Exit with error"
        exit 2
    }
    file delete -force $sortTmpFile
    file delete -force $ABCRefSeqFileFormattedGRCh37.tmp
    # ENSEMBL
    set sortTmpFile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_sort.tmp.bash"
    ReplaceTextInFile "#!/bin/bash" $sortTmpFile
    WriteTextInFile "# The locale specified by the environment can affects the traditional sort order. We need to use native byte values." $sortTmpFile
    WriteTextInFile "export LC_ALL=C" $sortTmpFile
    WriteTextInFile "sort -k1,1 -k2,2n $ABCENSEMBLfileFormattedGRCh37.tmp >> $ABCENSEMBLfileFormattedGRCh37" $sortTmpFile
    file attributes $sortTmpFile -permissions 0755
    if {[catch {eval exec bash $sortTmpFile} Message]} {
        puts "-- checkABCfiles --"
        puts "sort -k1,1 -k2,2n $ABCENSEMBLfileFormattedGRCh37.tmp >> $ABCENSEMBLfileFormattedGRCh37"
        puts "$Message"
        puts "Exit with error"
        exit 2
    }
    file delete -force $sortTmpFile
    file delete -force $ABCENSEMBLfileFormattedGRCh37.tmp
    
    # Delete the downloaded files
    if {[file exists $ABCRefSeqFileFormattedGRCh37] && [file exists $ABCENSEMBLfileFormattedGRCh37]} {
        file delete -force $DownloadedABCfile
    }
    
    set necessaryABCfile "$regElementsDir/$g_AnnotSV(genomeBuild)/ABC_${g_AnnotSV(tx)}_$g_AnnotSV(genomeBuild).sorted.bed"
    if {[file exists $necessaryABCfile]} {set g_AnnotSV(ABCann) 1}
}



##  MPRA
#################
## - Check if "MPRA_promoters.tsv" and "MPRA_enhancers.tsv" have been copy/paste
##   (https://kircherlab.bihealth.org/satMutMPRA)
##
## - Check and create if necessary:
##   - MPRA_RefSeq_GRCh37.sorted.bed
##   - MPRA_ENSEMBL_GRCh37.sorted.bed
##   - MPRA_RefSeq_GRCh38.sorted.bed
##   - MPRA_ENSEMBL_GRCh38.sorted.bed
proc checkMPRAfiles {} {
    
    global g_AnnotSV
    
    
    ## Check if MPRA files have been formatted
    ########################################
    set g_AnnotSV(MPRAann) 0
    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)"
    set regElementsDir "$extannDir/FtIncludedInSV/RegulatoryElements"
    set necessaryMPRAfile "$regElementsDir/$g_AnnotSV(genomeBuild)/MPRA_${g_AnnotSV(tx)}_$g_AnnotSV(genomeBuild).sorted.bed"
    if {[file exists $necessaryMPRAfile]} {set g_AnnotSV(MPRAann) 1; return}
    
    set MPRARefSeqFileFormattedGRCh37 "$regElementsDir/GRCh37/MPRA_RefSeq_GRCh37.sorted.bed"
    set MPRAENSEMBLfileFormattedGRCh37 "$regElementsDir/GRCh37/MPRA_ENSEMBL_GRCh37.sorted.bed"
    set MPRARefSeqFileFormattedGRCh38 "$regElementsDir/GRCh38/MPRA_RefSeq_GRCh38.sorted.bed"
    set MPRAENSEMBLfileFormattedGRCh38 "$regElementsDir/GRCh38/MPRA_ENSEMBL_GRCh38.sorted.bed"
    if {[file exists $MPRARefSeqFileFormattedGRCh37] && [file exists $MPRAENSEMBLfileFormattedGRCh37] && [file exists $MPRARefSeqFileFormattedGRCh38] && [file exists $MPRAENSEMBLfileFormattedGRCh38]} {
        return
    }
    
    ###############################################
    ## No formatted MPRA file available at this step
    ###############################################
    
    ## Check if the MPRA files have been downloaded
    ###############################################
    set DownloadedMPRAfile_promoters "$regElementsDir/MPRA_promoters.tsv" ;# Distributed with GRCh37/GRCh38 coordinates in the same file
    set DownloadedMPRAfile_enhancers "$regElementsDir/MPRA_enhancers.tsv" ;# Distributed with GRCh37/GRCh38 coordinates in the same file
    if {![file exists $DownloadedMPRAfile_promoters] || ![file exists $DownloadedMPRAfile_enhancers]} {
        return
    }
    
    ## Downloaded MPRA file is available
    ###################################
    ## MPRA_enhancers.tsv format:
    #  Name    Genomic coordinates (GRCh37)    Genomic coordinates (GRCh38)    Transcript      Associated Phenotype    Luciferase vector       MPRA vector     Cell line	...
    ## MPRA_enhancers.tsv format:
    #  Name    Genomic coordinates (GRCh37)    Genomic coordinates (GRCh38)    Associated Phenotype    Luciferase vector       MPRA vector     Cell line       Transf. time (hr)		...
    
    puts "\t...MPRA configuration ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    
    ## Creation of MPRA_*_*.sorted.bed
    ##########################################
    # Header: chr MPRAstart MPRAend MPRAtype MPRAgenes
    
    puts "\t\t...creation of [file tail $MPRARefSeqFileFormattedGRCh37], [file tail $MPRAENSEMBLfileFormattedGRCh37], [file tail $MPRARefSeqFileFormattedGRCh38] and [file tail $MPRAENSEMBLfileFormattedGRCh38]"
    puts "\t\t   (done only once)"
    
    file delete -force $MPRARefSeqFileFormattedGRCh37
    file delete -force $MPRAENSEMBLfileFormattedGRCh37
    file delete -force $MPRARefSeqFileFormattedGRCh38
    file delete -force $MPRAENSEMBLfileFormattedGRCh38
    
    # Writings:
    file delete -force $MPRARefSeqFileFormattedGRCh37.tmp
    file delete -force $MPRAENSEMBLfileFormattedGRCh37.tmp
    file delete -force $MPRARefSeqFileFormattedGRCh38.tmp
    file delete -force $MPRAENSEMBLfileFormattedGRCh38.tmp
    
    set L_toParse [LinesFromFile $DownloadedMPRAfile_promoters]
    lappend L_toParse {*}[LinesFromFile $DownloadedMPRAfile_enhancers]
    
    foreach L $L_toParse {
        set Ls [split $L "\t"]
        if {[regexp "^Name\t" $L]} {
            set i_coord37 [lsearch -exact $Ls "Genomic coordinates (GRCh37)"] ; #e.g. "chrX:138,612,622-138,612,924"
            set i_coord38 [lsearch -exact $Ls "Genomic coordinates (GRCh38)"]
            if {[regexp "Transcript" $L]} {
                set MPRAtype "MPRA_promoter"
            } else {
                set MPRAtype "MPRA_enhancer"
            }
            continue
        }
        regsub -all "," [lindex $Ls $i_coord37] "" coord37
        regsub -all "," [lindex $Ls $i_coord38] "" coord38
        if {![regexp "chr(\[0-9XYMT\]+):(\[0-9\]+)-(\[0-9\]+)" $coord37 match chrom37 start37 end37]} {
            puts "Bad format in $DownloadedMPRAfile_promoters"
            puts "$L"
            puts "Check this format."
        }
        if {![regexp "chr(\[0-9XYMT\]+):(\[0-9\]+)-(\[0-9\]+)" $coord38 match chrom38 start38 end38]} {
            puts "Bad format in $DownloadedMPRAfile_promoters"
            puts "$L"
            puts "Check this format."
        }
        set gene [lindex $Ls 0]
        regsub " \\(.+|\\+.+" $gene "" gene
        
        set type($chrom37\t$start37\t$end37) "$MPRAtype"
        set type($chrom38\t$start38\t$end38) "$MPRAtype"
        # RefSeq
        if {[isARefSeqGeneName $gene]} {
            lappend L_coordRefSeq37($chrom37) "$chrom37\t$start37\t$end37"
            lappend L_genesRefSeq37($chrom37\t$start37\t$end37) "$gene"
            lappend L_coordRefSeq38($chrom38) "$chrom38\t$start38\t$end38"
            lappend L_genesRefSeq38($chrom38\t$start38\t$end38) "$gene"
        }
        # ENSEMBL
        if {[isAnENSEMBLgeneName $gene]} {
            lappend L_coordENSEMBL37($chrom37) "$chrom37\t$start37\t$end37"
            lappend L_genesENSEMBL37($chrom37\t$start37\t$end37) "$gene"
            lappend L_coordENSEMBL38($chrom38) "$chrom38\t$start38\t$end38"
            lappend L_genesENSEMBL38($chrom38\t$start38\t$end38) "$gene"
        }
    }
    foreach chrom {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M MT} {
        foreach build {"37" "38"} {
            # RefSeq
            if {[info exists L_coordRefSeq[set build]($chrom)]} {
                set L_coordRefSeq[set build]($chrom) [lsort -command AscendingSortOnElement1 [lsort -command AscendingSortOnElement2 [lsort -unique [set L_coordRefSeq[set build]($chrom)]]]]
                set L_linesToWrite[set build] {}
                foreach coord [set L_coordRefSeq[set build]($chrom)] {
                    set L_genesRefSeq[set build]($coord) [lsort -unique [set L_genesRefSeq[set build]($coord)]]
                    lappend L_linesToWrite[set build] "$coord\t$type($coord)\t[join [set L_genesRefSeq[set build]($coord)] ";"]"
                }
                WriteTextInFile "[join [set L_linesToWrite[set build]] "\n"]" [set MPRARefSeqFileFormattedGRCh[set build]].tmp
            }
            # ENSEMBL
            if {[info exists L_coordENSEMBL[set build]($chrom)]} {
                set L_coordENSEMBL[set build]($chrom) [lsort -command AscendingSortOnElement1 [lsort -command AscendingSortOnElement2 [lsort -unique [set L_coordENSEMBL[set build]($chrom)]]]]
                set L_linesToWrite[set build] {}
                foreach coord [set L_coordENSEMBL[set build]($chrom)] {
                    set L_genesENSEMBL[set build]($coord) [lsort -unique [set L_genesENSEMBL[set build]($coord)]]
                    lappend L_linesToWrite[set build] "$coord\t$type($coord)\t[join [set L_genesENSEMBL[set build]($coord)] ";"]"
                }
                WriteTextInFile "[join [set L_linesToWrite[set build]] "\n"]" [set MPRAENSEMBLfileFormattedGRCh[set build]].tmp
            }
        }
    }
    
    
    foreach build {"37" "38"} {
        # Sorting of the bedfile:
        # Intersection with very large files can cause trouble with excessive memory usage.
        # A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm.
        # RefSeq
        set sortTmpFile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_sort.tmp.bash"
        ReplaceTextInFile "#!/bin/bash" $sortTmpFile
        WriteTextInFile "# The locale specified by the environment can affects the traditional sort order. We need to use native byte values." $sortTmpFile
        WriteTextInFile "export LC_ALL=C" $sortTmpFile
        WriteTextInFile "sort -k1,1 -k2,2n [set MPRARefSeqFileFormattedGRCh[set build]].tmp >> [set MPRARefSeqFileFormattedGRCh[set build]]" $sortTmpFile
        file attributes $sortTmpFile -permissions 0755
        if {[catch {eval exec bash $sortTmpFile} Message]} {
            puts "-- checkMPRAfiles --"
            puts "sort -k1,1 -k2,2n [set MPRARefSeqFileFormattedGRCh[set build]].tmp >> [set MPRARefSeqFileFormattedGRCh[set build]]"
            puts "$Message"
            puts "Exit with error"
            exit 2
        }
        file delete -force $sortTmpFile
        file delete -force [set MPRARefSeqFileFormattedGRCh[set build]].tmp
        # ENSEMBL
        set sortTmpFile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_sort.tmp.bash"
        ReplaceTextInFile "#!/bin/bash" $sortTmpFile
        WriteTextInFile "# The locale specified by the environment can affects the traditional sort order. We need to use native byte values." $sortTmpFile
        WriteTextInFile "export LC_ALL=C" $sortTmpFile
        WriteTextInFile "sort -k1,1 -k2,2n [set MPRAENSEMBLfileFormattedGRCh[set build]].tmp >> [set MPRAENSEMBLfileFormattedGRCh[set build]]" $sortTmpFile
        file attributes $sortTmpFile -permissions 0755
        if {[catch {eval exec bash $sortTmpFile} Message]} {
            puts "-- checkMPRAfiles --"
            puts "sort -k1,1 -k2,2n [set MPRAENSEMBLfileFormattedGRCh[set build]].tmp >> [set MPRAENSEMBLfileFormattedGRCh[set build]]"
            puts "$Message"
            puts "Exit with error"
            exit 2
        }
        file delete -force $sortTmpFile
        file delete -force [set MPRAENSEMBLfileFormattedGRCh[set build]].tmp
    }
    
    
    # Delete the downloaded files
    if {[file exists $MPRARefSeqFileFormattedGRCh37] && [file exists $MPRAENSEMBLfileFormattedGRCh37] && [file exists $MPRARefSeqFileFormattedGRCh38] && [file exists $MPRAENSEMBLfileFormattedGRCh38]} {
        file delete -force $DownloadedMPRAfile_promoters
        file delete -force $DownloadedMPRAfile_enhancers
    }
    
    set necessaryMPRAfile "$regElementsDir/$g_AnnotSV(genomeBuild)/MPRA_${g_AnnotSV(tx)}_$g_AnnotSV(genomeBuild).sorted.bed"
    if {[file exists $necessaryMPRAfile]} {set g_AnnotSV(MPRAann) 1}
}





##################################################################################################################
##################################################################################################################
##################################################################################################################




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
    lappend L_REfiles "$regElementsDir/ABC_${g_AnnotSV(tx)}_$g_AnnotSV(genomeBuild).sorted.bed"; #Activity-by-Contact (ABC) Model data
    lappend L_REfiles "$regElementsDir/MPRA_${g_AnnotSV(tx)}_$g_AnnotSV(genomeBuild).sorted.bed"; #Massive parallel reporter assays (MPRA) data (Kircher lab)
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
    
    ## Preparation of the phenotype-driven analysis (PhenoGenius + Exomiser)
    ## (to be able to access to the exomiser score of a gene with [ExomiserAnnotation $gName "score"])
    ## (to be able to access to the PhenoGenius specificity of a gene with [PhenoGeniusAnnotation $gName "specificity"])
    ####################################################################################################################
    set L_allGenes $L_allGenesOverlapped
    lappend L_allGenes {*}$L_allRegulatedGenes
    set L_allGenes [lsort -unique $L_allGenes]
	if {$g_AnnotSV(PhenoGenius)} {		
		set L_NCBI_ID {}
		foreach g $L_allGenes {
			set NCBI_ID [searchforGeneID $g]
			if {$NCBI_ID ne ""} {
				lappend L_NCBI_ID $NCBI_ID
			}
		}
		set L_NCBI_ID [lsort -unique $L_NCBI_ID]
		runPhenoGenius "$L_allGenes" "$L_NCBI_ID" "$g_AnnotSV(hpo)"
    }
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
		set re_priority($SV) ""
		set re_other($SV) ""
		set re_noann($SV) ""
        foreach gName "$g_genesReg($SV)" {
            if {$g_AnnotSV(REselect2)} {
                # AnnotSV restrict the report of the regulated genes to the ones not present in "Gene_name".
                if {[lsearch -exact $L_allGenesOverlapped $gName] ne "-1"} {continue}
            }
            set HI ""
            catch {set HI "$g_HI($gName)"}
            set TS ""
            catch {set TS "$g_TS($gName)"}
            if {$g_AnnotSV(PhenoGenius)} {
                set PhenoGeniusSpecificity "[PhenoGeniusAnnotation $gName "specificity"]"
            } else {set PhenoGeniusSpecificity ""}
            if {$g_AnnotSV(hpo) ne ""} {
                set exomiserScore "[ExomiserAnnotation $gName "score"]"
            } else {set exomiserScore ""}
            
            set lAnn {}
			set priority 0
            if {$g_AnnotSV(REselect1)} {
                # By default, only the genes entering in one of the following categories are reported:
                #  - OMIM morbid genes
                #  - HI genes (ClinGen HI = 3)
                #  - TS genes (ClinGen TS = 3)
                #  - Phenotype matched genes (PhenoGenius specificity = "A" or Exomiser gene score > 0.7)
                #  - User candidate genes
                if {$HI eq "3"} {lappend lAnn "HI=$HI"}
                if {$TS eq "3"} {lappend lAnn "TS=$TS"}
				if {$PhenoGeniusSpecificity eq "A"} {lappend lAnn "PG=A"; set priority 1}
                if {$exomiserScore ne "" && $exomiserScore > 0.7} {lappend lAnn "EX=$exomiserScore"; set priority 1}
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
                if {$PhenoGeniusSpecificity ne ""} {lappend lAnn "PG=A"}
                if {$exomiserScore ne "" && $exomiserScore ne "0.0000" && $exomiserScore ne "-1.0"} {lappend lAnn "EX=$exomiserScore"}
                if {[isMorbid $gName]} {lappend lAnn "morbid"}
                if {[isCandidate $gName]} {lappend lAnn "candidate"}
                lappend lAnn "RE=[join $g_reDB($SV,$gName) "+"]"
            }
            if {$lAnn ne ""} {
				if {$priority} {
					lappend re_priority($SV) "$gName ([join $lAnn "/"])"
				} else {
					lappend re_other($SV) "$gName ([join $lAnn "/"])"
				}
            } else {
                lappend re_noann($SV) "$gName"
            }
        }
	 
		# We prioritarily keep RE with "PhenoGenius specificity = "A" or Exomiser gene score > 0.7" (for the ranking)
        if {$re_priority($SV) ne ""} {
            set g_re($SV) $re_priority($SV)
        } elseif {$re_other($SV) ne ""} {
            lappend g_re($SV) {*}$re_other($SV)
        } elseif {$re_noann($SV) ne ""} {
            lappend g_re($SV) {*}$re_noann($SV)
		}

        # AnnotSV restricts the number of overlapping reported features to 50
        if {[llength $g_re($SV)] > 50} {
            set g_re($SV) "[join [lrange $g_re($SV) 0 49] ";"]..."
        } else {
            set g_re($SV) [join $g_re($SV) ";"]
		}
    }
    
    return
}

