############################################################################################################
# AnnotSV 3.1.1                                                                                            #
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                                             #
# GNU General Public License for more details.                                                             #
#                                                                                                          #
# You should have received a copy of the GNU General Public License                                        #
# along with this program; If not, see <http://www.gnu.org/licenses/>.                                     #
############################################################################################################


## - Check and create if necessary the "genes.RefSeq.sorted.bed" file
proc checkGenesRefSeqFile {} {
    
    global g_AnnotSV
    
    ## Check if the Genes file has been downloaded then formatted
    #############################################################
    set genesDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Genes/$g_AnnotSV(genomeBuild)"
    
    set genesFileDownloaded "[glob -nocomplain $genesDir/refGene.txt.gz]"
    set genesFileFormatted "[glob -nocomplain $genesDir/genes.RefSeq.sorted.bed]"
    
    if {$genesFileDownloaded eq "" && $genesFileFormatted eq ""} {
	puts "############################################################################"
	puts "\"$genesDir/refGene.txt.gz\" file doesn't exist"
	puts "Please check your install - Exit with error."
	puts "############################################################################"
	exit 2
    }
    
    if {$genesFileFormatted eq ""} {
	
	## Delete promoters files (need to be updated after the creation of new genes file)
	##################################################################################### 
	set promoterDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/FtIncludedInSV/Promoter/$g_AnnotSV(genomeBuild)"
	foreach promFile [glob -nocomplain "$promoterDir/promoter_*bp_RefSeq_$g_AnnotSV(genomeBuild).sorted.bed"] {
	    file delete -force $promFile
	}
	
	## - Create the "genes.RefSeq.sorted.bed"
	#####################################
	set genesFileFormatted "$genesDir/genes.RefSeq.sorted.bed"
	puts "\t...creation of $genesFileFormatted ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	puts "\t   (done only once during the first annotation)"
	
	# Removing non-standard contigs (other than the standard 1-22,X,Y,MT) and sorting the file in karyotypic order
	# -> create L_genesTXTsorted
	
	## Save the line of the genesTXT by chromosome (L_lines($chrom))
	foreach L [LinesFromGZFile $genesFileDownloaded] {
	    set Ls [split $L "\t"]
	    set chrom [lindex $Ls 2]
	    lappend L_lines($chrom) "$L"
	}
	## Sorting in karyotypic order
	set L_genesTXTsorted {}
	foreach val {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M MT} {
	    if {![info exists L_lines(chr$val)]} {continue}
	    lappend L_genesTXTsorted {*}"[lsort -command AscendingSortOnElement4 [lsort -command AscendingSortOnElement5 $L_lines(chr$val)]]"
	}
	# Creation of the $genesFileFormatted.
	# Chromosome nomenclature used: 1-22,X,Y,M,T (without "chr")
	## INPUT:    #bin name chrom strand txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds score name2 cdsStartStat cdsEndStat exonFrames
	## OUTPUT:   chrom txStart txEnd strand name2 name cdsStart cdsEnd exonStarts exonEnds
	# WARNING : NR_* are for non coding RNA. However, cdsStart=cdsEnd, => CDSlength=1
	foreach L $L_genesTXTsorted {
	    set Ls [split $L "\t"]
	    regsub "chr" [lindex $Ls 2] "" chrom
	    set line "$chrom\t[lindex $Ls 4]\t[lindex $Ls 5]\t[lindex $Ls 3]\t[lindex $Ls 12]\t[lindex $Ls 1]\t[lindex $Ls 6]\t[lindex $Ls 7]\t[lindex $Ls 9]\t[lindex $Ls 10]"
	    if {![info exists infos($line)]} {
		WriteTextInFile $line $genesFileFormatted.tmp
		set infos($line) 1
	    }
	}
	file delete -force $genesFileDownloaded
	# Sorting of the bedfile:
	# Intersection with very large files can cause trouble with excessive memory usage.
	# A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm.
	set sortTmpFile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_sort.tmp.bash"
	ReplaceTextInFile "#!/bin/bash" $sortTmpFile
	WriteTextInFile "# The locale specified by the environment can affects the traditional sort order. We need to use native byte values." $sortTmpFile
	WriteTextInFile "export LC_ALL=C" $sortTmpFile
	WriteTextInFile "sort -k1,1 -k2,2n $genesFileFormatted.tmp > $genesFileFormatted" $sortTmpFile
	file attributes $sortTmpFile -permissions 0755
	if {[catch {eval exec bash $sortTmpFile} Message]} {
	    puts "-- checkGenesRefSeqFile --"
	    puts "sort -k1,1 -k2,2n $genesFileFormatted.tmp > $genesFileFormatted"
	    puts "$Message"
	    puts "Exit with error"
	    exit 2
	}
	file delete -force $sortTmpFile 
	file delete -force $genesFileFormatted.tmp
	
    }
    
    # DISPLAY:
    ##########
    set g_AnnotSV(genesFile) $genesFileFormatted
    
}



## - Check the "genes.ENSEMBL.sorted.bed" file
proc checkGenesENSEMBLfile {} {
    
    global g_AnnotSV
    
    ## Check if the formatted ENSEMBL genes file is present
    #######################################################
    set genesDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Genes/$g_AnnotSV(genomeBuild)"
    
    set GenesENSEMBLfileFormatted [glob -nocomplain $genesDir/genes.ENSEMBL.sorted.bed]
    
    if {$GenesENSEMBLfileFormatted eq ""} {
        puts "############################################################################"
        puts "\"$genesDir/genes.ENSEMBL.sorted.bed\" doesn't exist"
        puts "Please check your install - Exit with error."
        puts "############################################################################"
        exit 2
    }
    
    # DISPLAY:
    ##########
    set g_AnnotSV(genesFile) $GenesENSEMBLfileFormatted
}



## Annotate the SV bedFile with the genes file.
## Keep only 1 transcript annotation by gene:
##   - the one selected by the user with the "-txFile" option
##   - the one with the most of "bp from CDS" (=overlappedCDSlength)
##   - if x transcript with same "bp from CDS", the one with the most of "bp from UTR, exon, intron" (=overlappedTxLength)
##
## Creation of FullAndSplitBedFile ($g_AnnotSV(outputDir)/$g_AnnotSV(outputFile).tmp)
## -> formatted and sorted
proc genesAnnotation {} {

    global g_AnnotSV
    global g_Lgenes


    # Check the -svtBEDcol and -samplesidBEDcol options
    # Create the -svtTSVcol variable
    #####################################################################################
    checksvtBEDcol $g_AnnotSV(bedFile)
    checksamplesidBEDcol $g_AnnotSV(bedFile)

    # Bedfile should be sorted and should not have "chr" in the first column
    ########################################################################
    # Removing non-standard contigs (other than the standard 1-22,X,Y,MT) and sorting the file in karyotypic order
    set f [open $g_AnnotSV(bedFile)]
    set test 0
    while {![eof $f]} {
	set L [gets $f]
	if {$L eq ""} {continue}
	
	# To avoid a bug with "bedtools intersect", the last value of a line could NOT be empty.
	# Exemple of bedtools error message:
	# "Error: line number 716 of file input.formatted.sorted.bed has 6 fields, but 4 were expected."
	regsub "\t$" $L "\t." L
	
	# To avoid a bug, the SV type value should be trimmed (should not be equal to "DEL ")
	# For each tab separated value, any leading or trailing white space is removed
	set Ls [split $L "\t"]
	set newLs ""
	foreach val $Ls {
	    lappend newLs [string trim $val]
	}
	set L [join $newLs "\t"]

	if {[regsub "chr" [lindex $L 0] "" chrom] ne 0} {
	    lappend L_Text($chrom) "[string range $L 3 end]"
	} else {
	    lappend L_Text($chrom) "$L"
	}
    }
    close $f
    ## Writing the bedfile (not sorted)
    set emptyBed 1
    regsub -nocase ".bed$" $g_AnnotSV(bedFile) ".formatted.sorted.bed" newBed
    set newBed "$g_AnnotSV(outputDir)/[file tail $newBed]"
    file delete -force $newBed
    foreach chrom {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M MT} {
	if {![info exists L_Text($chrom)]} {continue}
	set emptyBed 0
	WriteTextInFile [join $L_Text($chrom) "\n"] $newBed.tmp
	unset L_Text($chrom)
    }
    if {$emptyBed} {
	puts "The SV input file is empty or badly formatted."
	puts "Exit with error\n"
	exit 2
    }
    
    # Sorting of the bedfile:
    # Intersection with very large files can cause trouble with excessive memory usage.
    # A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm.
    set sortTmpFile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_sort.tmp.bash"
    ReplaceTextInFile "#!/bin/bash" $sortTmpFile
    WriteTextInFile "# The locale specified by the environment can affects the traditional sort order. We need to use native byte values." $sortTmpFile
    WriteTextInFile "export LC_ALL=C" $sortTmpFile
    WriteTextInFile "sort -k1,1 -k2,2n $newBed.tmp >> $newBed" $sortTmpFile
    file attributes $sortTmpFile -permissions 0755
    if {[catch {eval exec bash $sortTmpFile} Message]} {
	puts "-- genesAnnotation --"
	puts "sort -k1,1 -k2,2n $newBed.tmp >> $newBed"
	puts "$Message"
	puts "Exit with error\n"
	exit 2
    }
    file delete -force $sortTmpFile 
    file delete -force $newBed.tmp
    set g_AnnotSV(bedFile) $newBed

    # Used for the insertion of the "full/split" information
    set L_Bed [LinesFromFile $newBed]

    # OUTPUT:
    ###############
    set FullAndSplitBedFile "$g_AnnotSV(outputDir)/$g_AnnotSV(outputFile).tmp" ;# DON'T REMOVE IT! Used in AnnotSV-write.tcl
    set splitBedFile "$FullAndSplitBedFile.tmp"
    file delete -force $FullAndSplitBedFile
    file delete -force $splitBedFile
    

    # Intersect of the input SV bedfile with the genes annotations
    # + create the "FullAndSplitBedFile"
    ###################################################################################################

    set L_linesToWrite {}
    set i_toWrite 0
    ## -> Creation of the split annotations
    if {[catch {exec $g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(bedFile) -b $g_AnnotSV(genesFile) -wa -wb > $splitBedFile} Message]} {
	if {[catch {exec $g_AnnotSV(bedtools) intersect -a $g_AnnotSV(bedFile) -b $g_AnnotSV(genesFile) -wa -wb > $splitBedFile} Message]} {
	    puts "-- genesAnnotation --"
	    puts "$g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(bedFile) -b $g_AnnotSV(genesFile) -wa -wb > $splitBedFile"
	    puts "$Message"
	    puts "Exit with error"
	    exit 2
	}
    }
    if {[file size $splitBedFile] eq 0} {
	puts "...no intersection between SV and gene annotation"
	set n 0
	while {$n < [llength $L_Bed]} {
	    lappend L_linesToWrite "[lindex $L_Bed $n]\tfull"
	    incr i_toWrite
	    incr n
	}
	WriteTextInFile [join $L_linesToWrite "\n"] "$FullAndSplitBedFile"

	## Delete temporary file
	file delete -force $splitBedFile
	return
    }

    # List of the user selected transcripts
    #######################################
    set L_selectedTx {}
    if {$g_AnnotSV(txFile) ne ""} {
	foreach L [LinesFromFile $g_AnnotSV(txFile)] {
	    foreach tx [split $L " |\t"] {
		# Remove the version number if present
		regsub "\\.\[0-9\]+" $tx "" tx
		lappend L_selectedTx $tx
	    }
	}
    }

    # Parse the "splitBedFile" 
    ##########################
    set L_allGenesOverlapped {}
    set L_genes {}
    set L "[FirstLineFromFile $splitBedFile]"
    set Ls [split $L "\t"]
    set splitSVleft [lindex $Ls 1]
    set splitSVright [lindex $Ls 2]

    # oldsplitSV (a complete line from the bedfile)
    # Compulsory: In case of same SV presents several times in the bed file (but with different user annotation), a complete line should be kept
    if {$g_AnnotSV(svtBEDcol) ne "-1"} {
	set SVtype "\t[lindex $Ls "$g_AnnotSV(svtBEDcol)"]"
    } else {
	set SVtype ""
    }
    set shortOldSplitSV "[join [lrange $Ls 0 2] "\t"]"; # chrom, start, end
    set oldSplitSV "[join [lrange $Ls 0 end-10] "\t"]"

    # SVfromBED (a complete line from the bedfile)
    set n 0
    set Ls_fromBED "[split [lindex $L_Bed $n] "\t"]"
    if {$g_AnnotSV(svtBEDcol) ne "-1"} {
	set SVtypeFromBED "\t[lindex $Ls_fromBED "$g_AnnotSV(svtBEDcol)"]"
    } else {
	set SVtypeFromBED ""
    }
    set SVfromBED "[join $Ls_fromBED "\t"]"
	
    # Particular case: only 1 SV is associated with split lines AND it is not the first SV of the BED.
    while {$SVfromBED ne $oldSplitSV} {
	# Writing of the "full" SV line (not present in the $splitBedFile file, if not covering a gene)
	lappend L_linesToWrite "[lindex $L_Bed $n]\tfull"
	incr i_toWrite
	incr n
	set Ls_fromBED "[split [lindex $L_Bed $n] "\t"]"
	if {$g_AnnotSV(svtBEDcol) ne "-1"} {
	    set SVtypeFromBED "\t[lindex $Ls_fromBED "$g_AnnotSV(svtBEDcol)"]"
	} else {
	    set SVtypeFromBED ""
	}
	set SVfromBED "[join $Ls_fromBED "\t"]"
    }
   
    # We have several annotations for 1 SV: 1 by gene.
    #
    # Keep only 1 transcript by gene:
    #   - the one selected by the user with the "-txFile" option
    #   - the one with the most of "bp" of overlapped CDS (=overlappedCDSlength)
    #   - if x transcript with same "bp of overlapped CDS", the one with the most of "bp" of overlapped transcript (=overlapped "UTR, exon, intron" = overlappedTxLength)
    set f [open $splitBedFile]
    while {![eof $f]} {
	set L [gets $f]
	if {$L eq ""} {continue}
	set Ls [split $L "\t"]

	# splitSV (all the line from the bedfile)
	if {$g_AnnotSV(svtBEDcol) ne "-1"} {
	    set SVtype "\t[lindex $Ls "$g_AnnotSV(svtBEDcol)"]"
	} else {
	    set SVtype ""
	}
	set splitSV "[join [lrange $Ls 0 end-10] "\t"]"

	if {$splitSV ne $oldSplitSV} {;# new annotated SV line (all the split lines are done for the oldSV) => we write all information about the oldSV
	    while {$SVfromBED ne $splitSV} {
		# Writing of the "full" SV line (not present in the $splitBedFile file, if not covering a gene)
		lappend L_linesToWrite "[lindex $L_Bed $n]\tfull"
		incr i_toWrite
		if {$i_toWrite > 10000} {
		    WriteTextInFile [join $L_linesToWrite "\n"] "$FullAndSplitBedFile"
		    set L_linesToWrite {}
		    set i_toWrite 0
		}
		# Writing of 1 "split" line...
		if {$SVfromBED eq $oldSplitSV} {
		    set L_genes [lsort -unique $L_genes]		    
		    set g_Lgenes($shortOldSplitSV) [join $L_genes ";"]
		    # ...for each gene overlapped by the SV
		    foreach gene $L_genes {
			set bestCDSlov -1
			set bestCDSlcompl -1
			set bestTxL -1
			set bestAnn ""
			foreach CDSlov $L_overlappedCDSlength($gene) CDSlcompl $L_completeCDSlength($gene) txL $L_overlappedTxLength($gene) ann $L_annot($gene)  {
			    if {$CDSlov > $bestCDSlov} {
				set bestAnn $ann; set bestCDSlov $CDSlov; set bestCDSlcompl $CDSlcompl; set bestTxL $txL; continue
			    }
			    if {$CDSlov eq $bestCDSlov} {
				if {$txL > $bestTxL} {
				    set bestAnn $ann; set bestCDSlov $CDSlov; set bestCDSlcompl $CDSlcompl; set bestTxL $txL; continue
				}
			    }
 			}
			if {$bestCDSlcompl eq 0} {
			    set percentCDS 0
			} else {
			    set percentCDS [expr {$bestCDSlov*100/$bestCDSlcompl}]
			}
			lappend L_linesToWrite "$bestAnn\t$bestCDSlov\t$percentCDS\t$bestTxL\tsplit" 
			incr i_toWrite
		    }
		    set splitSVleft [lindex $Ls 1]
		    set splitSVright [lindex $Ls 2]
		    unset L_annot
		    unset L_overlappedTxLength
		    unset L_overlappedCDSlength
		    unset L_completeCDSlength
		    # check the catch
		    catch {unset Finish}
		    set L_genes {}	     
		    set oldSplitSV "$splitSV"
		    regexp "^\[^\t\]+\t\[^\t\]+\t\[^\t\]+" $oldSplitSV shortOldSplitSV
		}
		incr n
		set Ls_fromBED "[split [lindex $L_Bed $n] "\t"]"
		if {$g_AnnotSV(svtBEDcol) ne "-1"} {
		    set SVtypeFromBED "\t[lindex $Ls_fromBED "$g_AnnotSV(svtBEDcol)"]"
		} else {
		    set SVtypeFromBED ""
		}
		set SVfromBED "[join $Ls_fromBED "\t"]"
	    }
	}
	set txName [lindex $Ls end-4]
	set gene [lindex $Ls end-5]
	lappend L_allGenesOverlapped $gene

	# Look if the annotation line of a user selected transcript has already been registered for this gene
	# If yes, no need to recalculate the txStart, txEnd, CDSstart...
	if {[info exists Finish($gene)]} {continue}

	# No preferred transcript annotation registered at this step, continue the selection with the CDS/Tx lengths
	lappend L_genes $gene
	lappend L_annot($gene) "$L"

	set txStart [lindex $Ls end-8]
	set txEnd [lindex $Ls end-7]
	set CDSstart [lindex $Ls end-3]
	set CDSend [lindex $Ls end-2]
	set exonStarts [lindex $Ls end-1]
	set exonEnds [lindex $Ls end]

	# Calcul of overlappedTxLength:
	if {$txStart > $splitSVleft} {set splitSV_start $txStart} else {set splitSV_start $splitSVleft}
	if {$txEnd < $splitSVright} {set splitSV_end $txEnd} else {set splitSV_end $splitSVright}
	set overlappedTxLength [expr {$splitSV_end-$splitSV_start}]
	lappend L_overlappedTxLength($gene) "$overlappedTxLength"

	# Calcul of completeCDSlength:
	set completeCDSlength 0
	foreach A [split $exonStarts ","] B [split $exonEnds ","] {
	    if {$A eq "" || $B eq ""} {continue}
	    #set txEnd $B
	    # Remove exon(s) not in the CDS:
	    if {$A<$CDSstart} {
		if {$B<$CDSstart} {continue}
		incr completeCDSlength [expr {$B-$CDSstart}]
		continue
	    }
	    if {$B>$CDSend} {
		if {$A>$CDSend} {continue}
		incr completeCDSlength [expr {$CDSend-$A}]
		continue
	    }
	   incr completeCDSlength [expr {$B-$A}]
	}
	lappend L_completeCDSlength($gene) "$completeCDSlength"

	# Calcul of overlappedCDSlength:
	set overlappedCDSlength 0
	foreach A [split $exonStarts ","] B [split $exonEnds ","] {
	    if {$A eq "" || $B eq ""} {continue}
	    # Remove exon(s) not in the CDS:
	    if {$A<$CDSstart} {
		if {$B<$CDSstart} {continue}
		set A $CDSstart
	    }
	    if {$B>$CDSend} {
		if {$A>$CDSend} {continue}
		set B $CDSend
	    }
	    # Intersect coding exon(s) with the SV:
	    if {$splitSV_start<$A} {
		if {$splitSV_end<$A} {break}
		set i $A
	    } else {set i $splitSV_start}
	    if {$splitSV_start>$B} {continue}
	    if {$splitSV_end<$B} {set j $splitSV_end} else {set j $B}
	    incr overlappedCDSlength [expr {$j-$i}]
	}
	lappend L_overlappedCDSlength($gene) "$overlappedCDSlength"

	# Do we have a user selected transcript?
	if {"$txName" ne "" && [lsearch -exact $L_selectedTx "$txName"] ne -1} {
	    set L_annot($gene) {}
	    lappend L_annot($gene) "$L"
	    set L_overlappedTxLength($gene) "$overlappedTxLength"
	    set L_overlappedCDSlength($gene) "$overlappedCDSlength"
	    set L_completeCDSlength($gene) "$completeCDSlength"
	    set Finish($gene) 1
	}

    }
    close $f

    # Treatment for the last SV of the file
    ########################################
    
    # Writing of the "full" SV line (not present in the $splitBedFile file, if not covering a gene)
    # Insertion of the "full length" SV line
    lappend L_linesToWrite "[lindex $L_Bed $n]\tfull" 
    incr i_toWrite
    incr n

    # Treatment of the "split by gene" SV line
    set L_genes [lsort -unique $L_genes]
    set g_Lgenes($shortOldSplitSV) [join $L_genes ";"]
    foreach gene $L_genes {
	set bestCDSlov -1
	set bestCDSlcompl -1
	set bestTxL -1
	set bestAnn ""
	foreach CDSlov $L_overlappedCDSlength($gene) CDSlcompl $L_completeCDSlength($gene) txL $L_overlappedTxLength($gene) ann $L_annot($gene)  {
	    if {$CDSlov > $bestCDSlov} {
		set bestAnn $ann; set bestCDSlov $CDSlov; set bestCDSlcompl $CDSlcompl; set bestTxL $txL; continue
	    }
	    if {$CDSlov eq $bestCDSlov} {
		if {$txL > $bestTxL} {
		    set bestAnn $ann; set bestCDSlov $CDSlov; set bestCDSlcompl $CDSlcompl; set bestTxL $txL; continue
		}
	    }
	}
	if {$bestCDSlcompl eq 0} {
	    set percentCDS 0
	} else {
	    set percentCDS [expr {$bestCDSlov*100/$bestCDSlcompl}]
	}
	lappend L_linesToWrite "$bestAnn\t$bestCDSlov\t$percentCDS\t$bestTxL\tsplit"
	incr i_toWrite
    }
    
    # Insertion of the last "full length" SV lines
    while {$n < [llength $L_Bed]} {
	lappend L_linesToWrite "[lindex $L_Bed $n]\tfull"
	incr i_toWrite
	incr n
    }

    ## Writing
    WriteTextInFile [join $L_linesToWrite "\n"] "$FullAndSplitBedFile"

    ## Delete temporary file
    file delete -force $splitBedFile

    
    ## Annotate the SV bedFile with the regulatory elements files
    ## (Preparation of the phenotype-driven analysis (Exomiser) will be done there)
    set L_allGenesOverlapped [lsort -unique $L_allGenesOverlapped]
    regulatoryElementsAnnotation $L_allGenesOverlapped
    
    ## Memo:
    #
    #  The same SV can be annotated on several genes. Warning: annotations are not necessarily group by gene. Example:
    ##################################################################################################################
    #     1       144676654       144680028       DEL     SGT161364       NBPF8   NR_102405       0       3375    intron5-intron5 144676654       144680028
    #     1       144676654       144680028       DEL     SGT161364       NBPF9   NM_001277444    0       3375    intron8-intron8 144676654       144680028
    #     1       144676654       144680028       DEL     SGT161364       NBPF8   NR_102404       0       3375    intron6-intron6 144676654       144680028
    #
    # It's due to the alternance of NBPF8 and NBPF9 in the genes file:
    ##################################################################
    #     1       144614958       144830407       NBPF8   NR_102405       144830407       144830407       144614958,144618081,144619346,144619882,144621446,144813741,144814679,144815935,144816472,144817965,144823127,144823812,144824704,144825352,144826234,144826932,144827819,144828540,    144615303,144618296,144619419,144620094,144621656,144813844,144814894,144816008,144816678,144818017,144823179,144823985,144824756,144825525,144826286,144827105,144827928,144830407,
    #     1       144614958       145370304       NBPF9   NM_001277444    144615130       145368684       144614958,144615095,144615246,144617149,144618081,144619346,144619882,144621446,144813741,144814679,144815935,144816472,144817965,144821920,144823127,144823812,144824704,144825352,144826234,144826932,144827819,145313333,145314215,145314903,145315790,145368440,    144614998,144615246,144615303,144617252,144618296,144619419,144620094,144621656,144813844,144814894,144816008,144816678,144818017,144822084,144823179,144823985,144824756,144825525,144826286,144827105,144827928,145313506,145314267,145315076,145315899,145370304,
    #     1       144614958       144830407       NBPF8   NR_102404       144830407       144830407       144614958,144617149,144618081,144619346,144619882,144621446,144813741,144814679,144815935,144816472,144817965,144823127,144823812,144824704,144825352,144826234,144826932,144827819,144828540,  144615303,144617252,144618296,144619419,144620094,144621656,144813844,144814894,144816008,144816678,144818017,144823179,144823985,144824756,144825525,144826286,144827105,144827928,144830407,

    return
}


