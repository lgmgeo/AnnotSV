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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                                             #
# GNU General Public License for more details.                                                             #
#                                                                                                          #
# You should have received a copy of the GNU General Public License                                        #
# along with this program; If not, see <http://www.gnu.org/licenses/>.                                     #
############################################################################################################


## - Check and create if necessary the "refGene.sorted.bed" file.
proc checkRefGeneFile {} {

    global g_AnnotSV

    ## Check if the refGene file has been downloaded then formatted
    ##############################################################
    set refgeneDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/RefGene/$g_AnnotSV(genomeBuild)"

    set refgeneFileDownloaded "[glob -nocomplain $refgeneDir/refGene.txt.gz]"
    set refgeneFileFormatted "[glob -nocomplain $refgeneDir/refGene.sorted.bed]"

    if {$refgeneFileDownloaded eq "" && $refgeneFileFormatted eq ""} {
	puts "############################################################################"
	puts "\"$refgeneFileDownloaded\" doesn't exist"
	puts "Please check your install - Exit with error."
	puts "############################################################################"
	exit 2
    }

    if {$refgeneFileFormatted eq ""} {

	## Delete promoters files (need to be updated after the creation of new refGene file)
	##################################################################################### 
	set promoterDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/FtIncludedInSV/Promoter/$g_AnnotSV(genomeBuild)"
	foreach promFile [glob -nocomplain "$promoterDir/promoter.*bp.sorted.bed"] {
	    file delete -force $promFile
	}

	## - Create the "refGene.sorted.bed".
	#####################################
	set refgeneFileFormatted "$refgeneDir/refGene.sorted.bed"
	puts "...creation of $refgeneFileFormatted ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	puts "\t   (done only once during the first annotation)\n"

	# Removing non-standard contigs (other than the standard 1-22,X,Y,MT) and sorting the file in karyotypic order
	# -> create L_refGeneTXTsorted

	## Save the line of the refGeneTXT by chromosome (L_lines($chrom))
	foreach L [LinesFromGZFile $refgeneFileDownloaded] {
	    set Ls [split $L "\t"]
	    set chrom [lindex $Ls 2]
	    lappend L_lines($chrom) "$L"
	}
	## Sorting in karyotypic order
	set L_refGeneTXTsorted {}
	foreach val {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M MT} {
	    if {![info exists L_lines(chr$val)]} {continue}
	    lappend L_refGeneTXTsorted {*}"[lsort -command AscendingSortOnElement4 [lsort -command AscendingSortOnElement5 $L_lines(chr$val)]]"
	}
	# Creation of the $refgeneFileFormatted.
	# Chromosome nomenclature used: 1-22,X,Y,M,T (without "chr")
	## INPUT:    #bin name chrom strand txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds score name2 cdsStartStat cdsEndStat exonFrames
	## OUTPUT:   chrom txStart txEnd strand name2 name cdsStart cdsEnd exonStarts exonEnds
	# WARNING : NR_* are for non coding RNA. However, cdsStart=cdsEnd, => CDSlength=1
	foreach L $L_refGeneTXTsorted {
	    set Ls [split $L "\t"]
	    regsub "chr" [lindex $Ls 2] "" chrom
	    set line "$chrom\t[lindex $Ls 4]\t[lindex $Ls 5]\t[lindex $Ls 3]\t[lindex $Ls 12]\t[lindex $Ls 1]\t[lindex $Ls 6]\t[lindex $Ls 7]\t[lindex $Ls 9]\t[lindex $Ls 10]"
	    if {![info exists infos($line)]} {
		WriteTextInFile $line $refgeneFileFormatted.tmp
		set infos($line) 1
	    }
	}
	file delete -force $refgeneFileDownloaded
	# Sorting of the bedfile:
	# Intersection with very large files can cause trouble with excessive memory usage.
	# A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm.
	if {[catch {eval exec sort -k1,1 -k2,2n $refgeneFileFormatted.tmp > $refgeneFileFormatted} Message]} {
	    puts "-- checkRefGeneFile --"
	    puts "sort -k1,1 -k2,2n $refgeneFileFormatted.tmp > $refgeneFileFormatted"
	    puts "$Message"
	    puts "Exit with error"
	    exit 2
	}
	file delete -force $refgeneFileFormatted.tmp

    }

    # DISPLAY:
    ##########
    set g_AnnotSV(refGene) $refgeneFileFormatted

}



## Annotate the SV bedFile with the refGene file.
## Keep only 1 NM annotation by gene:
##   - the one selected by the user with the "-txFile" option
##   - the one with the most of "bp from CDS" (=CDSlength)
##   - if x NM with same "bp from CDS", the one with the most of "bp from UTR, exon, intron" (=txLength)
##
## Creation of tmpFullAndSplitBedFile ($g_AnnotSV(outputDir)/$g_AnnotSV(outputFile).tmp)
## -> formatted and sorted
proc refGeneAnnotation {} {

    global g_AnnotSV
    global g_Lgenes


    # Bedfile should be sorted and should not have "chr" in the first column
    ########################################################################
    # Removing non-standard contigs (other than the standard 1-22,X,Y,MT) and sorting the file in karyotypic order

    set f [open $g_AnnotSV(bedFile)]
    set test 0
    while {![eof $f]} {
	set L [gets $f]
	if {$L eq ""} {continue}
	if {[regsub "chr" [lindex $L 0] "" chrom] ne 0} {
	    lappend L_Text($chrom) "[string range $L 3 end]"
	} else {
	    lappend L_Text($chrom) "$L"
	}
    }
    close $f
    ## Writing the bedfile (not sorted)
    regsub -nocase ".bed$" $g_AnnotSV(bedFile) ".formatted.sorted.bed" newBed
    set newBed "$g_AnnotSV(outputDir)/[file tail $newBed]"
    file delete -force $newBed
    foreach chrom {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M MT} {
	if {![info exists L_Text($chrom)]} {continue}
	WriteTextInFile [join $L_Text($chrom) "\n"] $newBed.tmp
	unset L_Text($chrom)
    }
    # Sorting of the bedfile:
    # Intersection with very large files can cause trouble with excessive memory usage.
    # A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm.
    if {[catch {eval exec sort -k1,1 -k2,2n $newBed.tmp >> $newBed} Message]} {
	puts "-- refGeneAnnotation --"
	puts "sort -k1,1 -k2,2n $newBed.tmp >> $newBed"
	puts "$Message"
	    puts "Exit with error"
	exit 2
    }
    file delete -force $newBed.tmp
    set g_AnnotSV(bedFile) $newBed

    # Used for the insertion of the "full/split" information
    set L_Bed [LinesFromFile $newBed]


    # OUTPUT:
    ###############
    set tmpFullAndSplitBedFile "$g_AnnotSV(outputDir)/$g_AnnotSV(outputFile).tmp"
    set linesSplitByGene_File "$tmpFullAndSplitBedFile.tmp"
    file delete -force $tmpFullAndSplitBedFile
    file delete -force $linesSplitByGene_File


    # Intersect of the input SV bedfile with the genes annotations
    ##############################################################
    ## -> Creation of the split annotations
    if {[catch {exec $g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(bedFile) -b $g_AnnotSV(refGene) -wa -wb > $linesSplitByGene_File} Message]} {
	puts "-- refGeneAnnotation --"
	puts "$g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(bedFile) -b $g_AnnotSV(refGene) -wa -wb > $linesSplitByGene_File"
	puts "$Message"
	puts "Exit with error"
	exit 2
    }
    if {[file size $linesSplitByGene_File] eq 0} {
	puts "\tno intersection between SV and gene annotation"
	set n 0
	while {$n < [llength $L_Bed]} {
	    WriteTextInFile "[lindex $L_Bed $n]\tfull" "$tmpFullAndSplitBedFile"
	    incr n
	}
	## Delete temporary file
	file delete -force $linesSplitByGene_File
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

    # Parse
    ###############
    set L "[FirstLineFromFile $linesSplitByGene_File]"
    set Ls [split $L "\t"]
    set splitSVleft [lindex $Ls 1]
    set splitSVright [lindex $Ls 2]
    # SV/oldSV defined with: "chrom, start, end":
    set oldSplitSV "[join [lrange $Ls 0 2] "\t"]"
    set L_genes {}
    set n 0

    # Keep only 1 NM by gene:
    #   - the one selected by the user with the "-txFile" option
    #   - the one with the most of "bp from CDS" (=CDSlength)
    #   - if x NM with same "bp from CDS", the one with the most of "bp from UTR, exon, intron" (=txLength)
    # We have several annotations for 1 SV: 1 by gene.
    set f [open $linesSplitByGene_File]
    while {![eof $f]} {

	set L [gets $f]
	if {$L eq ""} {continue}
	set Ls [split $L "\t"]
	set splitSV "[join [lrange $Ls 0 2] "\t"]"
	if {$splitSV ne $oldSplitSV} {
	    while {![regexp "^$splitSV" [lindex $L_Bed $n]]} {
		# Insertion of the "full length" SV line (not present in the $linesSplitByGene_File file, if not covering a gene)
		WriteTextInFile "[lindex $L_Bed $n]\tfull" "$tmpFullAndSplitBedFile"
		if {[join [lrange [lindex $L_Bed $n] 0 2] "\t"] eq $oldSplitSV} {
		    # Treatment of the "split by gene" SV line
		    set L_genes [lsort -unique $L_genes]
		    set g_Lgenes($oldSplitSV) [join $L_genes "/"]
		    foreach gene $L_genes {
			set bestCDSl -1
			set bestTxL -1
			set bestAnn ""
			foreach CDSl $L_CDSlength($gene) txL $L_txLength($gene) ann $L_annot($gene)  {
			    if {$CDSl > $bestCDSl} {
				set bestAnn $ann; set bestCDSl $CDSl; set bestTxL $txL; continue
			    }
			    if {$CDSl eq $bestCDSl} {
				if {$txL > $bestTxL} {
				    set bestAnn $ann; set bestCDSl $CDSl; set bestTxL $txL; continue
				}
			    }
 			}
			WriteTextInFile "$bestAnn\t$bestCDSl\t$bestTxL\tsplit" "$tmpFullAndSplitBedFile"
		    }
		    set splitSVleft [lindex $Ls 1]
		    set splitSVright [lindex $Ls 2]
		    unset L_annot
		    unset L_txLength
		    unset L_CDSlength
			# check the catch
		    catch {unset Finish}
		    set L_genes {}	     
		    set oldSplitSV "$splitSV"
		}
		incr n
	    }
	}
	set txName [lindex $Ls end-4]
	set gene [lindex $Ls end-5]

	# Look if the annotation line of a user selected transcript has already been registered for this gene
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

	# Calcul of txLength:
	if {$txStart > $splitSVleft} {set start $txStart} else {set start $splitSVleft}
	if {$txEnd < $splitSVright} {set end $txEnd} else {set end $splitSVright}
	set txLength [expr {$end-$start+1}]
	lappend L_txLength($gene) "$txLength"

	# Calcul of CDSlength:
	set CDSlength 0
	foreach A [split $exonStarts ","] B [split $exonEnds ","] {
	    if {$A eq "" || $B eq ""} {continue}
	    set txEnd $B
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
	    if {$start<$A} {
		if {$end<$A} {break}
		set i $A
	    } else {set i $start}
	    if {$start>$B} {continue}
	    if {$end<$B} {set j $end} else {set j $B}
	    incr CDSlength [expr {$j-$i+1}]
	}
	lappend L_CDSlength($gene) "$CDSlength"

	# Do we have a user selected transcript?
	if {"$txName" ne "" && [lsearch -exact $L_selectedTx "$txName"] ne -1} {
	    set L_annot($gene) {}
	    lappend L_annot($gene) "$L"
	    set L_txLength($gene) "$txLength"
	    set L_CDSlength($gene) "$CDSlength"
	    set Finish($gene) 1
	}

    }
    close $f

    # Treatment for the last SV of the file
    ########################################
    # Insertion of the "full length" SV line
    while {![regexp "^$splitSV" [lindex $L_Bed $n]]} {
	WriteTextInFile "[lindex $L_Bed $n]\tfull" "$tmpFullAndSplitBedFile"
	incr n
    }
    WriteTextInFile "[lindex $L_Bed $n]\tfull" "$tmpFullAndSplitBedFile"
    incr n
    # Treatment of the "split by gene" SV line
    set L_genes [lsort -unique $L_genes]
    set g_Lgenes($oldSplitSV) [join $L_genes "/"]
    foreach gene $L_genes {
	set bestCDSl -1
	set bestTxL -1
	set bestAnn ""
	foreach CDSl $L_CDSlength($gene) txL $L_txLength($gene) ann $L_annot($gene)  {
	    if {$CDSl > $bestCDSl} {
		set bestAnn $ann; set bestCDSl $CDSl; set bestTxL $txL; continue
	    }
	    if {$CDSl eq $bestCDSl} {
		if {$txL > $bestTxL} {
		    set bestAnn $ann; set bestCDSl $CDSl; set bestTxL $txL; continue
		}
	    }
	}
	WriteTextInFile "$bestAnn\t$bestCDSl\t$bestTxL\tsplit" "$tmpFullAndSplitBedFile"
    }
    # Insertion of the last "full length" SV lines
    while {$n < [llength $L_Bed]} {
	WriteTextInFile "[lindex $L_Bed $n]\tfull" "$tmpFullAndSplitBedFile"
	incr n
    }

    ## Delete temporary file
    file delete -force $linesSplitByGene_File


    ## Memo:
    ## The same SV can be annotated on several genes. Warning: annotations are not necessarily group by gene. Example:
    #     1       144676654       144680028       DEL     SGT161364       NBPF8   NR_102405       0       3375    intron5-intron5 144676654       144680028
    #     1       144676654       144680028       DEL     SGT161364       NBPF9   NM_001277444    0       3375    intron8-intron8 144676654       144680028
    #     1       144676654       144680028       DEL     SGT161364       NBPF8   NR_102404       0       3375    intron6-intron6 144676654       144680028
    # It's due to the alternance of NBPF8 and NBPF9 in the refGene file:
    #     1       144614958       144830407       NBPF8   NR_102405       144830407       144830407       144614958,144618081,144619346,144619882,144621446,144813741,144814679,144815935,144816472,144817965,144823127,144823812,144824704,144825352,144826234,144826932,144827819,144828540,    144615303,144618296,144619419,144620094,144621656,144813844,144814894,144816008,144816678,144818017,144823179,144823985,144824756,144825525,144826286,144827105,144827928,144830407,
    #     1       144614958       145370304       NBPF9   NM_001277444    144615130       145368684       144614958,144615095,144615246,144617149,144618081,144619346,144619882,144621446,144813741,144814679,144815935,144816472,144817965,144821920,144823127,144823812,144824704,144825352,144826234,144826932,144827819,145313333,145314215,145314903,145315790,145368440,    144614998,144615246,144615303,144617252,144618296,144619419,144620094,144621656,144813844,144814894,144816008,144816678,144818017,144822084,144823179,144823985,144824756,144825525,144826286,144827105,144827928,145313506,145314267,145315076,145315899,145370304,
    #     1       144614958       144830407       NBPF8   NR_102404       144830407       144830407       144614958,144617149,144618081,144619346,144619882,144621446,144813741,144814679,144815935,144816472,144817965,144823127,144823812,144824704,144825352,144826234,144826932,144827819,144828540,  144615303,144617252,144618296,144619419,144620094,144621656,144813844,144814894,144816008,144816678,144818017,144823179,144823985,144824756,144825525,144826286,144827105,144827928,144830407,

    return
}
