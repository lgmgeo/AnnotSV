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


## - Check and create if necessary the "promoter.XXXbp.sorted.bed" file.
proc checkPromoterFile {} {

    global g_AnnotSV

    ## Check if the promoters file has been formatted
    #################################################
    set refgeneDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/RefGene/$g_AnnotSV(genomeBuild)"
    set promoterDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/FtIncludedInSV/Promoter/$g_AnnotSV(genomeBuild)"

    if {![file exists $promoterDir]} {file mkdir $promoterDir}

    set refgeneFileFormatted "$refgeneDir/refGene.sorted.bed"
    set promoterFormatted "$promoterDir/promoter.$g_AnnotSV(promoterSize)bp.sorted.bed"
    set g_AnnotSV(promAnn) 1
    if {[lsearch -exact "$g_AnnotSV(outputColHeader)" "promoters"] eq -1} {
	set g_AnnotSV(promAnn) 0; return
    }

    if {![file exists $promoterFormatted]} {
	
	# Creation of the promoterFormatted 
	foreach L [LinesFromFile $refgeneFileFormatted] {
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
		lappend L_allProm "$chrom\t$prom\t[join [lsort -unique $genesProm($chrom\t$prom)] ","]"
	    }
	}
	WriteTextInFile [join $L_allProm "\n"] $promoterFormatted.tmp
	# Sorting of the bedfile:
	# Intersection with very large files can cause trouble with excessive memory usage.
	# A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm. 
	if {[catch {eval exec sort -k1,1 -k2,2n $promoterFormatted.tmp > $promoterFormatted} Message]} {
	    puts "-- checkPromoterFile --"
	    puts "sort -k1,1 -k2,2n $promoterFormatted.tmp > $promoterFormatted"
	    puts "$Message"
	    puts "Exit with error"
	    exit 2
	}
	file delete -force $promoterFormatted.tmp

    }
}




proc promoterAnnotation {SVchrom SVstart SVend} {

    global g_AnnotSV
    global promoterText


    # headerOutput "promoters"
    set promoterText(Empty) ""


    set promoterDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/FtIncludedInSV/Promoter/$g_AnnotSV(genomeBuild)"
    set promoterBEDfile "$promoterDir/promoter.$g_AnnotSV(promoterSize)bp.sorted.bed"
    
    if {![info exists promoterText(DONE)]} {
	
	# Intersect
	regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".intersect.promoter" tmpFile
	set tmpFile "$g_AnnotSV(outputDir)/[file tail $tmpFile]"
	file delete -force $tmpFile	
	if {[catch {exec $g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(fullAndSplitBedFile) -b $promoterBEDfile -wa -wb > $tmpFile} Message]} {
	    puts "-- promoterAnnotation --"
	    puts "$g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(fullAndSplitBedFile) -b $promoterBEDfile -wa -wb > $tmpFile"
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

	    set promGenes  [lindex $Ls end]
	    set prom_end [lindex $Ls end-1]
	    set prom_start [lindex $Ls end-2]

	    # Select:
	    # - promoters that share > XX% length with the SV to annotate
	    # - promoters covering the SV that is an insertion
	    set prom_length [expr {$prom_end-$prom_start}]
	    set SVtoAnn_length [expr {$SVtoAnn_end - $SVtoAnn_start}]
	    if {$SVtoAnn_length>1} {
		# SV to annotate is not an insertion
		if {$SVtoAnn_start < $prom_start} {
		    set overlap_start $prom_start
		} else {
		    set overlap_start $SVtoAnn_start
		}
		if {$SVtoAnn_end < $prom_end} {
		    set overlap_end $SVtoAnn_end
		} else {
		    set overlap_end $prom_end
		}
		set overlap_length [expr {$overlap_end - $overlap_start}]
		
		# Keeping only promoter with > 70% (default) overlap with the SV
		if {[expr {$overlap_length*100.0/$prom_length}] < $g_AnnotSV(overlap)} {continue}
		# No reciprocal overlap used for promoter (it would not make sense)
	    }
	
	    lappend L_genes($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) {*}[split $promGenes ","]

	}
	

	# Loading promoter final annotation for each SV
	foreach SVtoAnn [array names L_genes] {
	    set promoterText($SVtoAnn) "[join [lsort -unique $L_genes($SVtoAnn)] "/"]"

	}

	set promoterText(DONE) 1	
	file delete -force $tmpFile
    }
    
    if {[info exist promoterText($SVchrom,$SVstart,$SVend)]} {
	return $promoterText($SVchrom,$SVstart,$SVend)
    } else {
	return $promoterText(Empty)
    }
}
