############################################################################################################
# AnnotSV 2.5                                                                                              #
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


## - Check and create if necessary the 'date'_1000g_SV.sorted.bed" file.
proc check1000gFile {} {

    global g_AnnotSV


    ## Check if 1000g file has been downloaded then formatted
    #########################################################
    set 1000gDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/SVincludedInFt/1000g/$g_AnnotSV(genomeBuild)"
    set 1000gFileDownloaded [glob -nocomplain "$1000gDir/ALL.wgs.mergedSV*.vcf.gz"]
    set 1000gFileFormattedAndSorted  [glob -nocomplain "$1000gDir/*_1000g_SV.sorted.bed"]

    if {$1000gFileDownloaded eq "" && $1000gFileFormattedAndSorted eq ""} {
	# No 1000g annotation
	set g_AnnotSV(1000gAnn) 0
	return
    } else {
	set g_AnnotSV(1000gAnn) 1
	# Check if the user asked for these annotations in the configfile
	if {[lsearch -exact "$g_AnnotSV(outputColHeader)" "1000g_event"] eq -1} {set g_AnnotSV(1000gAnn) 0; return}
    }

    if {[llength $1000gFileFormattedAndSorted]>1} {
	puts "Several 1000g files exist:"
	puts "$1000gFileFormattedAndSorted"
	puts "Keep only one: [lindex $1000gFileFormattedAndSorted end]\n"
	foreach 1000gF [lrange $1000gFileFormattedAndSorted 0 end-1] {
	    file rename -force $1000gF $1000gF.notused
	}
	return
    }

    if {$1000gFileFormattedAndSorted eq ""} {
	# The downloaded file exist but not the formatted.
	## Create:
	##   -> 'date'_1000g_SV.sorted.bed file     ; # Header: chr start end 1000g_event 1000g_AF 1000g_max_AF
	##   Only "1000g_event" is used for the annotations ("1000g_AF" and "1000g_max_AF" are kept in this file to come back to frequencies if needed by the user).
	set 1000gFileFormatted "$1000gDir/[clock format [clock seconds] -format "%Y%m%d"]_1000g_SV.bed"
	puts "...1000g configuration"

	# Split multiallelic sites into multiple rows
	set 1000gFileTmp "$1000gDir/tmp_1000g_SV.vcf"	
	puts "\t...split multiallelic sites into multiple rows for $1000gFileDownloaded ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	catch {eval exec $g_AnnotSV(bcftools) norm -m -both $1000gFileDownloaded > $1000gFileTmp} Message
	if {[file size $1000gFileTmp] eq 0} {
	    # we continue AnnotSV without splitting the multiallelic sites !
	    puts "\t   -- check1000gFile --"
	    puts "\t   $g_AnnotSV(bcftools) norm -m -both $1000gFileDownloaded > $1000gFileTmp"
	    puts "\t   $1000gFileTmp: file empty."
	    puts "\t   No multiallelic treatment done."
	    file delete -force $1000gFileTmp	   
	} 


	# Creation of $1000gFileFormatted
	puts "\t...creation of $1000gFileFormatted"
	puts "\t   (done only once during the first boundaries 1000g annotation)\n"
	if {[file exists $1000gFileTmp]} {
	    set f [open "$1000gFileTmp"]
	} else {
	    set f [open "| gzip -cd $1000gFileDownloaded"]
	}
	set L_TextToWrite {}
	while {![eof $f]} {
	    set L [gets $f]
	    if {[string index $L 0] eq "#" || $L eq ""} {continue}
	    # Line example:
	    # 1  645710  ALU_umary_ALU_2 A  <INS:ME:ALU>  .  .  AC=35;AF=0.00698882;AFR_AF=0;AMR_AF=0.0072;AN=5008;CS=ALU_umary;EAS_AF=0.0069;EUR_AF=0.0189;MEINFO=AluYa4_5,1,223,-;NS=2504;SAS_AF=0.0041;SITEPOST=0.9998;SVLEN=222;SVTYPE=ALU;TSD=null  GT  0|0  ...
	    set Ls [split $L "\t"]

	    set chrom [lindex $Ls 0]
	    regsub -nocase "chr" $chrom "" chrom
	    set pos [lindex $Ls 1]
 	    set L_infos [split [lindex $Ls 7] ";"]
	    set ref [lindex $Ls 3]
	    set alt [lindex $Ls 4]
	    # Consider only the SV (not the SNV/indel in the VCF file)
	    ##########################################################
	    # Example of SV:
	    # - Type1: ref="G" and alt="ACTGCTAACGATCCGTTTGCTGCTAACGATCTAACGATCGGGATTGCTAACGATCTCGGG"
	    # - Type2: "<INS>", "<DEL>", ...
	    # - Type3: complex rearrangements with breakends: "G]17:1584563]"
	    if {![regexp "<|\\\[|\\\]" $alt]} {
		set variantLength [expr {[string length $ref]-[string length $alt]}]
		if {[expr {abs($variantLength)}]<$g_AnnotSV(SVminSize)} {continue}; # it is not an SV
	    }
	    set max -1
	    catch {unset AF}
	    catch {unset END}
	    set SVTYPE "CNV"
	    foreach inf $L_infos {
		set inf [split $inf "="]
		set val [lindex $inf 0]
		set $val [lindex $inf 1]
		if {[regexp "AF$" $val]} {
		    if {[set $val]>$max} {set max [set $val]}
		}
	    }
	    # Consider only the frequent SV
	    if {![info exists AF] || $AF < 0.005} {continue}
	    if {![info exists END]} {
		# INS:ME (LINE1, ALU or SVA)
		set END [expr {$pos+1}]
	    }
	    # SVTYPE is defined in $L_infos
	    # In case of a single value for "alt", SVTYPE is correctly defined (e.g. SVTYPE=DUP).
	    # In case of multiple alt ("<CN0>,<CN2>"), SVTYPE has a unique value (often: SVTYPE=CNV)
	    #   => We replace CNV with a more precise value <CN0> or <CNV2>
	    if {$SVTYPE eq "CNV"} {set SVTYPE "$alt"}
	    # Sometimes, we have "SVTYPE=DUP" and alt=<CN0>
	    # or "SVTYPE=DEL" and alt=<CN2>
	    # It's due to multiallelic sites. The good value is in "alt"
	    if {[regexp "^<CN" $alt]} {
		set SVTYPE "$alt"
	    }
	    lappend L_TextToWrite "$chrom\t$pos\t$END\t$SVTYPE\t$AF\t$max"
	    incr i
	}
    
	close $f
	file delete -force $1000gFileDownloaded
	file delete -force $1000gFileTmp	   

	ReplaceTextInFile [join $L_TextToWrite "\n"] $1000gFileFormatted


	# Sorting of the bedfile:
	# Intersection with very large files can cause trouble with excessive memory usage.
	# A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm.
	regsub -nocase ".bed$" $1000gFileFormatted ".sorted.bed" 1000gFileFormattedAndSorted
	set sortTmpFile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_sort.tmp.bash"
	ReplaceTextInFile "#!/bin/bash" $sortTmpFile
	WriteTextInFile "# The locale specified by the environment can affects the traditional sort order. We need to use native byte values." $sortTmpFile
	WriteTextInFile "export LC_ALL=C" $sortTmpFile
	WriteTextInFile "sort -k1,1 -k2,2n $1000gFileFormatted > $1000gFileFormattedAndSorted" $sortTmpFile
	file attributes $sortTmpFile -permissions 0755
	if {[catch {eval exec bash $sortTmpFile} Message]} {
	    puts "-- check1000gFile --"
	    puts "sort -k1,1 -k2,2n $1000gFileFormatted > $1000gFileFormattedAndSorted"
	    puts "$Message"
	    puts "Exit with error"
	    exit 2
	}
	file delete -force $sortTmpFile 
	file delete -force $1000gFileFormatted
    }
}



proc 1000gAnnotation {SVchrom SVstart SVend} {

    global g_AnnotSV
    global 1000gText


    set 1000gDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/SVincludedInFt/1000g/$g_AnnotSV(genomeBuild)"
    set 1000gBEDfile [glob -nocomplain "$1000gDir/*_1000g_SV.sorted.bed"]

    if {![info exists 1000gText(DONE)]} {

	# headerOutput "1000g_event"
	set 1000gText(Empty) ""

	# Intersect
	regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".intersect.1000g" tmpFile
	set tmpFile "$g_AnnotSV(outputDir)/[file tail $tmpFile]"
	file delete -force $tmpFile
	if {[catch {exec $g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(fullAndSplitBedFile) -b $1000gBEDfile -wa -wb > $tmpFile} Message]} {
	    puts "-- 1000gAnnotation --"
	    puts "$g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(fullAndSplitBedFile) -b $1000gBEDfile -wa -wb > $tmpFile"
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

	    set 1000g_event [lindex $Ls end-2]
	    set 1000g_end [lindex $Ls end-3]
	    set 1000g_start [lindex $Ls end-4]

	    # Select:
	    # - 1000g that share > XX% length with the SV to annotate
	    # -> doesn't select the INSERTION (ALU, INS, LINE1, SVA: length < 50bp) if overlap < 70%
	    set 1000g_length [expr {$1000g_end-$1000g_start}]
	    set SVtoAnn_length [expr {$SVtoAnn_end-$SVtoAnn_start}]
	    # The 1000g SV is an insertion or a breakpoint
	    if {$1000g_length<1} {
		set 1000g_length 1
	    }
	    # The SV to annotate is an insertion or a breakpoint
	    if {$SVtoAnn_length<1} {
		set SVtoAnn_length 1
	    }

	    # Values of "1000g_event":
	    # ALU, <CN0>, <CN2>, <CN3>, <CN4>, <CN5>, <CN6>, <CN7>, <CN8>, <CN9>, DEL, DEL_ALU, DEL_HERV, DEL_LINE1, DEL_SVA, DUP, INS, INV, LINE1, SVA

	    # we keep only SV with > 70% (default) overlap with the SV to annotate
	    if {$SVtoAnn_start < $1000g_start} {
		set overlap_start $1000g_start
	    } else {
		set overlap_start $SVtoAnn_start
	    }
	    if {$SVtoAnn_end < $1000g_end} {
		set overlap_end $SVtoAnn_end
	    } else {
		set overlap_end $1000g_end
	    }
	    set overlap_length [expr {$overlap_end - $overlap_start}]

	    # Keeping only 1000g respecting the overlaps (reciprocal or not reciprocal)
	    if {[expr {$overlap_length*100.0/$SVtoAnn_length}] < $g_AnnotSV(overlap)} {continue}
	    if {$g_AnnotSV(reciprocal) eq "yes"} {
		if {[expr {$overlap_length*100.0/$1000g_length}] < $g_AnnotSV(overlap)} {continue}
	    }

	    lappend L_1000gAnn_event($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) $1000g_event

	}


	# Loading 1000g final annotation for each SV
	foreach SVtoAnn [array names L_1000gAnn_event] {
	    # For L_1000gAnn_event($SVtoAnn): we keep a unique list of the events
	    set L_1000gAnn_event($SVtoAnn) [lsort -unique $L_1000gAnn_event($SVtoAnn)]
	    set 1000gText($SVtoAnn) "[join $L_1000gAnn_event($SVtoAnn) "/"]"
	}
	set 1000gText(DONE) 1
	file delete -force $tmpFile
    }

    if {[info exist 1000gText($SVchrom,$SVstart,$SVend)]} {
	return $1000gText($SVchrom,$SVstart,$SVend)
    } else {
	return $1000gText(Empty)
    }
}
