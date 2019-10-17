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
	set test 0
	foreach col "1000g_event 1000g_AF 1000g_max_AF" {
	    if {[lsearch -exact "$g_AnnotSV(outputColHeader)" $col] ne -1} {set test 1;break}
	}
	if {$test eq 0} {set g_AnnotSV(1000gAnn) 0; return}
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
	##   - 'date'_1000g_SV.sorted.bed file     ; # Header: chr start end 1000g_event 1000g_AF 1000g_max_AF
	set 1000gFileFormatted "$1000gDir/[clock format [clock seconds] -format "%Y%m%d"]_1000g_SV.bed"

	puts "...1000g configuration"

	puts "\t...creation of $1000gFileFormatted"
	puts "\t   (done only once during the first boundaries 1000g annotation)\n"

	set L_TextToWrite {}
	set f [open "|gzip -cd $1000gFileDownloaded"]
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
	    set L_alt [split [lindex $Ls 4] ","]
	    set lengthAlt [llength $L_alt]
	    set i 0
	    foreach alt $L_alt {
		# Consider only the SV (not the SNV/indel in the VCF file)
		##########################################################
		# Example of SV:
		# - Type1: ref="G" and alt="ACTGCTAACGATCCGTTTGCTGCTAACGATCTAACGATCGGGATTGCTAACGATCTCGGG"
		# - Type2: "<INS>", "<DEL>", ...
		# - Type3: complex rearrangements with breakends: "G]17:1584563]"
		if {![regexp "<|\\\[|\\\]" $alt]} {
		    set variantLength [expr {[string length $ref]-[string length $alt]}]
		    if {[expr {abs($variantLength)}]<$g_AnnotSV(SVminSize)} {incr i; continue}; # it is not an SV
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
			set $val [lindex [split [set $val] ","] $i]
			if {[set $val]>$max} {set max [set $val]}
		    }
		}
		if {![info exists AF]} {incr i; continue}
		if {![info exists END]} {
		    # INS:ME (LINE1, ALU or SVA)
		    set END [expr {$pos+1}]
		}
		# In case of multiple alt ("<CN0>,<CN2>"), SVTYPE has a unique value (often: SVTYPE=CNV  else: SVTYPE=DUP).
		# => We replace CNV with a more precise value <CN0> or <CNV2>
		if {$lengthAlt>1} {set SVTYPE [lindex [split [lindex $Ls 4] ","] $i]}
		lappend L_TextToWrite "$chrom\t$pos\t$END\t$SVTYPE\t$AF\t$max"
		incr i
	    }
	}
	close $f
	file delete -force $1000gFileDownloaded

	ReplaceTextInFile [join $L_TextToWrite "\n"] $1000gFileFormatted


	# Sorting of the bedfile:
	# Intersection with very large files can cause trouble with excessive memory usage.
	# A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm.
	regsub -nocase ".bed$" $1000gFileFormatted ".sorted.bed" 1000gFileFormattedAndSorted
	if {[catch {eval exec sort -k1,1 -k2,2n $1000gFileFormatted > $1000gFileFormattedAndSorted} Message]} {
	    puts "-- check1000gFile --"
	    puts "sort -k1,1 -k2,2n $1000gFileFormatted > $1000gFileFormattedAndSorted"
	    puts "$Message"
	    puts "Exit with error"
	    exit 2
	}
	file delete -force $1000gFileFormatted
    }
}




proc 1000gAnnotation {SVchrom SVstart SVend L_i} {

    global g_AnnotSV
    global 1000gText


    set 1000gDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/SVincludedInFt/1000g/$g_AnnotSV(genomeBuild)"
    set 1000gBEDfile [glob -nocomplain "$1000gDir/*_1000g_SV.sorted.bed"]

    if {![info exists 1000gText(DONE)]} {

	# headerOutput "1000g_event 1000g_AF 1000g_max_AF"
	set L_1000gText(Empty) "{} {-1} {-1}"
	foreach i $L_i {
	    lappend 1000gText(Empty) "[lindex $L_1000gText(Empty) $i]"
	}
	set 1000gText(Empty) [join $1000gText(Empty) "\t"]

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

	    set 1000g_max_AF [lindex $Ls end]
	    set 1000g_AF [lindex $Ls end-1]
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
	    # if {![regexp "^(ALU)|(INS)|(LINE1)|(SVA)$" $1000g_event]} {}

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
	    lappend L_1000gAnn_AF($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) $1000g_AF
	    lappend L_1000gAnn_max_AF($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) $1000g_max_AF

	}


	# Loading 1000g final annotation for each SV
	foreach SVtoAnn [array names L_1000gAnn_event] {
	    # Change metrics from "." to ","
	    if {[set g_AnnotSV(metrics)] eq "fr"} {
		regsub -all {\.} $L_1000gAnn_AF($SVtoAnn) "," L_1000gAnn_AF($SVtoAnn)
		regsub -all {\.} $L_1000gAnn_max_AF($SVtoAnn) "," L_1000gAnn_max_AF($SVtoAnn)
	    }
	    # For L_1000gAnn_event($SVtoAnn): we keep a unique list of the events
	    set L_1000gAnn_event($SVtoAnn) [lsort -unique $L_1000gAnn_event($SVtoAnn)]
	    # For L_1000gAnn_AF($SVtoAnn) and L_1000gAnn_max_AF($SVtoAnn): we keep only the maximum value
	    set max_AF -1
	    foreach a $L_1000gAnn_AF($SVtoAnn) {if {$a > $max_AF} {set max_AF $a}}
	    set max2_AF -1
	    foreach a $L_1000gAnn_max_AF($SVtoAnn) {if {$a > $max2_AF} {set max2_AF $a}}

	    set L_1000gText($SVtoAnn) ""
	    lappend L_1000gText($SVtoAnn) "[join $L_1000gAnn_event($SVtoAnn) ";"]"
	    lappend L_1000gText($SVtoAnn) "$max_AF"
	    lappend L_1000gText($SVtoAnn) "$max2_AF"

	    # Keep only the user requested columns (defined in the configfile)
	    set 1000gText($SVtoAnn) ""
	    foreach i $L_i {
		lappend 1000gText($SVtoAnn) "[lindex $L_1000gText($SVtoAnn) $i]"
	    }
	    set 1000gText($SVtoAnn) [join $1000gText($SVtoAnn) "\t"]
	}
	catch {unset L_1000gText}
	set 1000gText(DONE) 1
	file delete -force $tmpFile
    }

    if {[info exist 1000gText($SVchrom,$SVstart,$SVend)]} {
	return $1000gText($SVchrom,$SVstart,$SVend)
    } else {
	return $1000gText(Empty)
    }
}
