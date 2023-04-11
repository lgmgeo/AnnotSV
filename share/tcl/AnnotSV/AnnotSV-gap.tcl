############################################################################################################
# AnnotSV 3.3.4                                                                                            #
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


## - Check if the downloaded "Gap.bed" file or the formatted "*_Gap.sorted.bed" file exist
proc checkGapFile {} {

    global g_AnnotSV

    ## Check if the gap file has been downloaded then formatted
    #############################################################
    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/BreakpointsAnnotations"
    set gapFileDownloaded [glob -nocomplain "$extannDir/Gap/$g_AnnotSV(genomeBuild)/Gap.bed"]   
    set gapFileFormatted [glob -nocomplain "$extannDir/Gap/$g_AnnotSV(genomeBuild)/*_Gap.sorted.bed"]   

    if {$gapFileDownloaded eq "" && $gapFileFormatted eq ""} {
	# No Gap annotation
	set g_AnnotSV(gapAnn) 0
	return
    } else {
	# Gap annotation
	set g_AnnotSV(gapAnn) 1
	# Check if the user asked for these annotations in the configfile
	set test 0
	foreach col "Gap_left Gap_right" {
	    if {[lsearch -exact "$g_AnnotSV(outputColHeader)" $col] ne -1} {set test 1;break}
	}
	if {$test eq 0} {set g_AnnotSV(gapAnn) 0; return}
    }

    if {[llength $gapFileFormatted]>1} {
	puts "Several gap files exist:"
	puts "$gapFileFormatted"
	puts "Keep only one: [lindex $gapFileFormatted end]\n"
	foreach gap [lrange $gapFileFormatted 0 end-1] {
	    file rename -force $gap $gap.notused
	}
	return
    } 

    if {$gapFileFormatted eq ""} {
	# The downloaded file exist but not the formatted.
	set gapFileFormatted "$extannDir/Gap/$g_AnnotSV(genomeBuild)/[clock format [clock seconds] -format "%Y%m%d"]_Gap.sorted.bed"
	puts "\t...creation of the \"[file tail $gapFileFormatted]\" file ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])\n\t   (done only once)"
	ReplaceTextInFile "#Chrom\tStart\tEnd\tGap" $gapFileFormatted.tmp

	set f [open $gapFileDownloaded]
	set i 0
	set L_Text {}
	while {![eof $f]} {
	    set L [gets $f]
	    if {$L eq ""} {continue}
            if {[regexp "^#" $L]} {continue}
	    
            set Ls [split $L "\t"]
            regsub "chr" [lindex $Ls 0] "" chrom
            set start [lindex $Ls 1]
            set end [lindex $Ls 2]
            set IDgap "${chrom}:${start}-${end}"

	    lappend L_Text "$chrom\t$start\t$end\t$IDgap"
	    if {$i>500000} {
		WriteTextInFile [join $L_Text "\n"] $gapFileFormatted.tmp
		set L_Text {}
		set i 0
	    }
            incr i
	}
	WriteTextInFile [join $L_Text "\n"] $gapFileFormatted.tmp

	# Sorting of the bedfile:
	# Intersection with very large files can cause trouble with excessive memory usage.
	# A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm. 
	set sortTmpFile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_sort.tmp.bash"
	ReplaceTextInFile "#!/bin/bash" $sortTmpFile
	WriteTextInFile "# The locale specified by the environment can affects the traditional sort order. We need to use native byte values." $sortTmpFile
	WriteTextInFile "export LC_ALL=C" $sortTmpFile
	WriteTextInFile "sort -k1,1 -k2,2n $gapFileFormatted.tmp > $gapFileFormatted" $sortTmpFile
	file attributes $sortTmpFile -permissions 0755
	if {[catch {eval exec bash $sortTmpFile} Message]} {
	    puts "-- checkGapFile --"
	    puts "sort -k1,1 -k2,2n $gapFileFormatted.tmp > $gapFileFormatted"
	    puts "$Message"
	    puts "Exit with error"
	    exit 2
	}
	file delete -force $sortTmpFile 
	file delete -force $gapFileFormatted.tmp
	file delete -force $gapFileDownloaded

	close $f
	
	return
    }
}



proc GapAnnotation {BreakpointChrom BreakpointPos} {

    global g_AnnotSV
    global g_Gap_type
    global g_Gap_coord
    global g_Gap


    # headerOutput "Gap_left" or "Gap_right"
    # (GapAnnotation is executed for each breakpoint)
    set g_Gap(Empty) ""

    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/BreakpointsAnnotations"
    set gapFileFormatted [glob -nocomplain "$extannDir/Gap/$g_AnnotSV(genomeBuild)/*_Gap.sorted.bed"]

    if {![info exists g_Gap(DONE)]} {
	
	# Creation of a bedfile with all breakpoints +/- 100 bp:
	# => done during the GC content annotation (if "GC content" annotation is in the output (user defined))  
	# FORMAT:  chrom   'breakpoint-100'   'breakpoint+100'   'breakpoint'
	regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".breakpoints.bed" tmpBreakpointsFile
	set tmpBreakpointsFile "$g_AnnotSV(outputDir)/[file tail $tmpBreakpointsFile]"
	if {![file exists $tmpBreakpointsFile]} {
	    file delete -force $tmpBreakpointsFile.tmp
	    set i 0
	    set L_LinesToWrite {}
	    set f [open $g_AnnotSV(bedFile)]
	    while {![eof $f]} {
		set L [gets $f]
		if {$L eq ""} {continue}
		set Ls [split $L "\t"]
		set chrom [lindex $Ls 0]
		set start [lindex $Ls 1]
		set end [lindex $Ls 2]
		set start_breakpoint1 [expr $start-100]; if {$start_breakpoint1 < 0} {set start_breakpoint1 0}
		set end_breakpoint1 [expr $start+100]  ; if {$end_breakpoint1 < 0} {set end_breakpoint1 0}
		set start_breakpoint2 [expr $end-100]  ; if {$start_breakpoint2 < 0} {set start_breakpoint2 0}
		set end_breakpoint2 [expr $end+100]    ; if {$end_breakpoint2 < 0} {set end_breakpoint2 0}
		lappend L_LinesToWrite "$chrom\t$start_breakpoint1\t$end_breakpoint1\t$start"
		lappend L_LinesToWrite "$chrom\t$start_breakpoint2\t$end_breakpoint2\t$end"
		incr i
		if {$i>100000} {
		    WriteTextInFile [join [lsort -unique $L_LinesToWrite] "\n"] $tmpBreakpointsFile.tmp
		    set L_LinesToWrite {}
		    set i 0
		}
	    }
	    close $f
	    WriteTextInFile [join [lsort -unique $L_LinesToWrite] "\n"] $tmpBreakpointsFile.tmp
	    unset L_LinesToWrite
	    
	    # Sorting of the bedfile:
	    # Intersection with very large files can cause trouble with excessive memory usage.
	    # A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm. 
	    set sortTmpFile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_sort.tmp.bash"
	    ReplaceTextInFile "#!/bin/bash" $sortTmpFile
	    WriteTextInFile "# The locale specified by the environment can affects the traditional sort order. We need to use native byte values." $sortTmpFile
	    WriteTextInFile "export LC_ALL=C" $sortTmpFile
	    WriteTextInFile "sort -k1,1 -k2,2n $tmpBreakpointsFile.tmp > $tmpBreakpointsFile" $sortTmpFile
	    file attributes $sortTmpFile -permissions 0755
	    if {[catch {eval exec bash $sortTmpFile} Message]} {
		puts "-- GapAnnotation, sort --"
		puts "sort -k1,1 -k2,2n $tmpBreakpointsFile.tmp > $tmpBreakpointsFile"
		puts "$Message"
		puts "Exit with error"
		exit 2
	    }
	    file delete -force $sortTmpFile 
	    file delete -force $tmpBreakpointsFile.tmp
	    
	}

	# Intersect
	regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".intersect.Gap.bed" tmpGapFile
	set tmpGapFile "$g_AnnotSV(outputDir)/[file tail $tmpGapFile]"
	file delete -force $tmpGapFile
	if {[catch {exec $g_AnnotSV(bedtools) intersect -sorted -a $tmpBreakpointsFile -b $gapFileFormatted -wa -wb > $tmpGapFile} Message]} {
	    if {[catch {exec $g_AnnotSV(bedtools) intersect -a $tmpBreakpointsFile -b $gapFileFormatted -wa -wb > $tmpGapFile} Message]} {
		puts "-- GapAnnotation, intersect --"
		puts "$g_AnnotSV(bedtools) intersect -sorted -a $tmpBreakpointsFile -b $gapFileFormatted -wa -wb > $tmpGapFile"
		puts "$Message"
		puts "Exit with error"
		exit 2
	    }
	}

	# Loading g_Gap for each SV breakpoint
	set f [open $tmpGapFile]
	while {![eof $f]} {
	    set L [gets $f]
	    if {$L eq ""} {continue}
	    set Ls [split $L "\t"]
	    
	    set chrom [lindex $Ls 0]
	    set breakpoint [lindex $Ls 3]
	    #set gapchrom [lindex $Ls 6]
	    #set gapStart [lindex $Ls 5]
	    #set gapEnd [lindex $Ls 6]
	    #set gapCoord [lindex $Ls 7]		    
	    lappend g_Gap_coord($chrom,$breakpoint) "[lindex $Ls 7]"
	}
	
	set g_Gap(DONE) 1	
	file delete -force $tmpGapFile
    }
    
    if {[info exist g_Gap_coord($BreakpointChrom,$BreakpointPos)]} {
	return "[join $g_Gap_coord($BreakpointChrom,$BreakpointPos) ";"]"
    } else {
	return $g_Gap(Empty)
    }
}
