############################################################################################################
# AnnotSV 3.2                                                                                              #
#                                                                                                          #
# AnnotSV: An integrated tool for Structural Variations annotation and ranking                             #
#                                                                                                          #
# Copyright (C) 2017-2022 Veronique Geoffroy (veronique.geoffroy@inserm.fr)                                #
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


## - Check if the downloaded "SegDup.bed" file or the formatted "*_SegDup.sorted.bed" file exist
proc checkSegDupFile {} {

    global g_AnnotSV

    ## Check if the SegDup file has been downloaded then formatted
    #############################################################
    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/BreakpointsAnnotations"
    set segdupFileDownloaded [glob -nocomplain "$extannDir/SegDup/$g_AnnotSV(genomeBuild)/SegDup.bed"]   
    set segdupFileFormatted [glob -nocomplain "$extannDir/SegDup/$g_AnnotSV(genomeBuild)/*_SegDup.sorted.bed"]   

    if {$segdupFileDownloaded eq "" && $segdupFileFormatted eq ""} {
	# No SegDup annotation
	set g_AnnotSV(segdupAnn) 0
	return
    } else {
	# SegDup annotation
	set g_AnnotSV(segdupAnn) 1
	# Check if the user asked for these annotations in the configfile
	set test 0
	foreach col "SegDup_left SegDup_right" {
	    if {[lsearch -exact "$g_AnnotSV(outputColHeader)" $col] ne -1} {set test 1;break}
	}
	if {$test eq 0} {set g_AnnotSV(segdupAnn) 0; return}
    }

    if {[llength $segdupFileFormatted]>1} {
	puts "Several segmental duplication files exist:"
	puts "$segdupFileFormatted"
	puts "Keep only one: [lindex $segdupFileFormatted end]\n"
	foreach segdup [lrange $segdupFileFormatted 0 end-1] {
	    file rename -force $segdup $segdup.notused
	}
	return
    } 

    if {$segdupFileFormatted eq ""} {
	# The downloaded file exist but not the formatted.
	set segdupFileFormatted "$extannDir/SegDup/$g_AnnotSV(genomeBuild)/[clock format [clock seconds] -format "%Y%m%d"]_SegDup.sorted.bed"
	puts "\t...creation of the \"[file tail $segdupFileFormatted]\" file ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])\n\t   (done only once)"
	ReplaceTextInFile "#Chrom\tStart\tEnd\tSegDup" $segdupFileFormatted.tmp

	set f [open $segdupFileDownloaded]
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
            set IDsegdup "${chrom}:${start}-${end}"

	    lappend L_Text "$chrom\t$start\t$end\t$IDsegdup"
	    if {$i>500000} {
		WriteTextInFile [join $L_Text "\n"] $segdupFileFormatted.tmp
		set L_Text {}
		set i 0
	    }
            incr i
	}
	WriteTextInFile [join $L_Text "\n"] $segdupFileFormatted.tmp

	# Sorting of the bedfile:
	# Intersection with very large files can cause trouble with excessive memory usage.
	# A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm. 
	set sortTmpFile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_sort.tmp.bash"
	ReplaceTextInFile "#!/bin/bash" $sortTmpFile
	WriteTextInFile "# The locale specified by the environment can affects the traditional sort order. We need to use native byte values." $sortTmpFile
	WriteTextInFile "export LC_ALL=C" $sortTmpFile
	WriteTextInFile "sort -k1,1 -k2,2n $segdupFileFormatted.tmp > $segdupFileFormatted" $sortTmpFile
	file attributes $sortTmpFile -permissions 0755
	if {[catch {eval exec bash $sortTmpFile} Message]} {
	    puts "-- checkSegDupFile --"
	    puts "sort -k1,1 -k2,2n $segdupFileFormatted.tmp > $segdupFileFormatted"
	    puts "$Message"
	    puts "Exit with error"
	    exit 2
	}
	file delete -force $sortTmpFile 
	file delete -force $segdupFileFormatted.tmp
	file delete -force $segdupFileDownloaded

	close $f
	
	return
    }
}



proc SegDupAnnotation {BreakpointChrom BreakpointPos} {

    global g_AnnotSV
    global g_SegDup_type
    global g_SegDup_coord
    global g_SegDup


    # headerOutput "SegDup_left" or "SegDup_right"
    # (SegDupAnnotation is executed for each breakpoint)
    set g_SegDup(Empty) ""

    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/BreakpointsAnnotations"
    set segdupFileFormatted [glob -nocomplain "$extannDir/SegDup/$g_AnnotSV(genomeBuild)/*_SegDup.sorted.bed"]

    if {![info exists g_SegDup(DONE)]} {
	
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
		puts "-- SegDupAnnotation, sort --"
		puts "sort -k1,1 -k2,2n $tmpBreakpointsFile.tmp > $tmpBreakpointsFile"
		puts "$Message"
		puts "Exit with error"
		exit 2
	    }
	    file delete -force $sortTmpFile 
	    file delete -force $tmpBreakpointsFile.tmp
	    
	}

	# Intersect
	regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".intersect.SegDup.bed" tmpSegDupFile
	set tmpSegDupFile "$g_AnnotSV(outputDir)/[file tail $tmpSegDupFile]"
	file delete -force $tmpSegDupFile
	if {[catch {exec $g_AnnotSV(bedtools) intersect -sorted -a $tmpBreakpointsFile -b $segdupFileFormatted -wa -wb > $tmpSegDupFile} Message]} {
	    if {[catch {exec $g_AnnotSV(bedtools) intersect -a $tmpBreakpointsFile -b $segdupFileFormatted -wa -wb > $tmpSegDupFile} Message]} {
		puts "-- SegDupAnnotation, intersect --"
		puts "$g_AnnotSV(bedtools) intersect -sorted -a $tmpBreakpointsFile -b $segdupFileFormatted -wa -wb > $tmpSegDupFile"
		puts "$Message"
		puts "Exit with error"
		exit 2
	    }
	}

	# Loading g_SegDup for each SV breakpoint
	set f [open $tmpSegDupFile]
	while {![eof $f]} {
	    set L [gets $f]
	    if {$L eq ""} {continue}
	    set Ls [split $L "\t"]
	    
	    set chrom [lindex $Ls 0]
	    set breakpoint [lindex $Ls 3]
	    #set segdupStart [lindex $Ls 5]
	    #set segdupEnd [lindex $Ls 6]
	    #set segdupCoord [lindex $Ls 7]
	    
	    lappend g_SegDup_coord($chrom,$breakpoint) "[lindex $Ls 7]"
	}
	
	set g_SegDup(DONE) 1	
	file delete -force $tmpSegDupFile
    }
    
    if {[info exist g_SegDup_coord($BreakpointChrom,$BreakpointPos)]} {
	return "[join $g_SegDup_coord($BreakpointChrom,$BreakpointPos) ";"]"
    } else {
	return $g_SegDup(Empty)
    }
}
