############################################################################################################
# AnnotSV 3.0.8                                                                                            #
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


## - Check if the downloaded "ENCODEblacklist.bed" file or the formatted "*_ENCODEblacklist.sorted.bed" file exist
proc checkENCODEblacklistFile {} {

    global g_AnnotSV

    ## Check if the encodeblacklist file has been downloaded then formatted
    #############################################################
    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/BreakpointsAnnotations"
    set encodeblacklistFileDownloaded [glob -nocomplain "$extannDir/ENCODEblacklist/$g_AnnotSV(genomeBuild)/ENCODEblacklist.bed"]   
    set encodeblacklistFileFormatted [glob -nocomplain "$extannDir/ENCODEblacklist/$g_AnnotSV(genomeBuild)/*_ENCODEblacklist.sorted.bed"]   

    if {$encodeblacklistFileDownloaded eq "" && $encodeblacklistFileFormatted eq ""} {
	# No ENCODEblacklist annotation
	set g_AnnotSV(ENCODEblacklistAnn) 0
	return
    } else {
	# ENCODEblacklist annotation
	set g_AnnotSV(ENCODEblacklistAnn) 1
	# Check if the user asked for these annotations in the configfile
	set test 0
	foreach col "ENCODE_blacklist_left ENCODE_blacklist_characteristics_left ENCODE_blacklist_right ENCODE_blacklist_characteristics_right" {
	    if {[lsearch -exact "$g_AnnotSV(outputColHeader)" $col] ne -1} {set test 1;break}
	}
	if {$test eq 0} {set g_AnnotSV(ENCODEblacklistAnn) 0; return}
    }

    if {[llength $encodeblacklistFileFormatted]>1} {
	puts "Several ENCODE blacklist files exist:"
	puts "$encodeblacklistFileFormatted"
	puts "Keep only one: [lindex $encodeblacklistFileFormatted end]\n"
	foreach encodeblacklist [lrange $encodeblacklistFileFormatted 0 end-1] {
	    file rename -force $encodeblacklist $encodeblacklist.notused
	}
	return
    } 

    if {$encodeblacklistFileFormatted eq ""} {
	# The downloaded file exist but not the formatted.
	set encodeblacklistFileFormatted "$extannDir/ENCODEblacklist/$g_AnnotSV(genomeBuild)/[clock format [clock seconds] -format "%Y%m%d"]_ENCODEblacklist.sorted.bed"
	ReplaceTextInFile "#Chrom\tStart\tEnd\tENCODE_blacklist\tENCODE_blacklist_characteristics" $encodeblacklistFileFormatted.tmp

	set f [open $encodeblacklistFileDownloaded]
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
            set IDencodeblacklist "${chrom}:${start}-${end}"
	    set characteristics [lindex $Ls 3]
	    
	    lappend L_Text "$chrom\t$start\t$end\t$IDencodeblacklist\t$characteristics"
	    if {$i>500000} {
		WriteTextInFile [join $L_Text "\n"] $encodeblacklistFileFormatted.tmp
		set L_Text {}
		set i 0
	    }
            incr i
	}
	WriteTextInFile [join $L_Text "\n"] $encodeblacklistFileFormatted.tmp

	# Sorting of the bedfile:
	# Intersection with very large files can cause trouble with excessive memory usage.
	# A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm. 
	set sortTmpFile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_sort.tmp.bash"
	ReplaceTextInFile "#!/bin/bash" $sortTmpFile
	WriteTextInFile "# The locale specified by the environment can affects the traditional sort order. We need to use native byte values." $sortTmpFile
	WriteTextInFile "export LC_ALL=C" $sortTmpFile
	WriteTextInFile "sort -k1,1 -k2,2n $encodeblacklistFileFormatted.tmp > $encodeblacklistFileFormatted" $sortTmpFile
	file attributes $sortTmpFile -permissions 0755
	if {[catch {eval exec bash $sortTmpFile} Message]} {
	    puts "-- checkENCODEblacklistFile --"
	    puts "sort -k1,1 -k2,2n $encodeblacklistFileFormatted.tmp > $encodeblacklistFileFormatted"
	    puts "$Message"
	    puts "Exit with error"
	    exit 2
	}
	file delete -force $sortTmpFile 
	file delete -force $encodeblacklistFileFormatted.tmp
	file delete -force $encodeblacklistFileDownloaded

	close $f
	
	return
    }
}



proc ENCODEblacklistAnnotation {BreakpointChrom BreakpointPos} {

    global g_AnnotSV
    global g_ENCODEblacklist_coord
    global g_ENCODEblacklist_charact
    global g_ENCODEblacklist


    # headerOutput "ENCODE_blacklists_coord_left ENCODE_blacklists_type_left" or "ENCODE_blacklists_coord_right and ENCODE_blacklists_type_right"
    # (ENCODEblacklistAnnotation is executed for each breakpoint)
    set g_ENCODEblacklist(Empty) "\t"

    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/BreakpointsAnnotations"
    set encodeblacklistFileFormatted [glob -nocomplain "$extannDir/ENCODEblacklist/$g_AnnotSV(genomeBuild)/*_ENCODEblacklist.sorted.bed"]

    if {![info exists g_ENCODEblacklist(DONE)]} {
	
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
		puts "-- ENCODEblacklistAnnotation, sort --"
		puts "sort -k1,1 -k2,2n $tmpBreakpointsFile.tmp > $tmpBreakpointsFile"
		puts "$Message"
		puts "Exit with error"
		exit 2
	    }
	    file delete -force $sortTmpFile 
	    file delete -force $tmpBreakpointsFile.tmp
	    
	}

	# Intersect
	regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".intersect.ENCODEblacklist.bed" tmpENCODEblacklistFile
	set tmpENCODEblacklistFile "$g_AnnotSV(outputDir)/[file tail $tmpENCODEblacklistFile]"
	file delete -force $tmpENCODEblacklistFile
	if {[catch {exec $g_AnnotSV(bedtools) intersect -sorted -a $tmpBreakpointsFile -b $encodeblacklistFileFormatted -wa -wb > $tmpENCODEblacklistFile} Message]} {
	    if {[catch {exec $g_AnnotSV(bedtools) intersect -a $tmpBreakpointsFile -b $encodeblacklistFileFormatted -wa -wb > $tmpENCODEblacklistFile} Message]} {
		puts "-- ENCODEblacklistAnnotation, intersect --"
		puts "$g_AnnotSV(bedtools) intersect -sorted -a $tmpBreakpointsFile -b $encodeblacklistFileFormatted -wa -wb > $tmpENCODEblacklistFile"
		puts "$Message"
		puts "Exit with error"
		exit 2
	    }
	}

	# Loading g_ENCODEblacklist for each SV breakpoint
	set f [open $tmpENCODEblacklistFile]
	while {![eof $f]} {
	    set L [gets $f]
	    if {$L eq ""} {continue}
	    set Ls [split $L "\t"]
	    
	    set chrom [lindex $Ls 0]
	    set breakpoint [lindex $Ls 3]
	    #set encodeblacklistChrom [lindex $Ls 4]
	    #set encodeblacklistStart [lindex $Ls 5]
	    #set encodeblacklistEnd [lindex $Ls 6]
	    #set encodeblacklistCoord [lindex $Ls 7]
	    #set encodeblacklistCharact [lindex $Ls 8]
  
	    lappend g_ENCODEblacklist_coord($chrom,$breakpoint) [lindex $Ls 7]
	    lappend g_ENCODEblacklist_charact($chrom,$breakpoint) [lindex $Ls 8]
	}
	
	set g_ENCODEblacklist(DONE) 1	
	file delete -force $tmpENCODEblacklistFile
    }
    
    if {[info exist g_ENCODEblacklist_coord($BreakpointChrom,$BreakpointPos)]} {
	return "[join $g_ENCODEblacklist_coord($BreakpointChrom,$BreakpointPos) ";"]\t[join $g_ENCODEblacklist_charact($BreakpointChrom,$BreakpointPos) ";"]"
    } else {
	return $g_ENCODEblacklist(Empty)
    }
}
