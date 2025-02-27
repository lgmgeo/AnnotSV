############################################################################################################
# AnnotSV 3.4.5                                                                                            #
#                                                                                                          #
# AnnotSV: An integrated tool for Structural Variations annotation and ranking                             #
#                                                                                                          #
# Copyright (C) 2017-present Veronique Geoffroy (veronique.geoffroy@inserm.fr)                             #
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


## - Check if the following file exist:
#    - 'date'_Repeat.bed
proc checkRepeatFile {} {
    
    global g_AnnotSV
    
    ## Check if the repeat file has been downloaded then formatted
    #############################################################
    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/BreakpointsAnnotations"
    set repeatFileDownloaded [glob -nocomplain "$extannDir/Repeat/$g_AnnotSV(genomeBuild)/Repeat.bed"]
    set repeatFileFormatted [glob -nocomplain "$extannDir/Repeat/$g_AnnotSV(genomeBuild)/*_Repeat.sorted.bed"]
    
    if {$repeatFileDownloaded eq "" && $repeatFileFormatted eq ""} {
        # No Repeat annotation
        set g_AnnotSV(repeatAnn) 0
        return
    } else {
        # Repeat annotation
        set g_AnnotSV(repeatAnn) 1
        # Check if the user asked for these annotations in the configfile
        set test 0
        foreach col "Repeat_coord_left Repeat_type_left Repeat_coord_right Repeat_type_right" {
            if {[lsearch -exact "$g_AnnotSV(outputColHeader)" $col] ne -1} {set test 1;break}
        }
        if {$test eq 0} {set g_AnnotSV(repeatAnn) 0; return}
    }
    
    if {[llength $repeatFileFormatted]>1} {
        puts "Several repeat files exist:"
        puts "$repeatFileFormatted"
        puts "Keep only one: [lindex $repeatFileFormatted end]\n"
        foreach repeat [lrange $repeatFileFormatted 0 end-1] {
            file rename -force $repeat $repeat.notused
        }
        return
    }
    
    if {$repeatFileFormatted eq ""} {
        # The downloaded file exist but not the formatted.
        set repeatFileFormatted "$extannDir/Repeat/$g_AnnotSV(genomeBuild)/[clock format [clock seconds] -format "%Y%m%d"]_Repeat.sorted.bed"
        file delete -force $repeatFileFormatted.tmp
        puts "\t...Repeats configuration ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
        puts "\t\t...creation of $repeatFileFormatted"
        puts "\t\t   (done only once during the first Repeats annotation)"
        
        set f [open $repeatFileDownloaded]
        set i 0
        set L_Text {}
        while {![eof $f]} {
            set L [gets $f]
            if {$L eq ""} {continue}
            regsub "^chr" $L "" L
            set chrom [lindex $L 0]
            if {[lsearch -exact {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 X Y M MT} $chrom] eq -1} {continue}
            lappend L_Text "$L"
            if {$i>500000} {
                WriteTextInFile [join $L_Text "\n"] $repeatFileFormatted.tmp
                set L_Text {}
                set i 0
            }
        }
        WriteTextInFile [join $L_Text "\n"] $repeatFileFormatted.tmp
        
        # Sorting of the bedfile:
        # Intersection with very large files can cause trouble with excessive memory usage.
        # A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm.
        set sortTmpFile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_sort.tmp.bash"
        ReplaceTextInFile "#!/bin/bash" $sortTmpFile
        WriteTextInFile "# The locale specified by the environment can affects the traditional sort order. We need to use native byte values." $sortTmpFile
        WriteTextInFile "export LC_ALL=C" $sortTmpFile
        WriteTextInFile "sort -k1,1 -k2,2n $repeatFileFormatted.tmp > $repeatFileFormatted" $sortTmpFile
        file attributes $sortTmpFile -permissions 0755
        if {[catch {eval exec bash $sortTmpFile} Message]} {
            puts "-- checkRepeatFile --"
            puts "sort -k1,1 -k2,2n $repeatFileFormatted.tmp > $repeatFileFormatted"
            puts "$Message"
            puts "Exit with error"
            exit 2
        }
        file delete -force $sortTmpFile
        file delete -force $repeatFileFormatted.tmp
        file delete -force $repeatFileDownloaded
        
        close $f
        
        return
    }
}



proc RepeatAnnotation {BreakpointChrom BreakpointPos} {
    
    global g_AnnotSV
    global g_Repeat_type
    global g_Repeat_coord
    global g_Repeat
    
    
    # headerOutput "Repeat_coord_left Repeat_type_left" or "Repeat_coord_right and Repeat_type_right"
    # (RepeatAnnotation is executed for each breakpoint)
    set g_Repeat(Empty) "\t"
    
    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/BreakpointsAnnotations"
    set repeatFileFormatted [glob -nocomplain "$extannDir/Repeat/$g_AnnotSV(genomeBuild)/*_Repeat.sorted.bed"]
    
    if {![info exists g_Repeat(DONE)]} {
        
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
                puts "-- RepeatAnnotation, sort --"
                puts "sort -k1,1 -k2,2n $tmpBreakpointsFile.tmp > $tmpBreakpointsFile"
                puts "$Message"
                puts "Exit with error"
                exit 2
            }
            file delete -force $sortTmpFile
            file delete -force $tmpBreakpointsFile.tmp
            
        }
        
        # Intersect
        regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".intersect.Repeat.bed" tmpRepeatFile
        set tmpRepeatFile "$g_AnnotSV(outputDir)/[file tail $tmpRepeatFile]"
        file delete -force $tmpRepeatFile
        if {[catch {exec $g_AnnotSV(bedtools) intersect -sorted -a $tmpBreakpointsFile -b $repeatFileFormatted -wa -wb > $tmpRepeatFile} Message]} {
            if {[catch {exec $g_AnnotSV(bedtools) intersect -a $tmpBreakpointsFile -b $repeatFileFormatted -wa -wb > $tmpRepeatFile} Message]} {
                puts "-- RepeatAnnotation, intersect --"
                puts "$g_AnnotSV(bedtools) intersect -sorted -a $tmpBreakpointsFile -b $repeatFileFormatted -wa -wb > $tmpRepeatFile"
                puts "$Message"
                puts "Exit with error"
                exit 2
            }
        }
        
        # Loading g_Repeat for each SV breakpoint
        set f [open $tmpRepeatFile]
        while {![eof $f]} {
            set L [gets $f]
            if {$L eq ""} {continue}
            set Ls [split $L "\t"]
            
            set chrom [lindex $Ls 0]
            set breakpoint [lindex $Ls 3]
            set repeatStart [lindex $Ls 5]
            set repeatEnd [lindex $Ls 6]
            set repeatType [lindex $Ls 7]
            
            lappend g_Repeat_coord($chrom,$breakpoint) "$chrom:$repeatStart-$repeatEnd"
            lappend g_Repeat_type($chrom,$breakpoint) "$repeatType"
        }
        
        set g_Repeat(DONE) 1
        file delete -force $tmpRepeatFile
        
    }
    
    if {[info exist g_Repeat_coord($BreakpointChrom,$BreakpointPos)]} {
        return "[join $g_Repeat_coord($BreakpointChrom,$BreakpointPos) ";"]\t[join $g_Repeat_type($BreakpointChrom,$BreakpointPos) ";"]"
    } else {
        return $g_Repeat(Empty)
    }
}

