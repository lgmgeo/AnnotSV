############################################################################################################
# AnnotSV 3.4.6                                                                                            #
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



proc closestGenesAnnotation {BreakpointChrom BreakpointPos Side} {
    
    global g_AnnotSV
    global g_closestGenesLeft
    global g_closestGenesRight

    
 
    if {![info exists g_closestGenesLeft(DONE)]} {
        
        # Creation of "$tmpCLosestLeftFile" with:
        #   SV Left coordinates: "start-5000000bp" "start"
        #   => To find the closest left gene
		#
        # Creation of "$tmpCLosestRightFile" with:
        #   SV Right coordinates: "end" "end+5000000bp"
        #   => To find the closest right gene
        
        regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".closestLeft.bed" tmpCLosestLeftFile
        regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".closestRight.bed" tmpCLosestRightFile
        
        set tmpCLosestLeftFile  "$g_AnnotSV(outputDir)/[file tail $tmpCLosestLeftFile]"
        set tmpCLosestRightFile "$g_AnnotSV(outputDir)/[file tail $tmpCLosestRightFile]"
        
        file delete -force $tmpCLosestLeftFile.tmp
        file delete -force $tmpCLosestRightFile.tmp
        
        set i 0
        set L_LinesToWriteLeft {}
        set L_LinesToWriteRight {}
        set f [open $g_AnnotSV(bedFile)]
        while {![eof $f]} {
            set L [gets $f]
            if {$L eq ""} {continue}
            set Ls [split $L "\t"]
            set chrom [lindex $Ls 0]
            set start [lindex $Ls 1]
            set end [lindex $Ls 2]
            set start_left [expr $start-5000000]; if {$start_left < 0} {set start_left 0}
            set end_right [expr $end+5000000]
            lappend L_LinesToWriteLeft "$chrom\t$start_left\t$start"
            lappend L_LinesToWriteRight "$chrom\t$end\t$end_right"
            incr i
            if {$i>100000} {
                WriteTextInFile [join [lsort -unique $L_LinesToWriteLeft] "\n"] $tmpCLosestLeftFile.tmp
                WriteTextInFile [join [lsort -unique $L_LinesToWriteRight] "\n"] $tmpCLosestRightFile.tmp
                set L_LinesToWriteLeft {}
                set L_LinesToWriteRight {}
                set i 0
            }
        }
        close $f
        WriteTextInFile [join [lsort -unique $L_LinesToWriteLeft] "\n"] $tmpCLosestLeftFile.tmp
        WriteTextInFile [join [lsort -unique $L_LinesToWriteRight] "\n"] $tmpCLosestRightFile.tmp

        unset L_LinesToWriteLeft
        unset L_LinesToWriteRight
        
        # Right sorting: $tmpCLosestLeftFile.tmp => ($bashTmpFileLeft) => $tmpCLosestLeftFile
		#
        # Intersection with very large files can cause trouble with excessive memory usage.
        # A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm.
        set bashTmpFileLeft "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_sort.left.tmp.bash"
        ReplaceTextInFile "#!/bin/bash" $bashTmpFileLeft
        WriteTextInFile "# The locale specified by the environment can affects the traditional sort order. We need to use native byte values." $bashTmpFileLeft
        WriteTextInFile "export LC_ALL=C" $bashTmpFileLeft
        WriteTextInFile "sort -k1,1 -k2,2n $tmpCLosestLeftFile.tmp > $tmpCLosestLeftFile" $bashTmpFileLeft
        file attributes $bashTmpFileLeft -permissions 0755
        if {[catch {eval exec bash $bashTmpFileLeft} Message]} {
            puts "-- closestGenesAnnotation, left sort --"
            puts "sort -k1,1 -k2,2n $tmpCLosestLeftFile.tmp > $tmpCLosestLeftFile"
            puts "$Message"
            puts "Exit with error"
            exit 2
        }
        file delete -force $bashTmpFileLeft
        file delete -force $tmpCLosestLeftFile.tmp
        
        # Right sorting: $tmpCLosestRightFile.tmp => ($bashTmpFileRight) => $tmpCLosestRightFile
        set bashTmpFileRight "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_sort.left.tmp.bash"
        ReplaceTextInFile "#!/bin/bash" $bashTmpFileRight
        WriteTextInFile "# The locale specified by the environment can affects the traditional sort order. We need to use native byte values." $bashTmpFileRight
        WriteTextInFile "export LC_ALL=C" $bashTmpFileRight
        WriteTextInFile "sort -k1,1 -k2,2n $tmpCLosestRightFile.tmp > $tmpCLosestRightFile" $bashTmpFileRight
        file attributes $bashTmpFileRight -permissions 0755
        if {[catch {eval exec bash $bashTmpFileRight} Message]} {
            puts "-- closestGenesAnnotation, right sort --"
            puts "sort -k1,1 -k2,2n $tmpCLosestRightFile.tmp > $tmpCLosestRightFile"
            puts "$Message"
            puts "Exit with error"
            exit 2
        }
        file delete -force $bashTmpFileRight
        file delete -force $tmpCLosestRightFile.tmp

    
		# RefSeq or ENSEMBL gene file
		set genesDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Genes/$g_AnnotSV(genomeBuild)"
	    set GenesFileFormatted $genesDir/genes.$g_AnnotSV(tx).sorted.bed
 
        # Intersect Left
        regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".intersect.closestGenesLeft.bed" tmpIntersectGeneFileLeft
        set tmpIntersectGeneFileLeft "$g_AnnotSV(outputDir)/[file tail $tmpIntersectGeneFileLeft]"
        file delete -force $tmpIntersectGeneFileLeft
        if {[catch {exec $g_AnnotSV(bedtools) intersect -sorted -a $tmpCLosestLeftFile -b $GenesFileFormatted -wa -wb > $tmpIntersectGeneFileLeft} Message]} {
            if {[catch {exec $g_AnnotSV(bedtools) intersect -a $tmpCLosestLeftFile -b $GenesFileFormatted -wa -wb > $tmpIntersectGeneFileLeft} Message]} {
                puts "-- closestGenesAnnotation, left intersect --"
                puts "$g_AnnotSV(bedtools) intersect -sorted -a $tmpCLosestLeftFile -b $GenesFileFormatted -wa -wb > $tmpIntersectGeneFileLeft"
                puts "$Message"
                puts "Exit with error"
                exit 2
            }
        }
        # Intersect Right
        regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".intersect.closestGenesRight.bed" tmpIntersectGeneFileRight
        set tmpIntersectGeneFileRight "$g_AnnotSV(outputDir)/[file tail $tmpIntersectGeneFileRight]"
        file delete -force $tmpIntersectGeneFileRight
        if {[catch {exec $g_AnnotSV(bedtools) intersect -sorted -a $tmpCLosestRightFile -b $GenesFileFormatted -wa -wb > $tmpIntersectGeneFileRight} Message]} {
            if {[catch {exec $g_AnnotSV(bedtools) intersect -a $tmpCLosestRightFile -b $GenesFileFormatted -wa -wb > $tmpIntersectGeneFileRight} Message]} {
                puts "-- closestGenesAnnotation, right intersect --"
                puts "$g_AnnotSV(bedtools) intersect -sorted -a $tmpCLosestRightFile -b $GenesFileFormatted -wa -wb > $tmpIntersectGeneFileRight"
                puts "$Message"
                puts "Exit with error"
                exit 2
            }
        }
 
        # Loading g_closestGenesLeft($chrom,$pos)
        set f [open $tmpIntersectGeneFileLeft]
        while {![eof $f]} {
            set L [gets $f]
            if {$L eq ""} {continue}
            set Ls [split $L "\t"]
            
            # SV Left coordinates: "chrom" "start-5000000bp" "start"
            set chrom [lindex $Ls 0]
            set pos [lindex $Ls 2] 
            
			if {![info exists g_closestGenesLeft($chrom,$pos)]} {
				lappend g_closestGenesLeft($chrom,$pos) [lindex $Ls 7]
			} else {
				if {[lsearch -exact $g_closestGenesLeft($chrom,$pos) [lindex $Ls 7]] eq -1} {
					lappend g_closestGenesLeft($chrom,$pos) [lindex $Ls 7]
				}
			}
        }
        file delete -force $tmpIntersectGeneFileLeft

        # Loading g_closestGenesRight($chrom,$pos)
        set f [open $tmpIntersectGeneFileRight]
        while {![eof $f]} {
            set L [gets $f]
            if {$L eq ""} {continue}
            set Ls [split $L "\t"]

	        # SV Right coordinates: "chrom" "end" "end+5000000bp"
            set chrom [lindex $Ls 0]
            set pos [lindex $Ls 1] 

            if {![info exists g_closestGenesRight($chrom,$pos)]} {
                lappend g_closestGenesRight($chrom,$pos) [lindex $Ls 7]
            } else {
 				if {[lsearch -exact $g_closestGenesRight($chrom,$pos) [lindex $Ls 7]] eq -1} {
					lappend g_closestGenesRight($chrom,$pos) [lindex $Ls 7]
				}
			}
        }
        file delete -force $tmpIntersectGeneFileRight


		# Memorisation done
        set g_closestGenesLeft(DONE) 1

		# Clean
        file delete -force $tmpCLosestLeftFile
        file delete -force $tmpCLosestRightFile

    }
   
	if {$Side eq "left"} { 
	    if {[info exist g_closestGenesLeft($BreakpointChrom,$BreakpointPos)]} {
	        return "$g_closestGenesLeft($BreakpointChrom,$BreakpointPos)"
	    } else {
	        return ""
	    }
	} elseif {$Side eq "right"} {
        if {[info exist g_closestGenesRight($BreakpointChrom,$BreakpointPos)]} {
            return "$g_closestGenesRight($BreakpointChrom,$BreakpointPos)"
        } else {
            return ""
        }
	} else {
		return ""
	}
}

