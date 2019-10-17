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



## - Check if the following chromFa files has been downloaded:
#    - *chromFa.tar.gz
#
## - Check and create if necessary the following file:
#    - 'date'_genomeBuild_chromFa.fasta
#    - 'date'_genomeBuild_chromFa.chromSizes
proc checkFASTAfiles {} {

    global g_AnnotSV

    ## Check if the FASTA file has been downloaded then formatted
    ############################################################
    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/BreakpointsAnnotations"
    set FASTAfileDownloaded [glob -nocomplain "$extannDir/GCcontent/$g_AnnotSV(genomeBuild)/*chromFa.tar.gz"]
    set FASTAfileFormatted [glob -nocomplain "$extannDir/GCcontent/$g_AnnotSV(genomeBuild)/$g_AnnotSV(genomeBuild)_chromFa.fasta"]

    if {$FASTAfileDownloaded eq "" && $FASTAfileFormatted eq ""} {
	# No GCcontent annotations requested by the user
	set g_AnnotSV(gcContentAnn) 0
	return
    } else {
 	set g_AnnotSV(gcContentAnn) 1
	# Check if the user asked for these annotations in the configfile
	set test 0
	foreach col "GCcontent_left GCcontent_right" {
	    if {[lsearch -exact "$g_AnnotSV(outputColHeader)" $col] ne -1} {set test 1;break}
	}
	if {$test eq 0} {set g_AnnotSV(gcContentAnn) 0; return}
    }

    if {$FASTAfileFormatted eq ""} {
	# The downloaded file exist but not the formatted:
	updateFASTAfiles
    }
}

proc updateFASTAfiles {} {

    global g_AnnotSV

    if {[catch {package require tar} Message]} {
	puts "Tcl package \"tar\" is required for GC content annotations."
	puts "$Message"
	puts "-> no GC content annotations.\n"
	set g_AnnotSV(gcContentAnn) 0
	return
    }

    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/BreakpointsAnnotations"
    set FASTAfileDownloaded [glob -nocomplain "$extannDir/GCcontent/$g_AnnotSV(genomeBuild)/*chromFa.tar.gz"]
    set FASTAfileFormatted "$extannDir/GCcontent/$g_AnnotSV(genomeBuild)/$g_AnnotSV(genomeBuild)_chromFa.fasta"
    set ChromSizesFile "$extannDir/GCcontent/$g_AnnotSV(genomeBuild)/$g_AnnotSV(genomeBuild)_chromFa.chromSizes"

    puts "...GC content configuration ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    puts "\t...creation of $FASTAfileFormatted"
    puts "\t   (done only once during the first GC content annotation)\n"

    # Extracting files from the .tar.gz file
    set chan [open "|gzip -cd $FASTAfileDownloaded"]
    ::tar::untar $chan -chan -dir "$extannDir/GCcontent/$g_AnnotSV(genomeBuild)"

    # Merging in a unique file ($FASTAfileFormatted)
    foreach chrom {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M MT} {
	if {$g_AnnotSV(genomeBuild) eq "GRCh37"} {
	    set PathGB "$g_AnnotSV(genomeBuild)"
	} elseif {$g_AnnotSV(genomeBuild) eq "GRCh38"} {
	    set PathGB "$g_AnnotSV(genomeBuild)/chroms"
	} elseif {$g_AnnotSV(genomeBuild) eq "mm9"} {
	    set PathGB "$g_AnnotSV(genomeBuild)"
	} elseif {$g_AnnotSV(genomeBuild) eq "mm10"} {
	    set PathGB "$g_AnnotSV(genomeBuild)"
	}
	if {![file exists $extannDir/GCcontent/$PathGB/chr$chrom.fa]} {continue}
	WriteTextInFile [exec sed "s/chr//" $extannDir/GCcontent/$PathGB/chr$chrom.fa] $FASTAfileFormatted
	if {![regexp "Y|M" $chrom]} {
	    set wcValue [exec wc $extannDir/GCcontent/$PathGB/chr$chrom.fa]
	    set nbChars [lindex $wcValue 2]
	    set nLignes [lindex $wcValue 0]
	    puts "wcValue $wcValue; nbChars $nbChars; nLignes $nLignes - length [string length ">chr$chrom"]"
	    WriteTextInFile "$chrom\t[expr {$nbChars-$nLignes-[string length ">chr$chrom"]}]" $ChromSizesFile ; # Removing of the header and "\n" characters
	    # WARNING for chrY:
	    # >chrY dna:chromosome chromosome:GRCh37:Y:2649521:59034049:1
	    # The chr Y doesn't begin at the 1th position => The length of the sequence doesn't correspond to the length of the chrom Y
	}
    }

    foreach tmpFile [glob -nocomplain "$extannDir/GCcontent/$PathGB/chr*.fa"] {
	file delete -force $tmpFile
    }
    file delete -force $extannDir/GCcontent/GRCh38/chroms

    # index file generation
    ReplaceTextInFile "1\t150000\t150200" $extannDir/GCcontent/test.bed
    catch {exec bedtools nuc $FASTAfileFormatted -bed $extannDir/GCcontent/test.bed} Message
    file delete -force $extannDir/GCcontent/test.bed

    file delete -force $FASTAfileDownloaded

}



proc GCcontentAnnotation {BreakpointChrom BreakpointPos} {

    global g_AnnotSV
    global g_GCcontent


    # headerOutput "GCcontent_left GCcontent_right"
    # GCcontentAnnotation is executed for each breakpoint:
    # - the left breakpoint -> GCcontent_left
    # - the right one -> GCcontent_right
    set g_GCcontent(Empty) "" ; # Ok, no "\t" because we return only one value for each procedure call

    # Checking with the size of the chrom is needed
    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/BreakpointsAnnotations"
    set ChromSizesFile "$extannDir/GCcontent/$g_AnnotSV(genomeBuild)/$g_AnnotSV(genomeBuild)_chromFa.chromSizes"
    foreach L [LinesFromFile $ChromSizesFile] {
	set sizeOf([lindex $L 0]) [lindex $L 1]
    }

    set FASTAfileFormatted [glob -nocomplain "$extannDir/GCcontent/$g_AnnotSV(genomeBuild)/$g_AnnotSV(genomeBuild)_chromFa.fasta"]


    if {![info exists g_GCcontent(DONE)]} {

	# Creation of a bedfile with all breakpoints +/- 100 bp:
	# FORMAT:  chrom   'breakpoint-100'   'breakpoint+100'   'breakpoint'
	regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".breakpoints.bed" tmpBreakpointsFile
	set tmpBreakpointsFile "$g_AnnotSV(outputDir)/[file tail $tmpBreakpointsFile]"
	file delete -force $tmpBreakpointsFile.tmp
	file delete -force $tmpBreakpointsFile
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
	    if {[info exists sizeOf($chrom)] && $end_breakpoint1 > $sizeOf($chrom)} {
		# The start of the SV is near to the end of the chromosome
		if {$start_breakpoint1 < $sizeOf($chrom)} {
		    set end_breakpoint1 $sizeOf($chrom)
		}
	    }
	    set start_breakpoint2 [expr $end-100]  ; if {$start_breakpoint2 < 0} {set start_breakpoint2 0}
	    set end_breakpoint2 [expr $end+100]    ; if {$end_breakpoint2 < 0} {set end_breakpoint2 0}
	    if {[info exists sizeOf($chrom)] && $end_breakpoint2 > $sizeOf($chrom)} {
		# The end of the SV corresponds (or is near) to the end of the chromosome
		if {$start_breakpoint2 < $sizeOf($chrom)} {
		    set end_breakpoint2 $sizeOf($chrom)
		}
	    }

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

	# Sorting of the bedfile
	if {[catch {exec $g_AnnotSV(bedtools) sort -i $tmpBreakpointsFile.tmp > $tmpBreakpointsFile} Message]} {
	    puts "-- GCcontentAnnotation, sort --"
	    puts "$g_AnnotSV(bedtools) sort -i $tmpBreakpointsFile.tmp > $tmpBreakpointsFile"
	    puts "$Message"
	}
	file delete -force $tmpBreakpointsFile.tmp

	# Intersect
	regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".GCcontent.txt" tmpNucFile
	set tmpNucFile "$g_AnnotSV(outputDir)/[file tail $tmpNucFile]"
	file delete -force $tmpNucFile
	if {[catch {exec $g_AnnotSV(bedtools) nuc -fi $FASTAfileFormatted -bed $tmpBreakpointsFile > $tmpNucFile} Message]} {
	    if {$Message ne "Warning: the index file is older than the FASTA file."} {
		puts "-- GCcontentAnnotation, nuc --"
		puts "$g_AnnotSV(bedtools) nuc -fi $FASTAfileFormatted -bed $tmpBreakpointsFile > $tmpNucFile"
		puts "$Message"
	    }
	}
	# file delete -force $tmpBreakpointsFile; "tmpBreakpointsFile" is used later for the repeat annotation!

	# Loading g_GCcontent for each SV breakpoint
	set f [open $tmpNucFile]
	while {![eof $f]} {
	    set L [gets $f]
	    if {$L eq ""} {continue}
	    set Ls [split $L "\t"]
	    if {[regexp "_pct_gc" $L]} {
		set i_gc [lsearch -regexp $Ls "_pct_gc"]
		set i_n [lsearch -regexp $Ls "_num_N"]
		set i_oth [lsearch -regexp $Ls "_num_oth"]
		continue
	    }
	    set chrom [lindex $Ls 0]
	    set breakpoint [lindex $Ls 3]
	    set gcPercent [format "%.3f" [lindex $Ls $i_gc]]
	    set numN [lindex $Ls $i_n]
	    set numOth [lindex $Ls $i_oth]
	    # Change metrics from "." to ","
	    if {[set g_AnnotSV(metrics)] eq "fr"} {
		regsub -all {\.} $gcPercent "," gcPercent
	    }
	    set g_GCcontent($chrom,$breakpoint) $gcPercent
	    if {$numN ne "0"} {
		if {$numOth ne "0"} {append g_GCcontent($chrom,$breakpoint) " ($numN N, $numOth Oth)"} else {append g_GCcontent($chrom,$breakpoint) " ($numN N)"}
	    } elseif {$numOth ne "0"} {append g_GCcontent($chrom,$breakpoint) " ($numOth Oth)"}
	}

	set g_GCcontent(DONE) 1
	file delete -force $tmpNucFile
    }

    if {[info exist g_GCcontent($BreakpointChrom,$BreakpointPos)]} {
	return $g_GCcontent($BreakpointChrom,$BreakpointPos)
    } else {
	return $g_GCcontent(Empty)
    }
}
