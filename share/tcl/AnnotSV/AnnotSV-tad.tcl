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


## - Check if the following TAD files has been downloaded:
#    - ENC*.bed
#  
## - Check and create if necessary the 'date'_boundariesTAD.sorted.bed file
proc checkTADfiles {} {

    global g_AnnotSV


    ## Check if TAD files has been downloaded then formatted
    #######################################################
    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)"
    set TADfilesDownloaded [glob -nocomplain "$extannDir/FtIncludedInSV/TAD/$g_AnnotSV(genomeBuild)/ENC*.bed"]
    set boundariesTADfileFormatted [glob -nocomplain "$extannDir/FtIncludedInSV/TAD/$g_AnnotSV(genomeBuild)/*_boundariesTAD.sorted.bed"] 

    if {$TADfilesDownloaded eq "" && $boundariesTADfileFormatted eq ""} {
	# No TAD annotation (with extAnn procedure)
	set g_AnnotSV(tadAnn) 0
	return
    } else {
	set g_AnnotSV(tadAnn) 1
	# Check if the user asked for these annotations in the configfile
	set test 0
	foreach col "TADcoordinates ENCODEexperiments" {
	    if {[lsearch -exact "$g_AnnotSV(outputColHeader)" $col] ne -1} {set test 1;break}
	}
	if {$test eq 0} {set g_AnnotSV(tadAnn) 0; return}
    }

    if {[llength $boundariesTADfileFormatted]>1} {
	puts "Several TAD files exist:"
	puts "$boundariesTADfileFormatted"
	puts "Keep only one: [lindex $boundariesTADfileFormatted end]\n"
	foreach tad [lrange $boundariesTADfileFormatted 0 end-1] {
	    file rename -force $tad $tad.notused
	}
	return
    } 

    if {$boundariesTADfileFormatted eq ""} {
	# The downloaded file exist but not the formatted.
	## Create:
	##   - 'date'_boundariesTAD.sorted.bed   ; # Header: chr boundaryStart boundaryEnd ENCODEexperiments correspondingTADs
	set boundariesTADfileFormatted "$extannDir/FtIncludedInSV/TAD/$g_AnnotSV(genomeBuild)/[clock format [clock seconds] -format "%Y%m%d"]_boundariesTAD.sorted.bed"

	puts "...TAD configuration"

	puts "\t...creation of $boundariesTADfileFormatted"
	puts "\t   (done only once during the first boundaries TAD annotation)\n"

	foreach TADfile [glob -nocomplain "$extannDir/FtIncludedInSV/TAD/$g_AnnotSV(genomeBuild)/ENC*.bed"] {
	    regsub -nocase ".bed$" $TADfile "" ENCODEexperiment
	    regsub ".*/" $ENCODEexperiment "" ENCODEexperiment
	    foreach L [LinesFromFile $TADfile] {
		regsub "chr" [lindex $L 0] "" chrom
		set TADstart [lindex $L 1]
		set TADend [lindex $L 2]
		if {[regexp "___boundary" $L]} {
		    set boundarieStart $TADstart
		    set boundariesEnd [expr $TADstart+40000]
		} else {
		    set boundarieStart [expr $TADstart-40000]
		    set boundariesEnd $TADstart
		    set TADstart $boundarieStart
		    set TADend [expr $TADend+40000]
		}
		lappend experiment($chrom\t$boundarieStart\t$boundariesEnd) "$ENCODEexperiment"
		lappend TAD($chrom\t$boundarieStart\t$boundariesEnd) "$chrom:$TADstart-$TADend"
		lappend L_allLines($chrom) "$boundarieStart\t$boundariesEnd"
	    }
	}
	
	foreach chrom {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M MT} {
	    if {![info exists L_allLines($chrom)]} {continue}
	    foreach bounds [lsort -command AscendingSortOnElement0 [lsort -command AscendingSortOnElement1 [lsort -unique $L_allLines($chrom)]]] {
		WriteTextInFile "$chrom\t$bounds\t[join [lsort -unique $experiment($chrom\t$bounds)] ","]\t[join [lsort -unique $TAD($chrom\t$bounds)] ","]" $boundariesTADfileFormatted.tmp
	    }
	}
	# Sorting of the bedfile:
	# Intersection with very large files can cause trouble with excessive memory usage.
	# A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm. 
	if {[catch {eval exec sort -k1,1 -k2,2n $boundariesTADfileFormatted.tmp > $boundariesTADfileFormatted} Message]} {
	    puts "-- checkTADfiles --"
	    puts "sort -k1,1 -k2,2n $boundariesTADfileFormatted.tmp > $boundariesTADfileFormatted"
	    puts "$Message"
	    puts "Exit with error"
	    exit 2
	}
	file delete -force $boundariesTADfileFormatted.tmp

	foreach TADfile $TADfilesDownloaded {
	    file delete -force $TADfile
	}

    }
}


proc TADannotation {SVchrom SVstart SVend L_i} {

    global g_AnnotSV
    global tadText


    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)"
    set boundariesTADfileFormatted [glob -nocomplain "$extannDir/FtIncludedInSV/TAD/$g_AnnotSV(genomeBuild)/*_boundariesTAD.sorted.bed"] 
    
    if {![info exists tadText(DONE)]} {
	
	# headerOutput "TADcoordinates\tENCODEexperiments"
	set L_tadText(Empty) "{} {}"
	foreach i $L_i {
	    lappend tadText(Empty) "[lindex $L_tadText(Empty) $i]"
	}
	set tadText(Empty) [join $tadText(Empty) "\t"]

	# Intersect
	regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".intersect.tad" tmpFile
	set tmpFile "$g_AnnotSV(outputDir)/[file tail $tmpFile]"
	file delete -force $tmpFile	
	if {[catch {exec $g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(bedFile) -b $boundariesTADfileFormatted -wa -wb > $tmpFile} Message]} {
	    puts "-- TADannotation --"
	    puts "$g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(bedFile) -b $boundariesTADfileFormatted -wa -wb > $tmpFile"
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

	    set TAD_coordinates  [lindex $Ls end]
	    set TAD_ENCODEexperiments [lindex $Ls end-1]
	    set boundary_end [lindex $Ls end-2]
	    set boundary_start [lindex $Ls end-3]

	    # Select:
	    # - TAD boundaries that share > XX % length with the SV to annotate
	    # - TAD boundaries overlapping a SV that is an insertion
	    set boundary_length [expr {$boundary_end-$boundary_start}]
	    set SVtoAnn_length [expr {$SVtoAnn_end - $SVtoAnn_start}]
	    if {$SVtoAnn_length<=1} {
		# SV to annotate is an insertion
		set SVtoAnn_length 1
	    } else {
		# SV to annotate is not an insertion
		if {$SVtoAnn_start < $boundary_start} {
		    set overlap_start $boundary_start
		} else {
		    set overlap_start $SVtoAnn_start
		}
		if {$SVtoAnn_end < $boundary_end} {
		    set overlap_end $SVtoAnn_end
		} else {
		    set overlap_end $boundary_end
		}
		set overlap_length [expr {$overlap_end - $overlap_start}]
		
		# Keeping only TAD with > 70% (default) overlap with the SV
		if {[expr {$overlap_length*100.0/$boundary_length}] < $g_AnnotSV(overlap)} {continue}
		# No reciprocal overlap used for TAD (it would not make sense)
	    }
	    
	    lappend L_experiments($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) {*}[split $TAD_ENCODEexperiments ","]
	    lappend L_TADcoordinates($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) {*}[split $TAD_coordinates ","]

	}
	

	# Loading TAD final annotation for each SV
	foreach SVtoAnn [array names L_experiments] {
	    # Very large SV (e.g. 30Mb) can overlap many TAD (e.g. more than 2600). 
	    # e.g. 1:2806107-107058351
	    # So the “full” line can appear truncated in the output if visualized in a spreadsheet.  
	    # In order to avoid such embarrassing glitch, AnnotSV gives access to the list of the 20 first TAD coordinates and their associated ENCODE experiments.
	    set L_tadText($SVtoAnn) ""
	    lappend L_tadText($SVtoAnn) [join [lrange [lsort -unique $L_TADcoordinates($SVtoAnn)] 0 20] ";"] 
	    lappend L_tadText($SVtoAnn) [join [lrange [lsort -unique $L_experiments($SVtoAnn)] 0 20] ":"]

	    # Keep only the user requested columns (defined in the configfile)
	    set tadText($SVtoAnn) ""
	    foreach i $L_i {
		lappend tadText($SVtoAnn) "[lindex $L_tadText($SVtoAnn) $i]"
	    }
	    set tadText($SVtoAnn) [join $tadText($SVtoAnn) "\t"]
	}
	catch {unset L_tadText}
	set tadText(DONE) 1	
	file delete -force $tmpFile
    }
    
    if {[info exist tadText($SVchrom,$SVstart,$SVend)]} {
	return $tadText($SVchrom,$SVstart,$SVend)
    } else {
	return $tadText(Empty)
    }
}
