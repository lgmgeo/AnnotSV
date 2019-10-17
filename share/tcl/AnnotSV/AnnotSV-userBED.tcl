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


# Check the Users BED regions annotations files (from .../Annotations_$g_AnnotSV(organism)/Users/GRCh*/*ncludedIn*/*.bed)
# Create : .../Annotations_$g_AnnotSV(organism)/Users/GRCh*/*ncludedIn*/*.formatted.sorted
# (After the formatting step, the copy and/or linked users file(s) are deleted)
proc checkUsersBED {} {

    global g_AnnotSV
    global g_numberOfAnnotationCol

    set userBEDdir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Users/$g_AnnotSV(genomeBuild)"
    foreach userBEDfile [glob -nocomplain $userBEDdir/*ncludedIn*/*.bed] {
	# Create a formatted and sorted bedfile if it doesn't exist yet
	# and create a *.header.tsv file associated
	regsub -nocase "(.formatted.sorted)?.bed$" $userBEDfile ".bed" notFormattedFile
	regsub -nocase "(.formatted.sorted)?.bed$" $userBEDfile ".formatted.bed" formattedFile
	regsub -nocase "(.formatted.sorted)?.bed$" $userBEDfile ".formatted.sorted.bed" formattedSortedFile
	regsub -nocase ".bed$" $notFormattedFile ".header.tsv" userHeaderFile
	if {[file exists $notFormattedFile]} {
	    # Create the formatted bedfile and check the header existence
	    checkBed $notFormattedFile [file dirname $notFormattedFile]
	    # Create the formatted and sorted bedfile
	    #   Sorting of the bedfile:
	    #   Intersection with very large files can cause trouble with excessive memory usage.
	    #   A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm. 
	    if {[catch {eval exec sort -k1,1 -k2,2n $formattedFile > $formattedSortedFile} Message]} {
		puts "-- checkUsersBED --"
		puts "sort -k1,1 -k2,2n $formattedFile > $formattedSortedFile"
		puts "$Message"
		puts "Exit with error"
		exit 2
	    }
	    file delete -force $formattedFile
	}

	# Check that the bed and header files have the same number of columns (and more than 3)
	set L1header   [split [FirstLineFromFile $userHeaderFile] "\t"]
	set L_headerColName [lrange $L1header 3 end]
	set L1bed      [split [FirstLineFromFile $formattedSortedFile] "\t"]
	set nColHeader [llength $L1header]
	set nColBed    [llength $L1bed]
	if {$nColHeader ne $nColBed} {
	    puts "WARNING: [file tail $userHeaderFile] has $nColHeader fields, but $nColBed were expected. Header bed file not used.\n"
	    file rename -force $userHeaderFile $userHeaderFile.bad
	    WriteTextInFile "[join [lrepeat $nColBed ""] "\t"]" $userHeaderFile
	    set L_headerColName [expr {$nColBed-3}]
	} elseif {$nColHeader <=3} {
	    puts "WARNING: [file tail $userHeaderFile] has $nColHeader fields. More than 3 are expected. Header bed file not used.\n"
	    file rename -force $userHeaderFile $userHeaderFile.bad
	    WriteTextInFile "[join [lrepeat $nColBed ""] "\t"]" $userHeaderFile
	    set L_headerColName [expr {$nColBed-3}]
	}
	
	# Number of annotation columns in the user bedfile (without the 3 columns "chrom start end")
	set g_numberOfAnnotationCol($formattedSortedFile) [llength $L_headerColName]
    }
}


# Return the annotation of a SV with the userBED file
proc userBEDannotation {formattedSortedUserBEDfile SVchrom SVstart SVend} {

    global g_AnnotSV
    global g_numberOfAnnotationCol
    global userBEDtext

    if {![info exists userBEDtext($formattedSortedUserBEDfile,DONE)]} {
	set userBEDtext($formattedSortedUserBEDfile,Empty) ""
	
	regsub -nocase ".formatted.sorted.bed$" $formattedSortedUserBEDfile ".header.tsv" userHeaderFile
	set L_headerColName [lrange [split [FirstLineFromFile $userHeaderFile] "\t"] 3 end]
	set nColHeader [llength $L_headerColName]
	set g_AnnotSV($formattedSortedUserBEDfile,userBEDAnn) 1
		    
	# no intersection with the SV 
	append userBEDtext($formattedSortedUserBEDfile,Empty) "\t[join [lrepeat $g_numberOfAnnotationCol($formattedSortedUserBEDfile) ""] "\t"]"

	# Intersect
	# $g_AnnotSV(bedFile) = formatted and sorted bedfile!
	regsub -nocase ".formatted.sorted.bed$" $g_AnnotSV(bedFile) ".intersect.userBED" tmpFile
	set tmpFile "$g_AnnotSV(outputDir)/[file tail $tmpFile]"
	file delete -force $tmpFile	
	if {[catch {exec $g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(fullAndSplitBedFile) -b $formattedSortedUserBEDfile -wa -wb > $tmpFile} Message]} {
	    puts "-- userBEDannotation --"
	    puts "$g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(fullAndSplitBedFile) -b $formattedSortedUserBEDfile -wa -wb > $tmpFile"
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
	    
	    set userBED_end [lindex $Ls end-$g_numberOfAnnotationCol($formattedSortedUserBEDfile)]
	    set userBED_start [lindex $Ls end-[expr $g_numberOfAnnotationCol($formattedSortedUserBEDfile)+1]]
	    
	    # Select:
	    # - userBED that share > XX% length with the SV to annotate
	    # (doesn't select the INSERTION) 
	    set userBED_length [expr {$userBED_end-$userBED_start}]
	    set SVtoAnn_length [expr {$SVtoAnn_end - $SVtoAnn_start}]
	    # The user SV is an insertion or a breakpoint
	    if {$userBED_length<1} {
		set userBED_length 1
	    } 	
	    # The SV to annotate is an insertion or a breakpoint
	    if {$SVtoAnn_length<1} {
		set SVtoAnn_length 1
	    } 	

	    if {$SVtoAnn_start < $userBED_start} {
		set overlap_start $userBED_start
	    } else {
		set overlap_start $SVtoAnn_start
	    }
	    if {$SVtoAnn_end < $userBED_end} {
		set overlap_end $SVtoAnn_end
	    } else {
		set overlap_end $userBED_end
	    }
	    set overlap_length [expr {$overlap_end - $overlap_start}]
	    
	    if {[regexp "FtIncludedInSV" $formattedSortedUserBEDfile]} {
		## FtIncludedInSV:
		## <=> Keeping only user regions (userRegion=Ft=Feature) included in the SV to annotate
		## <=> Keeping only Ft with > 70% (default) overlap with the SV
		if {[expr {$overlap_length*100.0/$userBED_length}] < $g_AnnotSV(overlap)} {continue}
	    } elseif {[regexp "SVincludedInFt" $formattedSortedUserBEDfile]} {
		## SVincludedInFt:
		## <=> Keeping user regions respecting the overlaps (reciprocal or not reciprocal)
		if {[expr {$overlap_length*100.0/$SVtoAnn_length}] < $g_AnnotSV(overlap)} {continue}
		if {$g_AnnotSV(reciprocal) eq "yes"} {			
		    if {[expr {$overlap_length*100.0/$userBED_length}] < $g_AnnotSV(overlap)} {continue}
		}
	    }
	   
	    
	    # Each SV to annotate can be overlapped with several users regions => use of "$i" to merge the annotations column by column
	    set i "$g_numberOfAnnotationCol($formattedSortedUserBEDfile)"
	    while {$i > 0} {
		lappend L_userBEDAnn_${i}($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) [lindex $Ls end-[expr {$i-1}]]
		incr i -1
	    }
	}
	
	# Loading userBED final annotation for each SV
	set i "$g_numberOfAnnotationCol($formattedSortedUserBEDfile)"
	while {$i > 0} {
	    foreach SVtoAnn [array names L_userBEDAnn_$i] {
		append userBEDtext($formattedSortedUserBEDfile,$SVtoAnn) "\t[join [set L_userBEDAnn_${i}($SVtoAnn)] "/"]"
	    }
	    incr i -1
	}
	
	## Final reajustement
	set i "$g_numberOfAnnotationCol($formattedSortedUserBEDfile)"
	foreach SVtoAnn [array names L_userBEDAnn_$i] {
	    regsub "^\t" $userBEDtext($formattedSortedUserBEDfile,$SVtoAnn) "" userBEDtext($formattedSortedUserBEDfile,$SVtoAnn)
	}
	regsub "^\t" $userBEDtext($formattedSortedUserBEDfile,Empty) "" userBEDtext($formattedSortedUserBEDfile,Empty)


	## Annotation done
	set userBEDtext($formattedSortedUserBEDfile,DONE) 1	
	set g_AnnotSV(userBEDAnn) 1
	file delete -force $tmpFile		
    }

    if {[info exist userBEDtext($formattedSortedUserBEDfile,$SVchrom,$SVstart,$SVend)]} {
	return $userBEDtext($formattedSortedUserBEDfile,$SVchrom,$SVstart,$SVend)
    } else {
	return $userBEDtext($formattedSortedUserBEDfile,Empty)
    }
}
