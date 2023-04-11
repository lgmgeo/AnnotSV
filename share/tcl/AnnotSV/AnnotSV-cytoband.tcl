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



# Check the cytoband annotations files 
# Create : .../Annotations_$g_AnnotSV(organism)/AnyOverlap/CytoBand/$g_AnnotSV(genomeBuild)/cytoBand_$g_AnnotSV(genomeBuild).formatted.sorted.bed
# (After the formatting step, the "cytoBand_$g_AnnotSV(genomeBuild).bed" is deleted)
proc checkCytoband {} {

    global g_AnnotSV
    global g_numberOfAnnotationCol

    set cytobandDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/AnyOverlap/CytoBand/$g_AnnotSV(genomeBuild)"
    set cytobandBEDfile     "$cytobandDir/cytoBand_$g_AnnotSV(genomeBuild).bed"
    set formattedFile       "$cytobandDir/cytoBand_$g_AnnotSV(genomeBuild).formatted.bed" 
    set formattedSortedFile "$cytobandDir/cytoBand_$g_AnnotSV(genomeBuild).formatted.sorted.bed"
    set cytobandHeaderFile  "$cytobandDir/cytoBand_$g_AnnotSV(genomeBuild).header.tsv" 
    
    if {[file exists $cytobandBEDfile]} {
	
	# Create the formatted bedfile and check the header existence
	checkBed $cytobandBEDfile $cytobandDir

	# Create the cytoband header file
	WriteTextInFile "#chrom\tstart\tend\tCytoBand" $cytobandHeaderFile
	
	# Create the formatted and sorted bedfile
	#   Sorting of the bedfile:
	#   Intersection with very large files can cause trouble with excessive memory usage.
	#   A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm. 
	set sortTmpFile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_sort.tmp.bash"
	ReplaceTextInFile "#!/bin/bash" $sortTmpFile
	WriteTextInFile "# The locale specified by the environment can affects the traditional sort order. We need to use native byte values." $sortTmpFile
	WriteTextInFile "export LC_ALL=C" $sortTmpFile
	WriteTextInFile "sort -k1,1 -k2,2n $formattedFile > $formattedSortedFile" $sortTmpFile
	file attributes $sortTmpFile -permissions 0755
	if {[catch {eval exec bash $sortTmpFile} Message]} {
	    puts "-- checkCytoband --"
	    puts "sort -k1,1 -k2,2n $formattedFile > $formattedSortedFile"
	    puts "$Message"
	    puts "Exit with error"
	    exit 2
	}
	file delete -force $sortTmpFile 
	file delete -force $formattedFile
	file delete -force $cytobandBEDfile
    }

    # Number of annotation columns (without the 3 columns "chrom start end")
    set g_numberOfAnnotationCol($formattedSortedFile) 1

    # Cytoband annotation?
    if {[file exists $formattedSortedFile] && [file exists $cytobandHeaderFile]} {
	set g_AnnotSV(cytoband) 1
    } else {
	set g_AnnotSV(cytoband) 0
    }
}
