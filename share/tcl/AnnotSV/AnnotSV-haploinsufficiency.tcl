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



## - Check if the following HI file has been downloaded:
#    - HI_Predictions_Version3.bed.gz
#  
## - Check and create if necessary the following file:
#    - 'date'_HI.tsv.gz
proc checkHIfile {} {

    global g_AnnotSV

    ## Check if the HI file has been downloaded then formatted
    #########################################################
    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Genes-based"
 
    set HIfileDownloaded [glob -nocomplain "$extannDir/DDD/HI_Predictions_Version3.bed.gz"]
    set HIfileFormattedGzip [glob -nocomplain "$extannDir/DDD/*_HI.tsv.gz"] 

    if {$HIfileDownloaded eq "" && $HIfileFormattedGzip eq ""} {
	# No "Haploinsufficiency" annotation
	return
    } 

    if {[llength $HIfileFormattedGzip]>1} {
	puts "Several Haploinsufficiency files exist:"
	puts "$HIfileFormattedGzip"
	puts "Keep only one: [lindex $HIfileFormattedGzip end]\n"
	foreach hi [lrange $HIfileFormattedGzip 0 end-1] {
	    file rename -force $hi $hi.notused
	}
	return
    } 

    if {$HIfileFormattedGzip eq ""} {    
	## Create : 'date'_HI.tsv     
	# Header: chr start end syn_z mis_z pLI 

	set HIfileFormatted "$extannDir/DDD/[clock format [clock seconds] -format "%Y%m%d"]_HI.tsv"

	puts "...HI configuration ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"

	puts "\t...creation of $HIfileFormatted.gz"
	puts "\t   (done only once during the first HI annotation)\n"

	set TexteToWrite {genes\tHI_DDDpercent}
	foreach L [LinesFromGZFile $HIfileDownloaded] {
	    if {[regexp  "^track name" $L]} {continue}
	    set Ls [split $L "\t"]
	    
	    # Line example:
	    # chrX    37850070        37850569        HYPM|0.000043663|100%   0.000043663     .       37850070        37850569        0,255,0
	    set infos [split [lindex $Ls 3] "|"] 
	    set gene [lindex $infos 0]   
	    regsub "%" [lindex $infos end] "" percent  
	    lappend TexteToWrite "$gene\t$percent"
	} 

	# Write outputfile
	WriteTextInFile [join $TexteToWrite "\n"] $HIfileFormatted
	if {[catch {exec gzip $HIfileFormatted} Message]} {
	    puts "-- checkHIfile --"
	    puts "gzip $HIfileFormatted"
	    puts "$Message\n"
	}

	file delete -force $HIfileDownloaded
    }
}

