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


## - Check and create if necessary the 'date'_dbVar_pathogenic_NR_SV.formatted.sorted.bed" file.
proc checkPathogenicNRSVfile {} {

    global g_AnnotSV


    ## Check if dbVar files has been formatted (if yes, nothing to be done)
    #######################################################################
    set NRSVdir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/FtIncludedInSV/dbVar_pathogenic_NR_SV/$g_AnnotSV(genomeBuild)"
    set dbVarBedFileFormattedSorted {}
    foreach f [glob -nocomplain "$NRSVdir/*_dbVar_pathogenic_NR_SV.formatted.sorted.bed"] {
	regsub -nocase ".formatted.sorted.bed$" $f ".header.tsv" newf
	if {![file exists $newf]} {continue}
	lappend dbVarBedFileFormattedSorted $f
    }

    set g_AnnotSV(NRSVann) 0
    set n [llength $dbVarBedFileFormattedSorted]
    if {$n>1} {
        puts "Several dbVar pathogenic NR SV files exist:"
        puts "$dbVarBedFileFormattedSorted"
        puts "Keep only one: [lindex [lsort $dbVarBedFileFormattedSorted] end]\n"
        foreach f [lrange $dbVarBedFileFormattedSorted 0 end-1] {
            file rename -force $f $f.notused
        }
	set g_AnnotSV(NRSVann) 1
	# Check if the user asked for these annotations in the configfile
	set test 0
	foreach col "dbVar_event dbVar_variant dbVar_status" {
	    if {[lsearch -exact "$g_AnnotSV(outputColHeader)" $col] ne -1} {set test 1;break}
	}
	if {$test eq 0} {set g_AnnotSV(NRSVann) 0; return}

	return 1
    } elseif {$n eq 1} {
	set g_AnnotSV(NRSVann) 1
	# Check if the user asked for these annotations in the configfile
	set test 0
	foreach col "dbVar_event dbVar_variant dbVar_status" {
	    if {[lsearch -exact "$g_AnnotSV(outputColHeader)" $col] ne -1} {set test 1;break}
	}
	if {$test eq 0} {set g_AnnotSV(NRSVann) 0; return}
	return 1
    }


    ## Formatting of the dbVar files
    ################################
    set NRSVdir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/FtIncludedInSV/dbVar_pathogenic_NR_SV/$g_AnnotSV(genomeBuild)"
    set dbVarFileDownloaded1 "$NRSVdir/$g_AnnotSV(genomeBuild).nr_deletions.tsv.gz"
    set dbVarFileDownloaded2 "$NRSVdir/$g_AnnotSV(genomeBuild).nr_duplications.tsv.gz"

    if {[file exists $dbVarFileDownloaded1] && [file exists $dbVarFileDownloaded2]} {
	## The 2 downloaded files exist but have not been formatted.
	## Create:
	##   - 'date'_dbVar_pathogenic_NR_SV.bed
	##   - 'date'_dbVar_pathogenic_NR_SV.header.tsv  
	set dbVarBedFile "$NRSVdir/[clock format [clock seconds] -format "%Y%m%d"]_dbVar_pathogenic_NR_SV.bed"
	set dbVarBedFileFormatted "$NRSVdir/[clock format [clock seconds] -format "%Y%m%d"]_dbVar_pathogenic_NR_SV.formatted.bed"
	set dbVarBedFileFormattedSorted "$NRSVdir/[clock format [clock seconds] -format "%Y%m%d"]_dbVar_pathogenic_NR_SV.formatted.sorted.bed"
	set dbVarHeaderFile "$NRSVdir/[clock format [clock seconds] -format "%Y%m%d"]_dbVar_pathogenic_NR_SV.header.tsv"

	puts "...dbVar configuration"

	puts "\t...creation of [file tail $dbVarBedFileFormattedSorted] and [file tail $dbVarHeaderFile]"
	puts "\t   (done only once during the first dbVar annotation)"

	## Header from tsv downloaded files:
	## #chr    outermost_start outermost_stop  variant_count   variant_type    method  analysis        platform        study   variant clinical_assertion      clinvar_accession       bin_size
	set L_bedTextToWrite {}
	set headerTextToWrite "#chr\toutermost_start\toutermost_stop\tdbVar_event\tdbVar_variant\tdbVar_status"
	foreach dbVarFile "$dbVarFileDownloaded1 $dbVarFileDownloaded2" {
	    puts "\t   ...reading $dbVarFile"
	    set i_type -1
	    set i_variant -1
	    set i_clinvar -1
	    set f [open "|gzip -cd $dbVarFile"]
	    while {![eof $f]} {

		set L [gets $f]
		set Ls [split $L "\t"] 

		if {[string range $L 0 3] eq "#chr"} {
		    set i_type [lsearch -exact "$Ls" "variant_type"]
		    if {$i_type eq -1} {
			puts "\tWARNING: dbVar header doesn't contain a \"variant_type\" annotation."
			puts "\tNo dbVar annotation can be done"
			return 1
		    }
		    set i_variant [lsearch -exact "$Ls" "variant"]
		    if {$i_variant eq -1} {
			puts "\tWARNING: dbVar header doesn't contain a \"variant\" annotation."
			puts "\tNo dbVar annotation can be done"
			return 1
		    }
		    set i_clinvar [lsearch -exact "$Ls" "clinical_assertion"]
		    if {$i_clinvar eq -1} {
			puts "\tWARNING: dbVar header doesn't contain a \"clinical_assertion\" annotation."
			puts "\tNo dbVar annotation can be done"
			return 1
		    }
		    continue
		}
		if {[string index $L 0] eq "#" || $L eq ""} {continue}
		
		set clinic [lindex $Ls $i_clinvar]
		if {![regexp -nocase "pathogenic" $clinic]} {continue}
		set type [lindex $Ls $i_type]
		set variant [lindex $Ls $i_variant]
		lappend L_bedTextToWrite "[join [lrange $Ls 0 2] "\t"]\t$type\t$variant\t$clinic"
		
	    }
	    close $f
	}

	ReplaceTextInFile "$headerTextToWrite" $dbVarHeaderFile
	# Creation of dbVarBedFile
	ReplaceTextInFile [join $L_bedTextToWrite "\n"] $dbVarBedFile	
	# Creation of dbVarBedFileFormatted
	checkBed $dbVarBedFile [file dirname $dbVarBedFile]
	file delete -force $dbVarBedFile
	# Creation of dbVarBedFileFormattedSorted
	# Sorting of the bedfile:
	# Intersection with very large files can cause trouble with excessive memory usage.
	# A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm. 
	if {[catch {eval exec sort -k1,1 -k2,2n $dbVarBedFileFormatted > $dbVarBedFileFormattedSorted} Message]} {
	    puts "-- checkPathogenicNRSVfile --"
	    puts "sort -k1,1 -k2,2n $dbVarBedFileFormatted > $dbVarBedFileFormattedSorted"
	    puts "$Message"
	    puts "Exit with error"
	    exit 2
	}
	file delete -force $dbVarBedFileFormatted

	puts ""

	file delete -force $dbVarFileDownloaded1 
	file delete -force $dbVarFileDownloaded2

	set g_AnnotSV(NRSVann) 1
    }
}




proc pathogenicNRSVannotation {SVchrom SVstart SVend L_i} {

    global g_AnnotSV
    global pathogenicNRSVtext

    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)"
    set pathogenicNRSVfileFormatted [glob -nocomplain "$extannDir/FtIncludedInSV/dbVar_pathogenic_NR_SV/$g_AnnotSV(genomeBuild)/*_dbVar_pathogenic_NR_SV.formatted.sorted.bed"] 
    
    if {![info exists pathogenicNRSVtext(DONE)]} {
	
	# headerOutput "dbVar_event\tdbVar_variant\tdbVar_status"
	set L_pathogenicNRSVtext(Empty) "{} {} {}"
	foreach i $L_i {
	    lappend pathogenicNRSVtext(Empty) "[lindex $L_pathogenicNRSVtext(Empty) $i]"
	}
	set pathogenicNRSVtext(Empty) [join $pathogenicNRSVtext(Empty) "\t"]
	
	# Intersect
	regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".intersect.NRSV" tmpFile
	set tmpFile "$g_AnnotSV(outputDir)/[file tail $tmpFile]"
	file delete -force $tmpFile	
	# -F .7 <=> at least 70% of B is overlapped
	# Keeping only NRSV with > 70% (default) overlaped with the SV
	# (if SV_to_annotate is an insertion, no overlap is retained)
	# No reciprocal overlap used for NR-SV (it would not make sense)
	set p [expr {$g_AnnotSV(overlap)*1.0/100}]
	if {[catch {exec $g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(fullAndSplitBedFile) -b $pathogenicNRSVfileFormatted -F $p -wa -wb > $tmpFile} Message]} {
	    puts "-- pathogenicNRSVannotation --"
	    puts "$g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(fullAndSplitBedFile) -b $pathogenicNRSVfileFormatted -F $p -wa -wb > $tmpFile"
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

	    set dbVar_status [lindex $Ls end]
	    set dbVar_variant [lindex $Ls end-1]
	    set dbVar_event  [lindex $Ls end-2]
	    set NRSV_end [lindex $Ls end-3]
	    set NRSV_start [lindex $Ls end-4]


##	    # Select:
##	    # - NRSV that share > XX % length with the SV to annotate
##	    set NRSV_length [expr {$NRSV_end-$NRSV_start}]
##	    set SVtoAnn_length [expr {$SVtoAnn_end - $SVtoAnn_start}]
##	    if {$SVtoAnn_length<=1} {
##		# SV to annotate is an insertion
##		set SVtoAnn_length 1
##	    } else {
##		# SV to annotate is not an insertion
##		if {$SVtoAnn_start < $NRSV_start} {
##		    set overlap_start $NRSV_start
##		} else {
##		    set overlap_start $SVtoAnn_start
##		}
##		if {$SVtoAnn_end < $NRSV_end} {
##		    set overlap_end $SVtoAnn_end
##		} else {
##		    set overlap_end $NRSV_end
##		}
##		set overlap_length [expr {$overlap_end - $overlap_start}]
##		
##		# Keeping only NRSV with > 70% (default) overlaped with the SV
##		if {[expr {$overlap_length*100.0/$NRSV_length}] < $g_AnnotSV(overlap)} {continue}
##	    }
##	    
	    lappend L_event($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end)   {*}[split $dbVar_event ";"]
	    lappend L_variant($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end) {*}[split $dbVar_variant ";"]
	    lappend L_status($SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end)  {*}[split $dbVar_status ";"]

	}
	

	# Loading NRSV final annotation for each SV
	foreach SVtoAnn [array names L_event] {
	    # For L_event($SVtoAnn): we keep a unique list of the events
	    set L_event($SVtoAnn) [lsort -unique $L_event($SVtoAnn)]
	    # For L_variant($SVtoAnn): we keep a unique list of the variant
	    set L_variant($SVtoAnn) [lsort -unique $L_variant($SVtoAnn)]
	    # For L_status($SVtoAnn): we keep a unique list of the status
	    set L_status($SVtoAnn) [lsort -unique $L_status($SVtoAnn)]

	    set L_pathogenicNRSVtext($SVtoAnn) ""
	    lappend L_pathogenicNRSVtext($SVtoAnn) "[join $L_event($SVtoAnn) ";"]"
	    lappend L_pathogenicNRSVtext($SVtoAnn) "[join $L_variant($SVtoAnn) ";"]"
	    lappend L_pathogenicNRSVtext($SVtoAnn) "[join $L_status($SVtoAnn) ";"]"

	    # Keep only the user requested columns (defined in the configfile)
	    set pathogenicNRSVtext($SVtoAnn) ""
	    foreach i $L_i {
		lappend pathogenicNRSVtext($SVtoAnn) "[lindex $L_pathogenicNRSVtext($SVtoAnn) $i]"
	    }
	    set pathogenicNRSVtext($SVtoAnn) [join $pathogenicNRSVtext($SVtoAnn) "\t"]
	}

	catch {unset L_pathogenicNRSVtext}
	set pathogenicNRSVtext(DONE) 1	
	file delete -force $tmpFile
    }
    
    if {[info exist pathogenicNRSVtext($SVchrom,$SVstart,$SVend)]} {
	return $pathogenicNRSVtext($SVchrom,$SVstart,$SVend)
    } else {
	return $pathogenicNRSVtext(Empty)
    }
}
