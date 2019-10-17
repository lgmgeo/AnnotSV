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


## - Check and create if necessary the following file:
##   'date'_gnomAD_DUP_SV.sorted.bed"
##   'date'_gnomAD_DEL_SV.sorted.bed" 
##   'date'_gnomAD_INV_SV.sorted.bed" 
##   'date'_gnomAD_INS_SV.sorted.bed" 
proc checkgnomADfile {} {

    global g_AnnotSV

    ## Check if gnomAD file has been downloaded then formatted
    #########################################################
    set gnomADdir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/SVincludedInFt/gnomAD/$g_AnnotSV(genomeBuild)"
    set gnomADfileDownloaded [glob -nocomplain "$gnomADdir/gnomad_v2_sv.sites.bed.gz"]
    set gnomADfileFormattedAndSorted  [glob -nocomplain "$gnomADdir/*_gnomAD_DUP_SV.sorted.bed"]

    if {$gnomADfileDownloaded eq "" && $gnomADfileFormattedAndSorted eq ""} {
	# No gnomAD annotation
	set g_AnnotSV(gnomADann) 0
	return
    } else {
	set g_AnnotSV(gnomADann) 1
	# Check if the user asked for these annotations in the configfile
	set test 0
	foreach col "GD_ID GD_AN GD_N_HET GD_N_HOMALT GD_AF GD_POPMAX_AF GD_ID_others" {
	    if {[lsearch -exact "$g_AnnotSV(outputColHeader)" $col] ne -1} {set test 1;break}
	}
	if {$test eq 0} {set g_AnnotSV(gnomADann) 0; return}
    }
    
    foreach SVTYPE {DUP DEL INV INS} {
	set gnomADfileFormattedAndSorted [glob -nocomplain "$gnomADdir/*_gnomAD_${SVTYPE}_SV.sorted.bed"]
	if {[llength $gnomADfileFormattedAndSorted]>1} {
	    puts "Several gnomAD files exist:"
	    puts "$gnomADfileFormattedAndSorted"
	    puts "Keep only one: [lindex $gnomADfileFormattedAndSorted end]\n"
	    foreach gnomADf [lrange $gnomADfileFormattedAndSorted 0 end-1] {
		file rename -force $gnomADf $gnomADf.notused
	    }
	}
    }

    if {$gnomADfileFormattedAndSorted eq ""} {
	# The downloaded file exist but not the formatted.
	## Create:
	##   'date'_gnomAD_DUP_SV.sorted.bed" 
	##   'date'_gnomAD_DEL_SV.sorted.bed" 
	##   'date'_gnomAD_INV_SV.sorted.bed" 
	##   'date'_gnomAD_INS_SV.sorted.bed" 
	set gnomAD_DUP_fileFormatted "$gnomADdir/[clock format [clock seconds] -format "%Y%m%d"]_gnomAD_DUP_SV.sorted.bed"
	set gnomAD_DEL_fileFormatted "$gnomADdir/[clock format [clock seconds] -format "%Y%m%d"]_gnomAD_DEL_SV.sorted.bed"
	set gnomAD_INV_fileFormatted "$gnomADdir/[clock format [clock seconds] -format "%Y%m%d"]_gnomAD_INV_SV.sorted.bed"
	set gnomAD_INS_fileFormatted "$gnomADdir/[clock format [clock seconds] -format "%Y%m%d"]_gnomAD_INS_SV.sorted.bed"

	puts "...gnomAD configuration"

	puts "\t...creation of $gnomADdir/[clock format [clock seconds] -format "%Y%m%d"]_gnomAD_*_SV.sorted.bed"
	puts "\t   (done only once during the first gnomAD annotation)\n"

	set L_TextToWrite(DUP) {}
	set L_TextToWrite(DEL) {}
	set L_TextToWrite(INV) {}
	set L_TextToWrite(INS) {}
	set f [open "| gzip -cd $gnomADfileDownloaded"]

	while {![eof $f]} {
	    set L [gets $f]
	    if {$L eq ""} {continue}
	    set Ls [split $L "\t"] 

	    if {[regexp "^#CHROM" $L]} {
		set i_chrom      [lsearch -exact $Ls "#CHROM"];       if {$i_chrom == -1} {puts "Bad header line syntax. #CHROM column not found - Exit with error"; exit 2}
		set i_start      [lsearch -exact $Ls "START"];        if {$i_start == -1} {puts "Bad header line syntax. START column not found - Exit with error"; exit 2}
		set i_end        [lsearch -exact $Ls "END"];          if {$i_end == -1} {puts "Bad header line syntax. END column not found - Exit with error"; exit 2}
		set i_svid       [lsearch -exact $Ls "NAME"];         if {$i_svid == -1} {puts "Bad header line syntax. NAME column not found - Exit with error"; exit 2}
		set i_svtype     [lsearch -exact $Ls "SVTYPE"];       if {$i_svtype == -1} {puts "Bad header line syntax. SVTYPE column not found - Exit with error"; exit 2}
		set i_an         [lsearch -exact $Ls "AN"];           if {$i_an == -1} {puts "Bad header line syntax. AN column not found - Exit with error"; exit 2}
		set i_nhet       [lsearch -exact $Ls "N_HET"];        if {$i_nhet == -1} {puts "Bad header line syntax. N_HET column not found - Exit with error"; exit 2}
		set i_nhomalt    [lsearch -exact $Ls "N_HOMALT"];     if {$i_nhomalt == -1} {puts "Bad header line syntax. N_HOMALT column not found - Exit with error"; exit 2}
		set i_af         [lsearch -exact $Ls "AF"];           if {$i_af == -1} {puts "Bad header line syntax. AF column not found - Exit with error"; exit 2}
		set i_popmaxaf   [lsearch -exact $Ls "POPMAX_AF"];    if {$i_popmaxaf == -1} {puts "Bad header line syntax. POPMAX_AF column not found - Exit with error"; exit 2}
		continue
	    }

	    set chrom      [lindex $Ls $i_chrom]
	    set start      [lindex $Ls $i_start]
	    set end        [lindex $Ls $i_end]
	    set svid       [lindex $Ls $i_svid]
	    set SVTYPE     [lindex $Ls $i_svtype]
	    # WARNING:
	    # - MCNV have multiple coma-separated-values for nhet, nhomalt, af and popmaxaf... 
	    #   => We need to keep only 1 value (for the ranking)
	    # - Look at the "gnomad_v2_sv.sites.bed" file:
	    #   445858 lines of which only 1148 are MCNV
	    #   => We remove these SV from the annotation files
	    #
	    # if {$SVTYPE == "MCNV"} {set SVTYPE "DUP"} ;# MCNV = multiallelic CNV = CNV2, CNV3...
	    if {[lsearch -exact {DEL DUP INV INS} "$SVTYPE"] eq -1} {continue}
	    set an         [lindex $Ls $i_an]
	    set nhet       [lindex $Ls $i_nhet]
	    set nhomalt    [lindex $Ls $i_nhomalt]
	    set af         [lindex $Ls $i_af]
	    set popmaxaf   [lindex $Ls $i_popmaxaf]
	    lappend L_TextToWrite($SVTYPE) "$chrom\t$start\t$end\t$svid\t$an\t$nhet\t$nhomalt\t$af\t$popmaxaf"
	}
	close $f
	ReplaceTextInFile [join $L_TextToWrite(DUP) "\n"] $gnomAD_DUP_fileFormatted.tmp
	ReplaceTextInFile [join $L_TextToWrite(DEL) "\n"] $gnomAD_DEL_fileFormatted.tmp
	ReplaceTextInFile [join $L_TextToWrite(INV) "\n"] $gnomAD_INV_fileFormatted.tmp
	ReplaceTextInFile [join $L_TextToWrite(INS) "\n"] $gnomAD_INS_fileFormatted.tmp

	file delete -force $gnomADfileDownloaded


	# Sorting of the bedfile:
	# Intersection with very large files can cause trouble with excessive memory usage.
	# A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm. 
	foreach SVTYPE [array names L_TextToWrite] {
	    if {[catch {eval exec sort -k1,1 -k2,2n [set gnomAD_${SVTYPE}_fileFormatted].tmp > [set gnomAD_${SVTYPE}_fileFormatted]} Message]} {
		puts "-- checkgnomADfile --"
		puts "sort -k1,1 -k2,2n [set gnomAD_${SVTYPE}_fileFormatted].tmp > [set gnomAD_${SVTYPE}_fileFormatted]"
		puts "$Message"
		puts "Exit with error"
		exit 2
	    }
	    file delete -force [set gnomAD_${SVTYPE}_fileFormatted].tmp
	}
    }
}



# Return the overlapped SV annotation only with the gnomAD SV of the same SVtype
# If no SVtype is provided, no gnomAD annotation is returned (but gnomAD SV ID of different SVtype are stored in GD_ID_others)
proc gnomADannotation {SVchrom SVstart SVend SVtype L_i} {
    
    global g_AnnotSV
    global gnomADtext
    
    # SVTYPE = DUP, DEL, INV or INS
    
    # SVtype in which category: DUP? DEL? INV? INS? None?
    set SVtype [normalizeSVtype $SVtype]


    set gnomADdir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/SVincludedInFt/gnomAD/$g_AnnotSV(genomeBuild)"

    if {![info exists gnomADtext(DONE)]} {
	
	# headerOutput "GD_ID	GD_AN	GD_N_HET	GD_N_HOMALT	GD_AF	GD_POPMAX_AF	GD_ID_others"
	set L_gnomADtext(Empty) "{} {} {} {} {-1} {-1} {}"
	foreach i $L_i {
	    lappend gnomADtext(Empty) "[lindex $L_gnomADtext(Empty) $i]"
	}
	set gnomADtext(Empty) "[join $gnomADtext(Empty) "\t"]"
	
	foreach svtype {"DUP" "DEL" "INV" "INS"} {
	    # Intersect
	    set gnomAD_BEDfile [glob -nocomplain "$gnomADdir/*_gnomAD_${svtype}_SV.sorted.bed"]
	    regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".intersect.gnomAD-$svtype" tmpFile
	    set tmpFile "$g_AnnotSV(outputDir)/[file tail $tmpFile]"
	    file delete -force $tmpFile	
	    if {[catch {exec $g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(fullAndSplitBedFile) -b $gnomAD_BEDfile -wa -wb > $tmpFile} Message]} {
		puts "-- gnomADannotation --"
		puts "$g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(fullAndSplitBedFile) -b $gnomAD_BEDfile -wa -wb > $tmpFile"
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

		set GD_POPMAX_AF   [lindex $Ls end]
		set GD_AF          [lindex $Ls end-1]
		set GD_N_HOMALT    [lindex $Ls end-2]
		set GD_N_HET       [lindex $Ls end-3]
		set GD_AN          [lindex $Ls end-4]
		set GD_ID          [lindex $Ls end-5]
		set gnomAD_end     [lindex $Ls end-6]
		set gnomAD_start   [lindex $Ls end-7]

		# Select:
		# - gnomAD that share > XX% length with the SV to annotate
		# -> doesn't select the INSERTION (ALU, INS, LINE1, SVA: length < 50bp) if overlap < 70%
		set gnomAD_length [expr {$gnomAD_end-$gnomAD_start}]
		set SVtoAnn_length [expr {$SVtoAnn_end-$SVtoAnn_start}]
		# The gnomAD SV is an insertion or a breakpoint
		if {$gnomAD_length<1} {
		    set gnomAD_length 1
		} 	
		# The SV to annotate is an insertion or a breakpoint
		if {$SVtoAnn_length<1} {
		    set SVtoAnn_length 1
		} 	
		
		if {$SVtoAnn_start < $gnomAD_start} {
		    set overlap_start $gnomAD_start
		} else {
		    set overlap_start $SVtoAnn_start
		}
		if {$SVtoAnn_end < $gnomAD_end} {
		    set overlap_end $SVtoAnn_end
		} else {
		    set overlap_end $gnomAD_end
		}
		set overlap_length [expr {$overlap_end - $overlap_start}]
		
		# Keeping only gnomAD respecting the overlaps (reciprocal or not reciprocal)
		if {[expr {$overlap_length*100.0/$SVtoAnn_length}] < $g_AnnotSV(overlap)} {continue}
		if {$g_AnnotSV(reciprocal) eq "yes"} {			
		    if {[expr {$overlap_length*100.0/$gnomAD_length}] < $g_AnnotSV(overlap)} {continue}
		}	 
		
		set SVtoAnn "$SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end"
		lappend L_allSVtoAnn                       $SVtoAnn
		lappend L_SVtoAnn($svtype)                 $SVtoAnn
		lappend L_GD_ID($SVtoAnn,$svtype)          $GD_ID
		lappend L_GD_AN($SVtoAnn,$svtype)          $GD_AN
		lappend L_GD_N_HET($SVtoAnn,$svtype)       $GD_N_HET
		lappend L_GD_N_HOMALT($SVtoAnn,$svtype)    $GD_N_HOMALT
		lappend L_GD_AF($SVtoAnn,$svtype)          $GD_AF
		lappend L_GD_POPMAX_AF($SVtoAnn,$svtype)   $GD_POPMAX_AF
	    }
	    file delete -force $tmpFile
	}

	# Loading gnomAD final annotation for each SV
	foreach svtype {"DUP" "DEL" "INV" "INS" "None"} {
	    if {[info exists L_allSVtoAnn]} {
		foreach SVtoAnn [lsort -unique $L_allSVtoAnn] {
		    if {[info exists L_SVtoAnn($svtype)] && [lsearch -exact $L_SVtoAnn($svtype) "$SVtoAnn"] ne -1} {
			# Change metrics from "." to ","
			if {[set g_AnnotSV(metrics)] eq "fr"} {
			    regsub -all {\.} $L_GD_POPMAX_AF($SVtoAnn,$svtype)   "," L_GD_POPMAX_AF($SVtoAnn,$svtype)
			    regsub -all {\.} $L_GD_AF($SVtoAnn,$svtype)          "," L_GD_AF($SVtoAnn,$svtype)
			}
			# For AF, we keep only the maximum value
			set max_AF -1
			set max_POPMAX -1
			foreach a $L_GD_AF($SVtoAnn,$svtype)          {if {$a ne "NA" && $a > $max_AF}     {set max_AF $a}}
			foreach a $L_GD_POPMAX_AF($SVtoAnn,$svtype)   {if {$a ne "NA" && $a > $max_POPMAX} {set max_POPMAX $a}}
			
			# headerOutput "GD_ID GD_AN GD_N_HET GD_N_HOMALT GD_AF GD_POPMAX_AF GD_ID_others"
			set L_gnomADtext($SVtoAnn,$svtype) ""
			lappend L_gnomADtext($SVtoAnn,$svtype) "[join $L_GD_ID($SVtoAnn,$svtype) ";"]"
			lappend L_gnomADtext($SVtoAnn,$svtype) "[join $L_GD_AN($SVtoAnn,$svtype) ";"]" 
			lappend L_gnomADtext($SVtoAnn,$svtype) "[join $L_GD_N_HET($SVtoAnn,$svtype) ";"]" 
			lappend L_gnomADtext($SVtoAnn,$svtype) "[join $L_GD_N_HOMALT($SVtoAnn,$svtype) ";"]"
			lappend L_gnomADtext($SVtoAnn,$svtype) "$max_AF" 
			lappend L_gnomADtext($SVtoAnn,$svtype) "$max_POPMAX"
		    } else {
			set L_gnomADtext($SVtoAnn,$svtype) "{} {} {} {} {-1} {-1}"
		    }
		    # GD_ID_other (identifiers of gnomAD SV with different SV type)
		    set L_GD_ID_other($SVtoAnn,$svtype) ""
		    foreach othersvtype {"DUP" "DEL" "INV" "INS"} {
			if {$othersvtype eq $svtype} {continue}
			catch {lappend L_GD_ID_other($SVtoAnn,$svtype) {*}$L_GD_ID($SVtoAnn,$othersvtype)}
		    }
		    lappend L_gnomADtext($SVtoAnn,$svtype) [join $L_GD_ID_other($SVtoAnn,$svtype) ";"]
		    # Keep only the user requested columns (defined in the configfile)
		    set gnomADtext($SVtoAnn,$svtype) ""
		    foreach i $L_i {
			lappend gnomADtext($SVtoAnn,$svtype) "[lindex $L_gnomADtext($SVtoAnn,$svtype) $i]"
		    }
		    set gnomADtext($SVtoAnn,$svtype) [join $gnomADtext($SVtoAnn,$svtype) "\t"]
		}
	    }
	}
	catch {unset L_gnomADtext}
	set gnomADtext(DONE) 1	
    }


    if {[info exist gnomADtext($SVchrom,$SVstart,$SVend,$SVtype)]} {
	return $gnomADtext($SVchrom,$SVstart,$SVend,$SVtype)
    } else {
	return $gnomADtext(Empty)
    }
}
