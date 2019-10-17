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


## - Check and create if necessary the following files:
##   'date'_IMH_DUP_SV.sorted.bed" 
##   'date'_IMH_DEL_SV.sorted.bed" 
##   'date'_IMH_INV_SV.sorted.bed" 
##   'date'_IMH_INS_SV.sorted.bed" 
proc checkIMHfile {} {

    global g_AnnotSV

    ## Check if IMH file has been downloaded then formatted
    #########################################################
    set IMHdir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/SVincludedInFt/IMH/$g_AnnotSV(genomeBuild)"
    set IMHfileDownloaded [glob -nocomplain "$IMHdir/B*.callset.public.bedpe"]
    set IMHfileFormattedAndSorted  [glob -nocomplain "$IMHdir/*_IMH_DUP_SV.sorted.bed"]

    if {$IMHfileDownloaded eq "" && $IMHfileFormattedAndSorted eq ""} {
	# No IMH annotation
	set g_AnnotSV(IMHann) 0
	return
    } else {
	set g_AnnotSV(IMHann) 1
	# Check if the user asked for these annotations in the configfile
	set test 0
	foreach col "IMH_ID IMH_AF IMH_ID_others" {
	    if {[lsearch -exact "$g_AnnotSV(outputColHeader)" $col] ne -1} {set test 1;break}
	}
	if {$test eq 0} {set g_AnnotSV(IMHann) 0; return}
    }

    foreach SVTYPE {DUP DEL INV INS} {
	set IMHfileFormattedAndSorted [glob -nocomplain "$IMHdir/*_IMH_${SVTYPE}_SV.sorted.bed"]
	if {[llength $IMHfileFormattedAndSorted]>1} {
	    puts "Several IMH files exist:"
	    puts "$IMHfileFormattedAndSorted"
	    puts "Keep only one: [lindex $IMHfileFormattedAndSorted end]\n"
	    foreach IMHf [lrange $IMHfileFormattedAndSorted 0 end-1] {
		file rename -force $IMHf $IMHf.notused
	    }
	}
    } 
 
    if {$IMHfileFormattedAndSorted eq ""} {
	# The downloaded file exist but not the formatted.
	## Create:
	##   'date'_IMH_DUP_SV.sorted.bed" 
	##   'date'_IMH_DEL_SV.sorted.bed" 
	##   'date'_IMH_INV_SV.sorted.bed" 
	##   'date'_IMH_INS_SV.sorted.bed" ; # Header: chr start end IMH_*_ID IMH_*_AF
	set IMH_DUP_fileFormatted "$IMHdir/[clock format [clock seconds] -format "%Y%m%d"]_IMH_DUP_SV.sorted.bed"
	set IMH_DEL_fileFormatted "$IMHdir/[clock format [clock seconds] -format "%Y%m%d"]_IMH_DEL_SV.sorted.bed"
	set IMH_INV_fileFormatted "$IMHdir/[clock format [clock seconds] -format "%Y%m%d"]_IMH_INV_SV.sorted.bed"
	set IMH_INS_fileFormatted "$IMHdir/[clock format [clock seconds] -format "%Y%m%d"]_IMH_INS_SV.sorted.bed"

	puts "...IMH configuration"

	puts "\t...creation of $IMHdir/[clock format [clock seconds] -format "%Y%m%d"]_IMH_*_SV.sorted.bed"
	puts "\t   (done only once during the first IMH annotation)\n"

	set L_TextToWrite(DUP) {}
	set L_TextToWrite(DEL) {}
	set L_TextToWrite(INV) {}
	set L_TextToWrite(INS) {}
	set f [open "$IMHfileDownloaded"]
	while {![eof $f]} {
	    set L [gets $f]
	    if {[string index $L 0] eq "#" || $L eq ""} {continue}
	    # Line example:
	    #1       10361   10590   15      102520301       102520447       2838    47911.3 -       -       BND     LOW     2838_1  N       [15:102520406[N 2838_2  N       [1:10567[N      SVTYPE=BND;POS=10567;STRANDS=--:97;IMPRECISE;CIPOS=-205,24;CIEND=-104,42;CIPOS95=-20,20;CIEND95=-6,6;MATEID=2838_2;EVENT=2838;SU=97;PE=97;SR=0;ALG=PROD;AF=0.07568;NSAMP=1274;MSQ=37.02;AN=16834;AC=1274;NS=8438;AC_Hom=0;AC_Het=1275;AC_Hemi=0;MAF=0.0755511;HWE=4.03332e-23;NFAM=1263     SVTYPE=BND;POS=102520406;STRANDS=--:97;IMPRECISE;CIPOS=-104,42;CIEND=-205,24;CIPOS95=-6,6;CIEND95=-20,20;MATEID=2838_1;EVENT=2838;SECONDARY;SU=97;PE=97;SR=0;ALG=PROD;AF=0.07568;NSAMP=1274;MSQ=37.02;AN=16834;AC=1274;NS=8438;AC_Hom=0;AC_Het=1275;AC_Hemi=0;MAF=0.0755511;HWE=4.03332e-23;NFAM=1263
	    set Ls [split $L "\t"] 
	    regsub "chr" [lindex $Ls 0] "" chrom
	    set start  [lindex $Ls 1]
	    regsub "chr" [lindex $Ls 3] "" chrom2 
	    set end    [lindex $Ls 5]
	    set SVTYPE [normalizeSVtype [lindex $Ls 10]]
	    if {[lsearch -exact [array names L_TextToWrite] "$SVTYPE"] eq -1} {continue}
	    set Annotations [lindex $Ls 18]
	    if {![regexp ";AF=(.+?);" $Annotations match AF]} {continue}
	    if {$chrom ne $chrom2} {continue}	    
	    lappend L_TextToWrite($SVTYPE) "$chrom\t$start\t$end\t$chrom:$start-${end}_$SVTYPE\t$AF"
	}
	close $f
	ReplaceTextInFile [join $L_TextToWrite(DUP) "\n"] $IMH_DUP_fileFormatted.tmp
	ReplaceTextInFile [join $L_TextToWrite(DEL) "\n"] $IMH_DEL_fileFormatted.tmp
	ReplaceTextInFile [join $L_TextToWrite(INV) "\n"] $IMH_INV_fileFormatted.tmp
	ReplaceTextInFile [join $L_TextToWrite(INS) "\n"] $IMH_INS_fileFormatted.tmp

	file delete -force $IMHfileDownloaded


	# Sorting of the bedfile:
	# Intersection with very large files can cause trouble with excessive memory usage.
	# A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm. 
	foreach SVTYPE [array names L_TextToWrite] {
	    if {[catch {eval exec sort -k1,1 -k2,2n [set IMH_${SVTYPE}_fileFormatted].tmp > [set IMH_${SVTYPE}_fileFormatted]} Message]} {
		puts "-- checkIMHfile --"
		puts "sort -k1,1 -k2,2n [set IMH_${SVTYPE}_fileFormatted].tmp > [set IMH_${SVTYPE}_fileFormatted]"
		puts "$Message"
		puts "Exit with error"
		exit 2
	    }
	    file delete -force [set IMH_${SVTYPE}_fileFormatted].tmp
	}
    }
}




# Return the overlapped SV annotation only with the IMH SV of the same SVtype
# If no SVtype is provided, no IMH annotation is returned (but IMH SV ID of different SVtype are stored in IMH_ID_others)
proc IMHannotation {SVchrom SVstart SVend SVtype L_i} {

    global g_AnnotSV
    global IMHtext

    # SVtype = DEL, DUP, INV or MEI (=INS)    
    # SVtype in which category: DUP? DEL? INV? INS? None?
    set SVtype [normalizeSVtype $SVtype]

    set IMHdir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/SVincludedInFt/IMH/$g_AnnotSV(genomeBuild)"
    
    if {![info exists IMHtext(DONE)]} {
	
	# headerOutput "IMH_ID IMH_AF IMH_ID_others"
	set L_IMHtext(Empty) "{} -1 {}"
	foreach i $L_i {
	    lappend IMHtext(Empty) "[lindex $L_IMHtext(Empty) $i]"
	}
	set IMHtext(Empty) [join $IMHtext(Empty) "\t"]

	foreach svtype {"DUP" "DEL" "INV" "INS"} {
	    # Intersect
	    set IMH_BEDfile [glob -nocomplain "$IMHdir/*_IMH_${svtype}_SV.sorted.bed"]
	    regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".intersect.IMH-$svtype" tmpFile
	    set tmpFile "$g_AnnotSV(outputDir)/[file tail $tmpFile]"
	    file delete -force $tmpFile	
	    if {[catch {exec $g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(fullAndSplitBedFile) -b $IMH_BEDfile -wa -wb > $tmpFile} Message]} {
		puts "-- IMHannotation, $svtype --"
		puts "$g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(fullAndSplitBedFile) -b $IMH_BEDfile -wa -wb > $tmpFile"
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

		set IMH_AF [lindex $Ls end]
		set IMH_ID [lindex $Ls end-1]
		set IMH_end [lindex $Ls end-2]
		set IMH_start [lindex $Ls end-3]

		# Select:
		# - IMH that share > XX% length with the SV to annotate
		# -> doesn't select the INSERTION (ALU, INS, LINE1, SVA: length < 50bp) if overlap < 70%
		set IMH_length [expr {$IMH_end-$IMH_start}]
		set SVtoAnn_length [expr {$SVtoAnn_end-$SVtoAnn_start}]
		# The IMH SV is an insertion or a breakpoint
		if {$IMH_length<1} {
		    set IMH_length 1
		} 	
		# The SV to annotate is an insertion or a breakpoint
		if {$SVtoAnn_length<1} {
		    set SVtoAnn_length 1
		} 	

		if {$SVtoAnn_start < $IMH_start} {
		    set overlap_start $IMH_start
		} else {
		    set overlap_start $SVtoAnn_start
		}
		if {$SVtoAnn_end < $IMH_end} {
		    set overlap_end $SVtoAnn_end
		} else {
		    set overlap_end $IMH_end
		}
		set overlap_length [expr {$overlap_end - $overlap_start}]
		
		# Keeping only IMH respecting the overlaps (reciprocal or not reciprocal)
		if {[expr {$overlap_length*100.0/$SVtoAnn_length}] < $g_AnnotSV(overlap)} {continue}
		if {$g_AnnotSV(reciprocal) eq "yes"} {			
		    if {[expr {$overlap_length*100.0/$IMH_length}] < $g_AnnotSV(overlap)} {continue}
		}	  
		
		set SVtoAnn "$SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end"
		lappend L_allSVtoAnn               $SVtoAnn
		lappend L_SVtoAnn($svtype)         $SVtoAnn
		lappend L_IMH_ID($SVtoAnn,$svtype) $IMH_ID
		lappend L_IMH_AF($SVtoAnn,$svtype) $IMH_AF
	    }
	    file delete -force $tmpFile
	}
	    
	# Loading IMH final annotation for each SV
	foreach svtype {"DUP" "DEL" "INV" "INS" "None"} {
	    if {[info exists L_allSVtoAnn]} {
		foreach SVtoAnn [lsort -unique $L_allSVtoAnn] {
		    if {[info exists L_SVtoAnn($svtype)] && [lsearch -exact $L_SVtoAnn($svtype) "$SVtoAnn"] ne -1} {
			# Change metrics from "." to ","
			if {[set g_AnnotSV(metrics)] eq "fr"} {
			    regsub -all {\.} $L_IMH_AF($SVtoAnn,$svtype) "," L_IMH_AF($SVtoAnn,$svtype)
			}
			# For L_IMH_AF($SVtoAnn): we keep only the maximum value
			set max_AF -1
			foreach a $L_IMH_AF($SVtoAnn,$svtype) {if {$a ne "NA" && $a > $max_AF} {set max_AF $a}}
			
			# headerOutput "IMH_ID IMH_AF IMH_ID_others"
			set L_IMHtext($SVtoAnn,$svtype) ""
			lappend L_IMHtext($SVtoAnn,$svtype) "[join $L_IMH_ID($SVtoAnn,$svtype) ";"]"
			lappend L_IMHtext($SVtoAnn,$svtype) "$max_AF"
		    } else {
			set L_IMHtext($SVtoAnn,$svtype) "{} -1"
		    }
		    # IMH_ID_other (identifiers of IMH SV with different SV type)
		    set L_IMH_ID_other($SVtoAnn,$svtype) ""
		    foreach othersvtype {"DUP" "DEL" "INV" "INS"} {
			if {$othersvtype eq $svtype} {continue}
			catch {lappend L_IMH_ID_other($SVtoAnn,$svtype) {*}$L_IMH_ID($SVtoAnn,$othersvtype)}
		    }
		    lappend L_IMHtext($SVtoAnn,$svtype) "[join $L_IMH_ID_other($SVtoAnn,$svtype) ";"]"

		    # Keep only the user requested columns (defined in the configfile)
		    set IMHtext($SVtoAnn,$svtype) ""
		    foreach i $L_i {
			lappend IMHtext($SVtoAnn,$svtype) "[lindex $L_IMHtext($SVtoAnn,$svtype) $i]"
		    }
		    set IMHtext($SVtoAnn,$svtype) [join $IMHtext($SVtoAnn,$svtype) "\t"]
		}
	    }
	}
	catch {unset L_IMHtext}
	set IMHtext(DONE) 1	
    }
    
    if {[info exist IMHtext($SVchrom,$SVstart,$SVend,$SVtype)]} {
	return $IMHtext($SVchrom,$SVstart,$SVend,$SVtype)
    } else {
	return $IMHtext(Empty)
    }
}
