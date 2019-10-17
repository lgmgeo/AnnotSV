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


## - Check if the following GH files has been downloaded:
##   - GeneHancer_elements.txt
##   - GeneHancer_gene_associations_scores.txt
##   - GeneHancer_hg19.txt
##   - GeneHancer_tissues.txt
##  
## - Check and create if necessary:
##   - 'date'_GH_GRCh37.sorted.bed
##   - 'date'_GH_GRCh38.sorted.bed
proc checkGHfiles {} {

    global g_AnnotSV


    ## Check if GH files has been formatted
    #######################################
    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)"
    set GHdir "$extannDir/FtIncludedInSV/GeneHancer"
    set GHfileFormatted [glob -nocomplain "$GHdir/$g_AnnotSV(genomeBuild)/*_GH_$g_AnnotSV(genomeBuild).sorted.bed"] 

    if {[llength $GHfileFormatted]>1} {
	set g_AnnotSV(GHann) 1
	puts "Several GH files exist:"
	puts "$GHfileFormatted"
	puts "Keep only one: [lindex $GHfileFormatted end]\n"
	foreach GH [lrange $GHfileFormatted 0 end-1] {
	    file rename -force $GH $GH.notused
	}
	return
    } 

    if {[llength $GHfileFormatted] eq 1} {
	set g_AnnotSV(GHann) 1
	return
    }

    ###############################################
    ## No formatted GH files available at this step
    ###############################################

    ## Check if GH file have been downloaded and unzipped
    #####################################################
    # Checks if the 4 unzipped GH files exist:
    set GHelementsF "$GHdir/GeneHancer_elements.txt"
    set GHassociationsF "$GHdir/GeneHancer_gene_associations_scores.txt"
    set GHhg19F "$GHdir/GeneHancer_hg19.txt"
    set GHtissuesF "$GHdir/GeneHancer_tissues.txt"
    foreach GHfile "$GHelementsF $GHassociationsF $GHhg19F $GHtissuesF" {
	if {![file exists $GHfile]} {
	    set g_AnnotSV(GHann) 0
	    if {$g_AnnotSV(organism) eq "Human"} {
		puts "\nWARNING: No GeneHancer annotations available."
		puts "(Please, see in the README file how to add these annotations. Users need to contact the GeneCards team.)\n"
	    }
	    return
	}
    }
    
    ## Downloaded GH files have been unzipped and are available
    ###########################################################
    set g_AnnotSV(GHann) 1
    # Check if the user asked for these annotations in the configfile
    set test 0
    foreach col "GHid_elite GHid_not_elite GHtype GHgene_elite GHgene_not_elite GHtissue" {
	if {[lsearch -exact "$g_AnnotSV(outputColHeader)" $col] ne -1} {set test 1;break}
    }
    if {$test eq 0} {set g_AnnotSV(GHann) 0; return}


    puts "...GH configuration"
    
    ## Creation of the 2 GH formatted files (GRCh37 and GRCh38)
    ##   - # Header: chr GHstart GHend GHid_elite GHid_not_elite GHtype GHgene_elite GHgene_not_elite GHtissue
    ##   - 'date'_GH_GRCh37.sorted.bed 
    ##   - 'date'_GH_GRCh38.sorted.bed
    set GHfileFormatted37 "$GHdir/GRCh37/[clock format [clock seconds] -format "%Y%m%d"]_GH_GRCh37.sorted.bed"
    set GHfileFormatted38 "$GHdir/GRCh38/[clock format [clock seconds] -format "%Y%m%d"]_GH_GRCh38.sorted.bed"

    puts "\t...creation of [file tail $GHfileFormatted37] and [file tail $GHfileFormatted38]"
    puts "\t   (done only once during the first GH annotation)\n"

    ReplaceTextInFile "#chr\tGHstart\tGHend\tGHid_elite\tGHid_not_elite\tGHtype\tGHgene_elite\tGHgene_not_elite\tGHtissue" $GHfileFormatted37
    ReplaceTextInFile "#chr\tGHstart\tGHend\tGHid_elite\tGHid_not_elite\tGHtype\tGHgene_elite\tGHgene_not_elite\tGHtissue" $GHfileFormatted38

    # $GHelementsF
    # Header: chr element_start element_end GHid is_elite regulatory_element_type
    foreach L [LinesFromFile $GHelementsF] {
	set Ls [split $L "\t"]
	if {[regexp "^chr" $L]} {
	    set i_GHstart [lsearch -exact $Ls "element_start"]
	    set i_GHend [lsearch -exact $Ls "element_end"]
	    set i_GHid [lsearch -exact $Ls "GHid"]
	    set i_GH_is_elite [lsearch -exact $Ls "is_elite"]
	    set i_GHtype [lsearch -exact $Ls "regulatory_element_type"]
	    continue
	}
	regsub "chr" [lindex $Ls 0] "" chrom
	set GHid [lindex $Ls $i_GHid]
	lappend L_GHid($chrom) "$GHid"
	set L_GHgene_elite($GHid) ""
	set L_GHgene_not_elite($GHid) ""

	set GHstart [lindex $Ls $i_GHstart]
	set GHend [lindex $Ls $i_GHend]
	set coordGRCh38($GHid) "$chrom\t$GHstart\t$GHend"
	set is_elite($GHid) [lindex $Ls $i_GH_is_elite]
	set type($GHid) [lindex $Ls $i_GHtype]
    }

    # $GHassociationsF
    # Header: GHid symbol combined_score is_elite
    foreach L [LinesFromFile $GHassociationsF] {
	set Ls [split $L "\t"]
	if {[regexp "^GHid" $L]} {
	    set i_GHid [lsearch -exact $Ls "GHid"]
	    set i_regulatedGene [lsearch -exact $Ls "symbol"]
	    set i_GH_association_is_elite  [lsearch -exact $Ls "is_elite"]
	    continue
	}
	set GHid [lindex $Ls $i_GHid]
	set elite [lindex $Ls $i_GH_association_is_elite]
	set gene [lindex $Ls $i_regulatedGene] 
	if {$elite} {
	    lappend L_GHgene_elite($GHid) $gene
	} else {
	    lappend L_GHgene_not_elite($GHid) $gene
	}
    }

    # $GHtissuesF
    # Header: GHid source tissue category
    foreach L [LinesFromFile $GHtissuesF] {
	set Ls [split $L "\t"]
	if {[regexp "^GHid" $L]} {
	    set i_GHid [lsearch -exact $Ls "GHid"]
	    set i_GH_tissue [lsearch -exact $Ls "tissue"]
	    continue
	}
	set GHid [lindex $Ls $i_GHid]
	lappend L_GH_tissue($GHid) [lindex $Ls $i_GH_tissue]
    }

    # $GHhg19F
    # Header: chromosome start end GHid
    foreach L [LinesFromFile $GHhg19F] {
	set Ls [split $L "\t"]
	if {[regexp "^chromosome" $L]} {
	    set i_GHstart [lsearch -exact $Ls "start"]
	    set i_GHend [lsearch -exact $Ls "end"]
	    set i_GHid [lsearch -exact $Ls "GHid"]
	    continue
	}
	set GHid [lindex $Ls $i_GHid]
	regsub "chr" [lindex $Ls 0] "" chrom
	set GHstart [lindex $Ls $i_GHstart]
	set GHend [lindex $Ls $i_GHend]
	set coordGRCh37($GHid) "$chrom\t$GHstart\t$GHend"
    }

    # Writings:
    file delete -force $GHfileFormatted37.tmp
    file delete -force $GHfileFormatted38.tmp
    foreach chrom {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M MT} {
	if {![info exists L_GHid($chrom)]} {continue}
	set L_linesToWriteGRCh37 {}
	set L_linesToWriteGRCh38 {}
	foreach GHid $L_GHid($chrom) {
	    ## Header: chr GHstart GHend GHid_elite GHid_not_elite GHtype GHgene_elite GHgene_not_elite GHtissue
	    if {$is_elite($GHid)} {
		set gh "$GHid\t"
	    } else {
		set gh "\t$GHid"
	    }
	    catch {lappend L_linesToWriteGRCh37 "$coordGRCh37($GHid)\t$gh\t$type($GHid)\t[join $L_GHgene_elite($GHid) ";"]\t[join $L_GHgene_not_elite($GHid) ";"]\t[join $L_GH_tissue($GHid) ";"]"}
	    catch {lappend L_linesToWriteGRCh38 "$coordGRCh38($GHid)\t$gh\t$type($GHid)\t[join $L_GHgene_elite($GHid) ";"]\t[join $L_GHgene_not_elite($GHid) ";"]\t[join $L_GH_tissue($GHid) ";"]"}
	}
	set L_linesToWriteGRCh37 [lsort -command AscendingSortOnElement1 [lsort -command AscendingSortOnElement2 $L_linesToWriteGRCh37]]
	set L_linesToWriteGRCh38 [lsort -command AscendingSortOnElement1 [lsort -command AscendingSortOnElement2 $L_linesToWriteGRCh38]]
	WriteTextInFile "[join $L_linesToWriteGRCh37 "\n"]" $GHfileFormatted37.tmp
	WriteTextInFile "[join $L_linesToWriteGRCh38 "\n"]" $GHfileFormatted38.tmp
    }
    
    # Sorting of the bedfile:
    # Intersection with very large files can cause trouble with excessive memory usage.
    # A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm. 
    if {[catch {eval exec sort -k1,1 -k2,2n $GHfileFormatted37.tmp >> $GHfileFormatted37} Message]} {
	puts "-- checkGHfiles --"
	puts "sort -k1,1 -k2,2n $GHfileFormatted37.tmp >> $GHfileFormatted37"
	puts "$Message"
	    puts "Exit with error"
	exit 2
    }
    file delete -force $GHfileFormatted37.tmp
    if {[catch {eval exec sort -k1,1 -k2,2n $GHfileFormatted38.tmp >> $GHfileFormatted38} Message]} {
	puts "-- checkGHfiles --"
	puts "sort -k1,1 -k2,2n $GHfileFormatted38.tmp >> $GHfileFormatted38"
	puts "$Message"
	    puts "Exit with error"
	exit 2
    }
    file delete -force $GHfileFormatted38.tmp

    # Delete the downloaded files
    if {[file exists $GHfileFormatted37] && [file exists $GHfileFormatted38]} {
	foreach GHfile "$GHelementsF $GHassociationsF $GHhg19F $GHtissuesF $GHdir/ReadMe.txt" {
	    file delete -force $GHfile
	}
    } 
}


proc GHannotation {SVchrom SVstart SVend L_i} {

    global g_AnnotSV
    global GHtext
    global EliteGene
    global NotEliteGene
   

    set extannDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)"
    set GHfileFormatted [glob -nocomplain "$extannDir/FtIncludedInSV/GeneHancer/$g_AnnotSV(genomeBuild)/*_GH_$g_AnnotSV(genomeBuild).sorted.bed"] 
    
    if {![info exists GHtext(DONE)]} {
		
	# headerOutput "GHid_elite GHid_not_elite GHtype GHgene_elite GHgene_not_elite GHtissue"
	set L_GHtext(Empty) "{} {} {} {} {} {}"
	foreach i $L_i {
	    lappend GHtext(Empty) "[lindex $L_GHtext(Empty) $i]"
	}
	set GHtext(Empty) [join $GHtext(Empty) "\t"]

	# Intersect
	regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".intersect.GH" tmpFile
	set tmpFile "$g_AnnotSV(outputDir)/[file tail $tmpFile]"
	file delete -force $tmpFile	
	if {[catch {exec $g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(bedFile) -b $GHfileFormatted -wa -wb > $tmpFile} Message]} {
	    puts "-- GHannotation --"
	    puts "$g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(bedFile) -b $GHfileFormatted -wa -wb > $tmpFile"
	    puts "$Message"
	    puts "Exit with error"
	    exit 2
	}

	# Parse the intersection
	set f [open $tmpFile]
	while {![eof $f]} {
	    set L [gets $f]
	    if {$L eq ""} {continue}
	    set Ls [split $L "\t"]
	    set SVtoAnn_chrom [lindex $Ls 0]
	    set SVtoAnn_start [lindex $Ls 1]
	    set SVtoAnn_end   [lindex $Ls 2]

	    #set GHtissue         [lindex $Ls end]
	    #set GHgene_not_elite [lindex $Ls end-1]
	    #set GHgene_elite     [lindex $Ls end-2]
	    #set GHtype           [lindex $Ls end-3]
	    #set GHid_not_elite   [lindex $Ls end-4]
	    #set GHid_elite       [lindex $Ls end-5]
	    set GHend            [lindex $Ls end-6]
	    set GHstart          [lindex $Ls end-7]

	    # Select:
	    # - GH that share > XX % length with the SV to annotate
	    # - GH overlapping a SV that is an insertion
	    set GH_length [expr {$GHend-$GHstart}]
	    if {$GH_length == 0} {set GH_length 1}
	    set SVtoAnn_length [expr {$SVtoAnn_end - $SVtoAnn_start}]
	    if {$SVtoAnn_length<=1} {
		# SV to annotate is an insertion
		set SVtoAnn_length 1
	    } else {
		# SV to annotate is not an insertion
		if {$SVtoAnn_start < $GHstart} {
		    set overlap_start $GHstart
		} else {
		    set overlap_start $GHstart
		}
		if {$SVtoAnn_end < $GHend} {
		    set overlap_end $SVtoAnn_end
		} else {
		    set overlap_end $GHend
		}
		set overlap_length [expr {$overlap_end - $overlap_start}]
		
		# Keeping only GH with > 70% (default) overlap with the SV
		if {[expr {$overlap_length*100.0/$GH_length}] < $g_AnnotSV(overlap)} {continue}
		# No reciprocal overlap used for GH (it would not make sense)
	    }

	    # SV and GH overlap with > 70%,
	    # initialisation of the variables
	    set SVtoAnn "$SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end"
	    if {![info exists L_GHtissue($SVtoAnn)]} {
		set L_GHid_elite($SVtoAnn) ""
		set L_GHid_not_elite($SVtoAnn) ""
		#set L_GHtype($SVtoAnn) ""
		set L_GHgene_elite($SVtoAnn) ""
		set L_GHgene_not_elite($SVtoAnn) ""
		set L_GHtissue($SVtoAnn) ""
	    }
	    set GHtissue         [lindex $Ls end]
	    set GHgene_not_elite [lindex $Ls end-1]
	    set GHgene_elite     [lindex $Ls end-2]
	    set GHtype           [lindex $Ls end-3]
	    set GHid_not_elite   [lindex $Ls end-4]
	    set GHid_elite       [lindex $Ls end-5]

	    if {$GHtissue ne ""} {
		lappend L_GHtissue($SVtoAnn) {*}[split $GHtissue ";"]
	    }
	    if {$GHgene_not_elite ne ""} {
		lappend L_GHgene_not_elite($SVtoAnn) {*}[split $GHgene_not_elite ";"]
	    }
	    if {$GHgene_elite ne ""} {
		lappend L_GHgene_elite($SVtoAnn) {*}[split $GHgene_elite ";"]
	    }
	    lappend L_GHtype($SVtoAnn) {*}[split $GHtype "/"]   ; # = "Enhancer" or "Promoter" or "Promoter/Enhancer"
	    if {$GHid_not_elite ne ""} {
		lappend L_GHid_not_elite($SVtoAnn) $GHid_not_elite
	    }
	    if {$GHid_elite ne ""} {
		lappend L_GHid_elite($SVtoAnn) $GHid_elite
	    }
	}
	

	# Loading GH final annotation for each SV
	foreach SVtoAnn [array names L_GHtissue] {
	    # Very large SV (e.g. 30Mb) can overlap many GH (e.g. more than 2600). 
	    # e.g. 1:2806107-107058351
	    # So the “full” line can appear truncated in the output if visualized in a spreadsheet.  
	    # In order to avoid such embarrassing glitch, AnnotSV gives access to the list of the 20 first annotations.
	    set GHid_elite "[lsort -unique $L_GHid_elite($SVtoAnn)]"
	    set L_GHtext($SVtoAnn) ""
	    if {[llength $GHid_elite] > 20} {
		lappend L_GHtext($SVtoAnn) "[join [lrange $GHid_elite 0 19] ";"]..."
	    } else {
		lappend L_GHtext($SVtoAnn) "[join $GHid_elite ";"]"
	    }

	    set GHid_not_elite "[lsort -unique $L_GHid_not_elite($SVtoAnn)]"
	    if {[llength $GHid_not_elite] > 20} {
		lappend L_GHtext($SVtoAnn) "[join [lrange $GHid_not_elite 0 19] ";"]..."
	    } else {
		lappend L_GHtext($SVtoAnn) "[join $GHid_not_elite ";"]"
	    }

	    lappend L_GHtext($SVtoAnn) "[join [lsort -unique $L_GHtype($SVtoAnn)] ";"]"

	    set GHgene_elite "[lsort -unique $L_GHgene_elite($SVtoAnn)]"
	    if {[llength $GHgene_elite] > 20} {
		lappend L_GHtext($SVtoAnn) "[join [lrange $GHgene_elite 0 19] ";"]..."
		set EliteGene($SVtoAnn) $GHgene_elite  ;# Need of this information during the ranking
	    } else {
		lappend L_GHtext($SVtoAnn) "[join $GHgene_elite ";"]"
	    }

	    set GHgene_not_elite "[lsort -unique $L_GHgene_not_elite($SVtoAnn)]"
	    if {[llength $GHgene_not_elite] > 20} {
		lappend L_GHtext($SVtoAnn) "[join [lrange $GHgene_not_elite 0 19] ";"]..."
		set NotEliteGene($SVtoAnn) $GHgene_not_elite  ;# Need of this information during the ranking
	    } else {
		lappend L_GHtext($SVtoAnn) "[join $GHgene_not_elite ";"]"
	    }

	    set GHtissue "[lsort -unique $L_GHtissue($SVtoAnn)]"
	    if {[llength $GHtissue] > 20} {
		lappend L_GHtext($SVtoAnn) "[join [lrange $GHtissue 0 19] ";"]..."
	    } else {
		lappend L_GHtext($SVtoAnn) "[join $GHtissue ";"]"	
	    }

	    # Keep only the user requested columns (defined in the configfile)
	    set GHtext($SVtoAnn) ""
	    foreach i $L_i {
		lappend GHtext($SVtoAnn) "[lindex $L_GHtext($SVtoAnn) $i]"
	    }
	    set GHtext($SVtoAnn) [join $GHtext($SVtoAnn) "\t"]
	}
	catch {unset L_GHtext}
	set GHtext(DONE) 1	
	file delete -force $tmpFile
    }
    
    if {[info exist GHtext($SVchrom,$SVstart,$SVend)]} {
	return $GHtext($SVchrom,$SVstart,$SVend)
    } else {
	return $GHtext(Empty)
    }
}
