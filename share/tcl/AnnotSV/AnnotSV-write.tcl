############################################################################################################
# AnnotSV 2.5.2                                                                                            #
#                                                                                                          #
# AnnotSV: An integrated tool for Structural Variations annotation and ranking                             #
#                                                                                                          #
# Copyright (C) 2017-2020 Veronique Geoffroy (veronique.geoffroy@inserm.fr)                                #
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                                             #
# GNU General Public License for more details.                                                             #
#                                                                                                          #
# You should have received a copy of the GNU General Public License                                        #
# along with this program; If not, see <http://www.gnu.org/licenses/>.                                     #
############################################################################################################



proc OrganizeAnnotation {} {

    global g_AnnotSV
    global g_Lgenes
    global VCFheader
    global g_numberOfAnnotationCol
    global headerFileToRemove 
    global L_Candidates
    global g_SVLEN
    global g_ExtAnnotation
    global g_rankingExplanations
    global g_re
    global g_HITS
    
    # OUTPUT
    ###############
    set FullAndSplitBedFile "$g_AnnotSV(outputDir)/$g_AnnotSV(outputFile).tmp" ;# created in AnnotSV-genes.tcl
    set outputFile "$g_AnnotSV(outputDir)/$g_AnnotSV(outputFile)" 


    #########################################################################################
    ################### Writing of the header (first line of the output) ####################
    #########################################################################################
    set headerOutput "AnnotSV ID\tSV chrom\tSV start\tSV end\tSV length"
    if {[info exist VCFheader]} {
	# SVinputFile = VCF
	append headerOutput "\t$VCFheader"
	set theBEDlength [expr {[llength [split $headerOutput "\t"]]-2}] ; # we remove "2" for the 2 columns not in the input file ("AnnotSV ID" and "SV length")
	if {$g_AnnotSV(SVinputInfo)} {
	    set g_AnnotSV(formatColNumber) 10 ;# used in AnnotSV-filteredVCF.tcl
	} else {
	    set g_AnnotSV(formatColNumber) 6
	}
	if {[lsearch -exact $g_AnnotSV(outputColHeader) "Samples_ID"] ne "-1"} {
	    incr g_AnnotSV(formatColNumber) +1
	}
	set g_AnnotSV(SVinputInfo) 1 ; # The created bedfile contains only the info to report
    } else {
	# SVinputFile = BED
	if {$g_AnnotSV(SVinputInfo)} {
	    # The user wants to keep all the columns from the SV BED input file
	    regsub -nocase ".bed$" $g_AnnotSV(SVinputFile) ".header.tsv" headerSVinputFile
	    set theBEDlength [llength [split [FirstLineFromFile $g_AnnotSV(bedFile)] "\t"]] ;# theBEDlength = number of col in the SV input BED file

	    if {[file exists $headerSVinputFile]} {
		# The user has given a header for the SV BED input file
		set headerFromTheUser [split [FirstLineFromFile $headerSVinputFile] "\t"]	
		if {[llength $headerFromTheUser] ne $theBEDlength} { ;# A stupid error can be to have a "\t" at the end of the header line
		    puts "Numbers of columns from $g_AnnotSV(SVinputFile) and $headerSVinputFile are different ($theBEDlength != [llength $headerFromTheUser])" 
		    puts "=> Can not report: $headerFromTheUser"
		} else {
		    if {$g_AnnotSV(svtBEDcol) ne -1} {
			set headerFromTheUser [lreplace $headerFromTheUser $g_AnnotSV(svtBEDcol) $g_AnnotSV(svtBEDcol) "SV type"]
		    }
		    if {$g_AnnotSV(samplesidBEDcol) ne -1} {
			set headerFromTheUser [lreplace $headerFromTheUser $g_AnnotSV(samplesidBEDcol) $g_AnnotSV(samplesidBEDcol) "Samples_ID"]
		    } 
		    append headerOutput "\t[join [lrange $headerFromTheUser 3 end] "\t"]"		
		}	
	    } else {
		# No header given by the user
		set i 5 ; # headerOutput = "AnnotSV ID   SV chrom    SV start	   SV end	SV length"
		set j [expr {$theBEDlength+2}]
		while {$i < $j} {
		    if {$i eq $g_AnnotSV(svtTSVcol)} {
			append headerOutput "\tSV type"
		    } elseif {$i eq $g_AnnotSV(samplesidTSVcol)} {
			append headerOutput "\tSamples_ID"
		    } else {
			append headerOutput "\t"
		    }
		    incr i
		}
	    }
	} else {
	    # At least the "SV type" and the "samples_ID" columns should be reported for the ranking
	    if {$g_AnnotSV(svtTSVcol) ne -1} {
		# The user doesn't want to keep all the columns from the SV BED input file. We keep only the "SV type" (for the ranking) 
		append headerOutput "\tSV type"
	    }
	    if {$g_AnnotSV(samplesidTSVcol) ne -1} {
		append headerOutput "\tSamples_ID"
	    }
	}
    }
    append headerOutput "\tAnnotSV type\tGene name\tNumber of genes\ttx\ttx start\ttx end\toverlapped tx length\toverlapped CDS length\tframeshift\tNumber of exons\tlocation\tlocation2\tdistNearestSS\tnearestSStype\tintersectStart\tintersectEnd\tRE_gene_exomiser"
    	      

    ### Search for "ref" and "alt" information (to define the AnnotSV_ID)
    set i_ref [lsearch -exact [split $headerOutput "\t"] "REF"]
    set i_alt [lsearch -exact [split $headerOutput "\t"] "ALT"]
    set i_ref [expr {$i_ref-2}] ;# (-2) because of the insertion of the END and SVTYPE values.
    set i_alt [expr {$i_alt-2}] ;# (-2) because of the insertion of the END and SVTYPE values.

    ####### "Benign SV header"
    foreach svtype "gain loss ins inv" {
	if {[lsearch -regexp "$g_AnnotSV(outputColHeader)" "^B_[string tolower ${svtype}]_"] eq -1} { continue }
	append headerOutput "\tB_${svtype}_source\tB_${svtype}_coord"
    }

    set usersDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Users/$g_AnnotSV(genomeBuild)"
    ####### usersDir: "SVincludedInFt header"
    if {[glob -nocomplain $usersDir/SVincludedInFt/*.formatted.sorted.bed] ne ""} {
	foreach formattedUserBEDfile [glob -nocomplain $usersDir/SVincludedInFt/*.formatted.sorted.bed] {
	    regsub -nocase ".formatted.sorted.bed$" $formattedUserBEDfile ".header.tsv" userHeaderFile
	    set L1header [split [FirstLineFromFile $userHeaderFile] "\t"]
	    set L_headerColName [lrange $L1header 3 end]
	    lappend SVincludedInFtHeader "[join $L_headerColName "\t"]"
	}
	append headerOutput "\t[join $SVincludedInFtHeader "\t"]"
    }

    ####### "Pathogenic SV header"
    foreach svtype "gain loss ins inv" {
	if {[lsearch -regexp "$g_AnnotSV(outputColHeader)" "^P_[string tolower ${svtype}]_"] eq -1} { continue }
	append headerOutput "\tP_${svtype}_phen\tP_${svtype}_hpo\tP_${svtype}_source\tP_${svtype}_coord"
    }

    ####### "TAD header"
    if {$g_AnnotSV(tadAnn)} {
	set g_AnnotSV(tadAnn_i) ""
	set j 0
	foreach col "TADcoordinates ENCODEexperiments" {
	    if {[lsearch -exact "$g_AnnotSV(outputColHeader)" $col] ne -1} {
		append headerOutput "\t$col"
		lappend g_AnnotSV(tadAnn_i) $j
	    }
	    incr j
	}
	if {$g_AnnotSV(tadAnn_i) eq ""} {set g_AnnotSV(tadAnn) 0}
    }

    ####### "VCF files header"
    if {$g_AnnotSV(snvIndelFiles) ne ""} {
	foreach sample $g_AnnotSV(snvIndelSamples) {
	    append headerOutput "\t#hom($sample)\t#htz($sample)\t#htz/allHom($sample)\t#htz/total(cohort)"
	}
	append headerOutput "\t#total(cohort)"
    }
    if {$g_AnnotSV(candidateSnvIndelFiles) ne ""} {
	foreach sample $g_AnnotSV(candidateSnvIndelSamples) {
	    append headerOutput "\tcompound-htz($sample)"
	}
    }

    ####### "FtIncludedInSV header"
    if {[glob -nocomplain $usersDir/FtIncludedInSV/*.formatted.sorted.bed] ne ""} {
	foreach formattedUserBEDfile [glob -nocomplain $usersDir/FtIncludedInSV/*.formatted.sorted.bed] {
	    regsub -nocase ".formatted.sorted.bed$" $formattedUserBEDfile ".header.tsv" userHeaderFile
	    set L1header   [split [FirstLineFromFile $userHeaderFile] "\t"]
	    set L_headerColName [lrange $L1header 3 end]
	    lappend FtIncludedInSVHeader "[join $L_headerColName "\t"]"
	}
	append headerOutput "\t[join $FtIncludedInSVHeader "\t"]"
    }

    ####### "Breakpoints header"
    if {$g_AnnotSV(gcContentAnn)} {
	append headerOutput "\tGCcontent_left\tGCcontent_right"
    }
    if {$g_AnnotSV(repeatAnn)} {
	append headerOutput "\tRepeats_coord_left\tRepeats_type_left\tRepeats_coord_right\tRepeats_type_right"
    }
    if {$g_AnnotSV(gapAnn)} {
	append headerOutput "\tGap_left\tGap_right"
    }
    if {$g_AnnotSV(segdupAnn)} {
	append headerOutput "\tSegDup_left\tSegDup_right"
    }
    if {$g_AnnotSV(ENCODEblacklistAnn)} {
	append headerOutput "\tENCODE_blacklist_left\tENCODE_blacklist_characteristics_left\tENCODE_blacklist_right\tENCODE_blacklist_characteristics_right"
    }

    ####### "Genes-based header"
    if {$g_AnnotSV(genesBasedAnn)} {
	set L_allGenesBasedHeader ""
	if {[info exists g_AnnotSV(extann)] && $g_AnnotSV(extann) ne ""} { ; # g_AnnotSV(extann) defined in the main
	    ExternalAnnotations
	    foreach F [ExternalAnnotations L_Files] {
		# we remove the first "genes" column (not reported in output, used as an ID)
		regsub "\[^\t\]+\t" [ExternalAnnotations $F Header] "" extannHeader
		lappend L_allGenesBasedHeader {*}[split $extannHeader "\t"]
	    }
	}
	set g_AnnotSV(genesBasedAnn_i) ""
	set j 0
	foreach col "$L_allGenesBasedHeader" {
	    if {[lsearch -exact "$g_AnnotSV(outputColHeader)" $col] ne -1} {
		append headerOutput "\t$col"
		lappend g_AnnotSV(genesBasedAnn_i) $j
	    }
	    incr j
	}
	if {$g_AnnotSV(genesBasedAnn_i) eq ""} {set g_AnnotSV(genesBasedAnn) 0}
    }
 
    # Preparation for the ranking (from benign to pathogenic)
    SVprepareRanking $headerOutput    ; # svtTSVcol (for VCF input file) is defined there

    ####### "Exomiser header"
    if {$g_AnnotSV(hpo) ne ""} {
	append headerOutput "\tEXOMISER_GENE_PHENO_SCORE\tHUMAN_PHENO_EVIDENCE\tMOUSE_PHENO_EVIDENCE\tFISH_PHENO_EVIDENCE"
    }
    
    ####### "Ranking header"
    if {$g_AnnotSV(svtTSVcol) eq -1 && $g_AnnotSV(organism) eq "Human"} { ; # SV type is required for the ranking of human SV
	puts "\nWARNING: AnnotSV requires the SV type (duplication, deletion...) to classify the SV"
	puts "Not provided (svtBEDcol = -1)"
	puts "=> No SV ranking"
	set g_AnnotSV(ranking) 0
    }  
    if {$g_AnnotSV(ranking)} {
	if {[lsearch -exact "$g_AnnotSV(outputColHeader)" "ranking decision criteria"] ne -1} {
	    append headerOutput "\tAnnotSV ranking\tranking decision criteria"
	} else {
	    append headerOutput "\tAnnotSV ranking"
	}
    }
    
    
    ReplaceTextInFile $headerOutput $outputFile



    ##################################################################################
    ################### Display of the annotations to realize ########################
    ##################################################################################
    puts "...listing of the annotations to realized ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])" 
    puts "\t...Genes annotation" 
    puts "\t(with $g_AnnotSV(genesFile))"
    
    ####### "Regulatory elements annotations"
    puts "\t...Regulatory elements annotations"
    if {$g_AnnotSV(promAnn)} {puts "\t\t...Promoter annotations"}
    if {$g_AnnotSV(EAann)} {puts "\t\t...EnhancerAtlas annotations"}
    if {$g_AnnotSV(GHAnn)} {puts "\t\t...GeneHancer annotations"}
    
    #######  Annotations with benign genes or genomic regions (SVincludedInFt)
    puts "\t...Annotations with benign genes or genomic regions"
    puts "\t\t...gnomAD annotation"
    puts "\t\t...ClinVar annotation"
    puts "\t\t...ClinGen annotation"
    puts "\t\t...DGV annotation"
    puts "\t\t...DDD annotation"
    puts "\t\t...1000g annotation"
    puts "\t\t...Ira M. Hall's lab annotation"
    # "SVincludedInFt"
    foreach formattedUserBEDfile [glob -nocomplain $usersDir/SVincludedInFt/*.formatted.sorted.bed] {
	puts "\t\t...[file tail $formattedUserBEDfile]"
	regsub -nocase ".formatted.sorted.bed$" $formattedUserBEDfile ".header.tsv" userHeaderFile
	set L1header [split [FirstLineFromFile $userHeaderFile] "\t"]
	set L_headerColName [lrange $L1header 3 end]
	# Number of columns of annotation in the user bedfile (without the 3 columns "chrom start end")
	set nColHeader [llength $L_headerColName]
	puts "\t\t($nColHeader annotations columns: [join $L_headerColName ", "])"
    }
    #######  Annotations with pathogenic genes or genomic regions (FtIncludedInSV)
    puts "\t...Annotations with pathogenic genes or genomic regions"
    puts "\t\t...dbVar annotation"
    puts "\t\t...ClinVar annotation"
    puts "\t\t...ClinGen annotation"
    ####### "FtIncludedInSV"
    puts "\t...Annotations with features overlapped with the SV ($g_AnnotSV(overlap) %)"
    if {$g_AnnotSV(tadAnn)} {puts "\t\t...TAD annotation"}
    foreach formattedUserBEDfile [glob -nocomplain $usersDir/FtIncludedInSV/*.formatted.sorted.bed] {
	puts "\t\t...[file tail $formattedUserBEDfile]"
	regsub -nocase ".formatted.sorted.bed$" $formattedUserBEDfile ".header.tsv" userHeaderFile
	set L1header   [split [FirstLineFromFile $userHeaderFile] "\t"]
	set L_headerColName [lrange $L1header 3 end]
	append FtIncludedInSVHeader "\t[join $L_headerColName "\t"]"
	# Number of columns of annotation in the user bedfile (without the 3 columns "chrom start end")
	set nColHeader [llength $L_headerColName]
	puts "\t\t($nColHeader annotations columns: [join $L_headerColName ", "])"
    }
    ####### "Breakpoints annotations"
    puts "\t...Breakpoints annotations"
    if {$g_AnnotSV(gcContentAnn)} {puts "\t\t...GC content annotation"}
    if {$g_AnnotSV(repeatAnn)} {puts "\t\t...Repeat annotation"}
    if {$g_AnnotSV(gapAnn)} {puts "\t\t...Gap annotation"}
    if {$g_AnnotSV(segdupAnn)} {puts "\t\t...Segmental duplication annotation"}
    if {$g_AnnotSV(ENCODEblacklistAnn)} {puts "\t\t...ENCODE blacklist annotation"}

    ####### "Genes-based annotations"
    puts "\t...Genes-based annotations"
    if {$g_AnnotSV(genesBasedAnn)} {
	puts "[join $g_ExtAnnotation(display) "\n"]"
    }
    ####### "Exomiser annotation"
    if {$g_AnnotSV(hpo) ne ""} {puts "\t...Exomiser annotations"}


    

    ########################################################################
    ################### Display: annotation in progress ####################
    ########################################################################
    puts "\n...annotation in progress ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])" 


    
    #######################################################################################################
    ##   Prepare the intersection between "SV / genes" and "annotation BED files" (benign SV...):
    ##
    ##   => creation of $g_AnnotSV(fullAndSplitBedFile),
    ##      a bedfile with the coordinates (only 3 columns: chrom, start and end) of:
    ##      - the SV to annotate (full)
    ##      - the intersections with genes (split)
    #######################################################################################################
    set g_AnnotSV(fullAndSplitBedFile) "$g_AnnotSV(outputDir)/[file tail $g_AnnotSV(bedFile)].users.bed"
    file delete -force $g_AnnotSV(fullAndSplitBedFile)
    set L_UsersText {}
    set f [open "$FullAndSplitBedFile"]
    while {! [eof $f]} {
        set L [gets $f]
	if {$L eq ""} {continue}
	set Ls [split $L "\t"]
	set AnnotSVtype [lindex $Ls end]
	if {$AnnotSVtype eq "split"} {
	    set SVchrom [lindex $Ls 0]
	    set SVleft  [lindex $Ls 1]
	    set SVright [lindex $Ls 2]
	    set exonStarts [lindex $Ls end-4]
	    set exonEnds   [lindex $Ls end-3]
	    set tx_left [lindex [split $exonStarts ","] 0]
	    set tx_right [lindex [split $exonEnds ","] end-1]
	    if {$SVleft<$tx_left} {set intersectStart "$tx_left"} else {set intersectStart "$SVleft"}
	    if {$SVright<$tx_right} {set intersectEnd "$SVright"} else {set intersectEnd "$tx_right"}
	    lappend L_UsersText "$SVchrom\t$intersectStart\t$intersectEnd"
	} else {
	    lappend L_UsersText "[join [lrange $Ls 0 2] "\t"]"
	}
    }
    close $f
    set L_UsersText [lsort -unique $L_UsersText]
    WriteTextInFile [join $L_UsersText "\n"] $g_AnnotSV(fullAndSplitBedFile)
    checkBed $g_AnnotSV(fullAndSplitBedFile)
    regsub -nocase ".bed$" $g_AnnotSV(fullAndSplitBedFile) ".formatted.bed" newFullAndSplitBedFile

    # Sorting of the bedfile:
    # Intersection with very large files can cause trouble with excessive memory usage.
    # A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm. 
    set g_AnnotSV(fullAndSplitBedFile) "$g_AnnotSV(outputDir)/[file tail $g_AnnotSV(bedFile)].users.sorted.bed"
    set sortTmpFile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_sort.tmp.bash"
    ReplaceTextInFile "#!/bin/bash" $sortTmpFile
    WriteTextInFile "# The locale specified by the environment can affects the traditional sort order. We need to use native byte values." $sortTmpFile
    WriteTextInFile "export LC_ALL=C" $sortTmpFile
    WriteTextInFile "sort -k1,1 -k2,2n $newFullAndSplitBedFile > $g_AnnotSV(fullAndSplitBedFile)" $sortTmpFile
    file attributes $sortTmpFile -permissions 0755
    if {[catch {eval exec bash $sortTmpFile} Message]} {
	puts "-- OrganizeAnnotation --"
	puts "sort -k1,1 -k2,2n $newFullAndSplitBedFile > $g_AnnotSV(fullAndSplitBedFile)"
	puts "$Message"
	puts "Exit with error"
	exit 2
    }
    file delete -force $sortTmpFile 
    file delete -force $newFullAndSplitBedFile
    unset L_UsersText

    
    
    ########################################################################
    ################### Parse the "FullAndSplitBedFile" ####################
    ########################################################################

    set L_TextToWrite {}
    set f [open "$FullAndSplitBedFile"]
    while {! [eof $f]} {
        set L [gets $f]
	if {$L eq ""} {continue}
	set Ls [split $L "\t"]

	# Full + split
	set SVchrom   [lindex $Ls 0]
	set SVleft    [lindex $Ls 1]
	set SVright   [lindex $Ls 2]
	set AnnotSVtype   [lindex $Ls end]                 ;# full or split
	if {$g_AnnotSV(svtBEDcol) ne -1} { 
	    set SVtype [lindex $Ls "$g_AnnotSV(svtBEDcol)"]                     ;# DEL, DUP, <CN0>...
	} else {
	    set SVtype ""
	}
	if {$g_AnnotSV(samplesidBEDcol) ne -1} { 
	    set Samplesid [lindex $Ls "$g_AnnotSV(samplesidBEDcol)"]                    
	} else {
	    set Samplesid ""
	}

	# For the use of the -typeOfAnnotation option: keep only the corresponding lines (full or split or both)
	if {$g_AnnotSV(typeOfAnnotation) eq "full" && $AnnotSVtype eq "split"} {continue}
	if {$g_AnnotSV(typeOfAnnotation) eq "split" && $AnnotSVtype eq "full"} {continue}


	if {$AnnotSVtype eq "split"} {
	    # split
	    set txStart    [lindex $Ls end-11]
	    set txEnd      [lindex $Ls end-10]
	    set strand     [lindex $Ls end-9]
	    set geneName   [lindex $Ls end-8]
	    set NbGenes    ""
	    set transcript [lindex $Ls end-7]
	    set CDSstart   [lindex $Ls end-6]
	    set CDSend     [lindex $Ls end-5]
	    set exonStarts [lindex $Ls end-4]
	    set exonEnds   [lindex $Ls end-3]
	    set CDSl       [lindex $Ls end-2]
	    set txL        [lindex $Ls end-1]
	    set locationStart ""
	    set locationEnd ""		
	    set distNearestSSright ""
	    set distNearestSSrightA ""
	    set distNearestSSrightB ""
	    set distNearestSSleft ""
	    set distNearestSSleftA ""
	    set distNearestSSleftB ""
	    set distNearestSS ""
	    set nearestSStypeRight ""
	    set nearestSStypeLeft ""
	    set nearestSStype ""
	    set regionStart ""
	    set regionEnd ""
	    set intersectStart ""
	    set intersectEnd ""
	    set nbExons [expr {[llength [split $exonStarts ","]]-1}]
    	    if {[string is integer $CDSl]} {
    	        if {[expr {$CDSl%3}] eq 0} {
                    set frameshift "no"
                } else {
            	    set frameshift "yes"
                } 
            } else {set frameshift ""}
	} else {
	    set SV "[join [lrange $Ls 0 2] "\t"]"
	    # full
	    if {[info exists g_Lgenes($SV)]} {
		set geneName "$g_Lgenes($SV)"  ; # SV/oldSV defined with: "chrom, start, end, SVtype":	    set transcript         ""
		set NbGenes "[llength [split $g_Lgenes($SV) "/"]]"
	    } else {
		set geneName ""
		set NbGenes 0
	    }
	    set transcript ""
	    set CDSl       ""
	    set txStart    ""
	    set txEnd      ""
	    set txL        ""
	    set location   ""
	    set location2  ""
	    set distNearestSS ""
	    set nearestSStype ""
	    set intersect  "\t"
	    set nbExons    ""
	    set frameshift ""
	}	    

	# To select only the SV "split" annotations overlapping a gene from the "candidateGenesFile"
	if {$g_AnnotSV(candidateGenesFiltering) eq "yes"} {
	    set test 0
	    if {$geneName eq ""} {
		set test 1
	    } else {
		foreach g [split $geneName "/"] {
		    if {[lsearch -exact $L_Candidates $g] eq -1} {
			set test 1
		    }
		}	
	    }
	    if {$test} {continue}
	}


	# Definition of "locationStart" and "locationEnd" variables
	# Definition of "distNearestSS" and "nearestSStype" variables
	# NearestSStype:
	#   5' = splice donor site (on the right of the exon if strand +; else on the left) 
	#   3' = splice acceptor site (on the left of the exon if strand +; else on the right)
	if {$AnnotSVtype eq "split"} {
	    set nbEx [expr {[llength [split $exonStarts ","]]-1}] ; # Example: "1652370,1657120,1660664,1661968," --> 1652370 1657120 1660664 1661968 {}
	    
	    set tx_left [lindex [split $exonStarts ","] 0]
	    set tx_right [lindex [split $exonEnds ","] end-1]

	    if {$SVleft<=$tx_left} {         ; # SV begins before tx start (or tx end for strand "-")
		if {$strand eq "+"} {
		    set locationStart "txStart"
		} else {
		    set locationEnd "txEnd"
		}
	    }

	    set i 0
	    set previousB "[lindex [split $exonEnds ","] 0]"
	    foreach A [split $exonStarts ","] B [split $exonEnds ","] {
		if {$A eq "" || $B eq ""} {continue}
		incr i
		# SV left
		if {$SVleft<$B} {
		    if {$SVleft>$A} {
			# in an exon
			set distNearestSSleftA "[expr {$SVleft-$A}]"
			set distNearestSSleftB "[expr {$B-$SVleft}]"
			if {$i eq 1} { ;# Only 1 splice site on the first exon 
			    set distNearestSSleft $distNearestSSleftB
			    if {$strand eq "+"} {set nearestSStypeLeft "5'"} else {set nearestSStypeLeft "3'"}
			} elseif {$i eq $nbEx} { ;# Only 1 splice site on the last exon 
			    set distNearestSSleft $distNearestSSleftA
			    if {$strand eq "+"} {set nearestSStypeLeft "3'"} else {set nearestSStypeLeft "5'"}
			} elseif {$distNearestSSleftA < $distNearestSSleftB} {
			    set distNearestSSleft $distNearestSSleftA
			    if {$strand eq "+"} {set nearestSStypeLeft "3'"} else {set nearestSStypeLeft "5'"}
			} else {
			    set distNearestSSleft $distNearestSSleftB
			    if {$strand eq "+"} {set nearestSStypeLeft "5'"} else {set nearestSStypeLeft "3'"}
			}
			if {$strand eq "+"} {
			    set locationStart "exon$i"
			} else {
			    set locationEnd "exon[expr {$nbEx-$i+1}]" ; # gene on the strand "-"
			}
		    } elseif {$SVleft>=$previousB} {
			# in an intron
			set distNearestSSleftB "[expr {$SVleft-$previousB}]"
			set distNearestSSleftA "[expr {$A-$SVleft}]"
			if {$distNearestSSleftB < $distNearestSSleftA} {
			    set distNearestSSleft $distNearestSSleftB
			    if {$strand eq "+"} {set nearestSStypeLeft "5'"} else {set nearestSStypeLeft "3'"}
			} else {
			    set distNearestSSleft $distNearestSSleftA
			    if {$strand eq "+"} {set nearestSStypeLeft "3'"} else {set nearestSStypeLeft "5'"}
			}
			if {$strand eq "+"} {
			    set locationStart "intron[expr {$i-1}]"
			} else {
			    set locationEnd "intron[expr {$nbEx-$i+1}]"
			}
		    }
		}
		# SV right	
		if {$SVright<$B} {
		    if {$SVright>$A} {
			# in an exon
			set distNearestSSrightA "[expr {$SVright-$A}]"
			set distNearestSSrightB "[expr {$B-$SVright}]"
			if {$i eq 1} { ;# Only 1 splice site on the first exon 
			    set distNearestSSright $distNearestSSrightB
			    if {$strand eq "+"} {set nearestSStypeRight "5'"} else {set nearestSStypeRight "3'"}
			} elseif {$i eq $nbEx} { ;# Only 1 splice site on the last exon 
			    set distNearestSSright $distNearestSSrightA
			    if {$strand eq "+"} {set nearestSStypeRight "3'"} else {set nearestSStypeRight "5'"}
			} elseif {$distNearestSSrightA < $distNearestSSrightB} {
			    set distNearestSSright $distNearestSSrightA
			    if {$strand eq "+"} {set nearestSStypeRight "3'"} else {set nearestSStypeRight "5'"}
			} else {
			    set distNearestSSright $distNearestSSrightB
			    if {$strand eq "+"} {set nearestSStypeRight "5'"} else {set nearestSStypeRight "3'"}
			}
			if {$strand eq "+"} {
			    set locationEnd "exon$i"
			} else {
			    set locationStart "exon[expr {$nbEx-$i+1}]"
			}
			break
		    } elseif {$SVright>=$previousB} {
			# in an intron
			set distNearestSSrightB "[expr {$SVright-$previousB}]"
			set distNearestSSrightA "[expr {$A-$SVright}]"
			if {$distNearestSSrightB < $distNearestSSrightA} {
			    set distNearestSSright $distNearestSSrightB
			    if {$strand eq "+"} {set nearestSStypeRight "5'"} else {set nearestSStypeRight "3'"}
			} else {
			    set distNearestSSright $distNearestSSrightA
			    if {$strand eq "+"} {set nearestSStypeRight "3'"} else {set nearestSStypeRight "5'"}
			}
			if {$strand eq "+"} {
			    set locationEnd "intron[expr {$i-1}]" 
			} else {
			    set locationStart "intron[expr {$nbEx-$i+1}]"
			}
			break
		    }
		}
		set previousB "$B"
	    }	    
	    if {$locationEnd eq "" || $locationStart eq ""} {     ; # SV finishes after tx end
		if {$strand eq "+"} {
		    set locationEnd "txEnd"
		} else {
		    set locationStart "txStart"
		}
	    }

	    set location "${locationStart}-${locationEnd}"
	    if {$location ne "txStart-txEnd"} {
		if {$distNearestSSright eq ""} {
		    if {$distNearestSSleft ne ""} {
			set distNearestSS "$distNearestSSleft"
			set nearestSStype "$nearestSStypeLeft"
		    }
		} elseif {$distNearestSSleft eq ""} {
		    set distNearestSS "$distNearestSSright"
		    set nearestSStype "$nearestSStypeRight"
		} elseif {$distNearestSSright < $distNearestSSleft} {  
		    set distNearestSS "$distNearestSSright"
		    set nearestSStype "$nearestSStypeRight"
		} else {
		    set distNearestSS "$distNearestSSleft"
		    set nearestSStype "$nearestSStypeLeft"
		}
	    }
	    
	    # Definition of "regionStart" and "regionEnd" variables
	    # Warning: some genes have CDSstart=CDSend=txEnd (<=> no CDS)
	    if {$CDSstart eq $CDSend} {
		set regionLeft "UTR"
		set regionRight "UTR"
	    } else {
		# SVleft
		if {$SVleft < $CDSstart} {     
		    if {$strand eq "+"} {
			set regionLeft "5'UTR"
		    } else {
			set regionLeft "3'UTR"
		    }  
		} elseif {$SVleft < $CDSend} {
		    set regionLeft "CDS"
		} else {
		    if {$strand eq "+"} {
			set regionLeft "3'UTR"
		    } else {
			set regionLeft "5'UTR"
		    }  
		}
		# SVright
		if {$SVright < $CDSstart} {     
		    if {$strand eq "+"} {
			set regionRight "5'UTR"
		    } else {
			set regionRight "3'UTR"
		    }  
		} elseif {$SVright < $CDSend} {
		    set regionRight "CDS"
		} else {
		    if {$strand eq "+"} {
			set regionRight "3'UTR"
		    } else {
			set regionRight "5'UTR"
		    }  
		}
	    }
	    if {$strand eq "+"} {
		set regionStart $regionLeft
		set regionEnd   $regionRight
	    } else {
		set regionStart $regionRight
		set regionEnd   $regionLeft	
	    }
	    set location2 "$regionStart-$regionEnd"
	    if {[regexp "(\[^-\]+)-(\[^-\]+)" $location2 match titi tutu]} {
		if {$titi eq $tutu} {set location2 $titi}
	    }
	} 

	# Definition of "intersectStart" and "intersectEnd" variables
	if {$AnnotSVtype eq "split"} {
	    if {$SVleft<$tx_left} {set intersectStart "$tx_left"} else {set intersectStart "$SVleft"}
	    if {$SVright<$tx_right} {set intersectEnd "$SVright"} else {set intersectEnd "$tx_right"}
	    set intersect "$intersectStart\t$intersectEnd"
	} 

	# Regulatory elements annotation (only for the full lines)
	set reText ""
	if {$AnnotSVtype eq "full"} {
	    if {[info exists g_re($SVchrom\t$SVleft\t$SVright)]} {
		set L_regulatedGenes $g_re($SVchrom\t$SVleft\t$SVright)
		foreach gName "$L_regulatedGenes" {
		    set HITS ""
		    catch {set HITS "$g_HITS($gName)"}
		    if {$g_AnnotSV(hpo) ne ""} {
			set exomiserScore "EX=[ExomiserAnnotation $gName "score"]"
		    } else {set exomiserScore ""}

		    set lAnn {}
		    if {$HITS ne ""} {lappend lAnn $HITS}
		    if {$exomiserScore ne ""} {lappend lAnn $exomiserScore}
		    if {$lAnn ne ""} {
			lappend reText "$gName ([join $lAnn "; "])"
		    } else {
			lappend reText "$gName"
		    }
		}
		set reText [join $reText "; "]
	    } 
	} 
	
	# Annotations with benign genes or genomic regions (SVincludedInFt)
	if {$AnnotSVtype eq "split"} {
	    set benignText "[benignSVannotation $SVchrom $intersectStart $intersectEnd]"
	} else {
	    set benignText "[benignSVannotation $SVchrom $SVleft $SVright]"
	} 

	# User SVincludedInFt BED annotations. 
	set L_SVincludedInFtText {}
   	foreach formattedUserBEDfile [glob -nocomplain $usersDir/SVincludedInFt/*.formatted.sorted.bed] {
	    if {$AnnotSVtype eq "split"} {
		lappend L_SVincludedInFtText "[userBEDannotation $formattedUserBEDfile $SVchrom $intersectStart $intersectEnd]"
	    } else {
		lappend L_SVincludedInFtText "[userBEDannotation $formattedUserBEDfile $SVchrom $SVleft $SVright]"
	    } 
	}
	set SVincludedInFTtext [join $L_SVincludedInFtText "\t"]

	# Annotations with pathogenic genes or genomic regions (FtIncludedInSV)
	if {$AnnotSVtype eq "split"} {
	    set pathogenicText "[pathogenicSVannotation $SVchrom $intersectStart $intersectEnd]"
	} else {
	    set pathogenicText "[pathogenicSVannotation $SVchrom $SVleft $SVright]"
	} 

	# User FtIncludedInSV BED annotations. 
	set L_FtIncludedInSVtext {}
   	foreach formattedUserBEDfile [glob -nocomplain $usersDir/FtIncludedInSV/*.formatted.sorted.bed] {
	    if {$AnnotSVtype eq "split"} {
		lappend L_FtIncludedInSVtext "[userBEDannotation $formattedUserBEDfile $SVchrom $intersectStart $intersectEnd]"
	    } else {
		lappend L_FtIncludedInSVtext "[userBEDannotation $formattedUserBEDfile $SVchrom $SVleft $SVright]"
	    } 
	}
	set FtIncludedInSVtext [join $L_FtIncludedInSVtext "\t"]

	# Genes-based annotations.
	if {$g_AnnotSV(genesBasedAnn)} {
	    #   -> Number of columns from each Genes-based file
	    if {[info exists g_AnnotSV(extann)] && $g_AnnotSV(extann) ne ""} {		    
		foreach F [ExternalAnnotations L_Files] {
		    set nbColumns($F) [expr {[llength [split [ExternalAnnotations $F Header] "\t"]]-1}]
		}
	    }
	    #   -> Annotations
	    set L_genesBasedText {}
	    if {$AnnotSVtype eq "split"} {
		######### split lines ###################################
		if {[info exists g_AnnotSV(extann)] && $g_AnnotSV(extann) != ""} {		    
		    foreach F [ExternalAnnotations L_Files] {
			set AnnotFound "[ExternalAnnotations $F $geneName]"		
			if {$AnnotFound eq ""} {
			    lappend L_genesBasedText {*}"[lrepeat $nbColumns($F) ""]"
			} else {		
			    # Change metrics from "." to ","
			    if {[set g_AnnotSV(metrics)] eq "fr"} {
				foreach valueByColumn [split $AnnotFound "\t"] {
				    if {[regexp "^(-)?\[0-9\]+\\.\[0-9\]+$" $valueByColumn]} { ;# Only for the values like "0.5", "-25.3" 
					regsub -all {\.} $valueByColumn "," valueByColumn				    
				    }
				    lappend L_genesBasedText "$valueByColumn"
				} 
			    } else {
				lappend L_genesBasedText {*}[split $AnnotFound "\t"]
			    }
			}
		    }
		}
	    } else {
		######### full lines ###################################
		if {[info exists g_AnnotSV(extann)] && $g_AnnotSV(extann) != ""} {	
		    if {$geneName eq ""} {
			foreach F [ExternalAnnotations L_Files] {
			    lappend L_genesBasedText {*}"[lrepeat $nbColumns($F) ""]"
			}
		    } else {
			set allGenesFromFullLine [split $geneName "/"]
			foreach F [ExternalAnnotations L_Files] {
			    
			    # First, we search for the annotation of each gene that we merge with a "/"
			    set L_AnnotFound {} 
			    foreach g $allGenesFromFullLine {
				set AnnotFound "[ExternalAnnotations $F $g]"
				if {$AnnotFound eq ""} {
				    set AnnotFound "[join [lrepeat $nbColumns($F) ""] "\t"]"
				} 
				lappend L_AnnotFound "$AnnotFound"		
			    }
			    set tutu [MergeAnnotation $L_AnnotFound $nbColumns($F)]

			    # Second (for full lines), we doesn't keep the annotation of all genes (value is set to empty), 
			    # except for "decimal" scores and percentages where we keep the max value
			    # (in order not to have value such "1.5/////-0.2////")
			    # and except for OMIM number where we keep all the numbers
			    set L_newGenesBasedText ""
			    if {$tutu eq ""} {lappend L_genesBasedText ""} ; #if {$tutu eq ""} => doesn't enter in the foreach
			    foreach valueByColumn [split $tutu "\t"] { 
				if {$valueByColumn ne ""} {
				    set isAScore 1
				    set max -1000
				    foreach valueByGene [split $valueByColumn "/"] {
					if {[regexp "ClinGenAnnotations.tsv" $F]} {
					    set max ""
					    # For ClinGen file (HI + TS), the values are ordered as follow: 
					    # 3 > 2 > 1 > 0 > 40 > 30 > Not yet evaluated
					    # We only report 3, 2 or 1 in a full line
					    if {[lsearch -exact {3 2 1} $valueByGene] ne -1} {
						if {$valueByGene eq 3} {set max 3; break}
						if {$valueByGene eq 2} {
						    set max 2
						} elseif {$max ne 2} {
						    set max 1
						} 
					    }
					} elseif {[regexp "OMIM-1-annotations.tsv" $F]} {
					    if {$valueByGene ne ""} {if {$max eq "-1000"} {set max "$valueByGene"} else {append max ";$valueByGene"}}
					} elseif {[regexp "morbid" $F]} {
					    if {[regexp "yes" $valueByGene]} {set max "yes"}
					} else {
					    if {[regexp "^(-)?\[0-9\]{1,3}\\.\[0-9\]+(e-\[0-9\]+)?$" $valueByGene]} { ;# Only for the values like "0.25", "-25.3" 
						if {$valueByGene > $max} {set max $valueByGene}
					    } elseif {$valueByGene ne ""} {set isAScore 0; break}
					} 
				    }
				    if {$isAScore} {
					# Change metrics from "." to ","
					if {[set g_AnnotSV(metrics)] eq "fr"} {
					    regsub -all {\.} $max "," max
					}
					lappend L_newGenesBasedText $max
				    } else {
					lappend L_newGenesBasedText ""
				    }
				} else {
				    lappend L_newGenesBasedText ""
				}
			    }
			    lappend L_genesBasedText {*}$L_newGenesBasedText
			}
		    }
		}
	    }
	    set genesBasedText ""
	    foreach i $g_AnnotSV(genesBasedAnn_i) {
		lappend genesBasedText [lindex $L_genesBasedText $i]
	    }	    
	    set genesBasedText [join $genesBasedText "\t"]
	}

	# GC content annotation
	if {$g_AnnotSV(gcContentAnn)} {
	    if {$AnnotSVtype eq "full"} {	
		set gcContentText "[GCcontentAnnotation $SVchrom $SVleft]"
		append gcContentText "\t[GCcontentAnnotation $SVchrom $SVright]"
	    } else {set gcContentText "\t"}
	}

	# Repeat annotation
	if {$g_AnnotSV(repeatAnn)} {
	    if {$AnnotSVtype eq "full"} {	
		set repeatText "[RepeatAnnotation $SVchrom $SVleft]"
		append repeatText "\t[RepeatAnnotation $SVchrom $SVright]"
	    } else {set repeatText "\t\t\t"}
	}

	# Gap annotation 
	if {$g_AnnotSV(gapAnn)} {
	    if {$AnnotSVtype eq "full"} {	
		set gapText "[GapAnnotation $SVchrom $SVleft]"
		append gapText "\t[GapAnnotation $SVchrom $SVright]"
	    } else {set gapText "\t"}
	}

	# Segmental duplication annotation
	if {$g_AnnotSV(segdupAnn)} {
	    if {$AnnotSVtype eq "full"} {	
		set segdupText "[SegDupAnnotation $SVchrom $SVleft]"
		append segdupText "\t[SegDupAnnotation $SVchrom $SVright]"
	    } else {set segdupText "\t"}
	}

	# ENCODE blacklist annotation
	if {$g_AnnotSV(ENCODEblacklistAnn)} {
	    if {$AnnotSVtype eq "full"} {	
		set ENCODEblacklistText "[ENCODEblacklistAnnotation $SVchrom $SVleft]"
		append ENCODEblacklistText "\t[ENCODEblacklistAnnotation $SVchrom $SVright]"
	    } else {set ENCODEblacklistText "\t\t\t"}
	}

	# TAD annotation
	if {$g_AnnotSV(tadAnn)} {
	    if {$AnnotSVtype eq "split"} {
		set tadText "[TADannotation $SVchrom $intersectStart $intersectEnd $g_AnnotSV(tadAnn_i)]"
	    } else {
		set tadText "[TADannotation $SVchrom $SVleft $SVright $g_AnnotSV(tadAnn_i)]"
	    } 
	}

	# Calculate among the SNV/indel (only for deletion) the:
	# - #hom(sample), #htz(sample), #htz/allHom(sample), #htz/total(sample) variables for each sample
	# - #total(cohort) variable
	set HomHtz ""
	if {$g_AnnotSV(snvIndelFiles) ne ""} {
	    if {$AnnotSVtype eq "split"} {
		set HomHtz "[VCFannotation $SVchrom $intersectStart $intersectEnd $SVtype]"
	    } else {
		set HomHtz "[VCFannotation $SVchrom $SVleft $SVright $SVtype]"
	    } 
	} 

	# Calculate the compound-htz variable
	set compound ""
	if {$g_AnnotSV(candidateSnvIndelFiles) ne ""} {
	    if {$AnnotSVtype eq "split"} {
		set compound [filteredVCFannotation $SVchrom $txStart $txEnd $Ls $headerOutput]
	    } else {
		set compound [filteredVCFannotation "FULL" "" "" "" ""]
	    }
	}

	# "bestAnn" annotation order: 
	# chrom txStart txEnd name2 name cdsStart cdsEnd exonStarts exonEnds
	#
	# headerOutput:
	#  "Gene name\tNumber of genes\ttx\ttx Start\ttx End\ttx length\tnb exons\tCDS length\tframeshift\tlocation\tlocation2\tdistNearestSS\tnearestSStype\tintersectStart\tintersectEnd\tOMIM ID\tOMIM phenotype\tOMIM inheritance\t#hom\t#htz"

	# Insertion of the SV length in the fourth column:
	set SVchrom [lindex $Ls 0]
	set SVstart [lindex $Ls 1]
	set SVend [lindex $Ls 2]
	regsub "chr" [lindex $Ls $i_ref] "" ref
	regsub "chr" [lindex $Ls $i_alt] "" alt

	# Creation of the AnnotSV ID (chrom_start_end_SVtype_i)
	set AnnotSV_ID [settingOfTheAnnotSVID "${SVchrom}_${SVstart}_${SVend}_$SVtype" "$ref" "$alt"]
	# Report of the SV length
	#set SVlength [expr {$SVend-$SVstart}] ; # No! Wrong for an insertion, a BND or a translocation
	if {[info exists g_SVLEN($AnnotSV_ID)]} {
	    set SVlength $g_SVLEN($AnnotSV_ID)
	} else {
	    if {[regexp "DEL" [normalizeSVtype $SVtype]]} { ;# DEL
		set SVlength [expr {$SVstart-$SVend}]
	    } elseif {[normalizeSVtype $SVtype] eq "INV" || $SVtype eq "DUP" || [regexp -nocase "<CN(\[0-9\]+)>" $SVtype]} { ;# DUP or INV
		set SVlength [expr {$SVend-$SVstart}]
	    } else {set SVlength ""}
	}
	
	####### "Exomiser annotation"
	if {$g_AnnotSV(hpo) ne ""} {
	    if {$AnnotSVtype eq "split"} {
		set exomiserText [ExomiserAnnotation $geneName "all"]
	    } else {
		set bestScore "-1.0"
		foreach g [split $geneName "/"] {		    
		    set score [ExomiserAnnotation $g "score"]
		    if {$score > $bestScore} {set bestScore $score}
		}
		set exomiserText "$bestScore\t\t\t"
	    }
	}

	
	#################################################################################
	################### creation of the "L_TextToWrite" variable ####################
	#################################################################################
	
	set TextToWrite ""

	####### "Basic SV annotations"
	if {$g_AnnotSV(SVinputInfo)} {
	    set toadd [lrange $Ls 0 [expr {$theBEDlength-1}]]
	    set toadd [linsert $toadd 3 $SVlength]
	    append TextToWrite "$AnnotSV_ID\t[join $toadd "\t"]"
	} else {
	    append TextToWrite "$AnnotSV_ID\t[join [lrange $Ls 0 2] "\t"]\t$SVlength"
	    if {$g_AnnotSV(svtBEDcol) ne -1} { ; # SV type is required for the ranking
		append TextToWrite "\t$SVtype"
		if {$g_AnnotSV(samplesidBEDcol) ne -1} {
		    append TextToWrite "\t$Samplesid"  
		}
	    } elseif {$g_AnnotSV(samplesidBEDcol) ne -1} {
		append TextToWrite "\t$Samplesid"
	    } 
	}
	append TextToWrite "\t$AnnotSVtype"

	####### "Basic gene annotations"
	append TextToWrite "\t$geneName\t$NbGenes\t$transcript\t$txStart\t$txEnd\t$txL\t$CDSl\t$frameshift\t$nbExons\t$location\t$location2\t$distNearestSS\t$nearestSStype\t$intersect"

	####### "Regulatory elements annotations"
	append TextToWrite "\t$reText"
	
	#######  "Annotations with benign genes or genomic regions (SVincludedInFt)"
	append TextToWrite "\t$benignText"
	
	#######  "Users: SVincludedInFt"
	if {[glob -nocomplain $usersDir/SVincludedInFt/*.formatted.sorted.bed] ne ""} { ; # Don't put {$SVincludedInFTtext ne ""}: the user BED could have only 1 annotation column, and so $UserText can be equel to "" (without "\t")
	    append TextToWrite "\t$SVincludedInFTtext"
	}
	
	#######  "Annotations with pathogenic genes or genomic regions (FtIncludedInSV)"
	append TextToWrite "\t$pathogenicText"
	
	####### "FtIncludedInSV"
	if {$g_AnnotSV(tadAnn)} {
	    append TextToWrite "\t$tadText"
	}
	if {$HomHtz ne ""} {
	   append TextToWrite "\t$HomHtz"
	} 
	if {$g_AnnotSV(candidateSnvIndelFiles) ne ""} {
	   append TextToWrite "\t$compound"
	}
	if {[glob -nocomplain $usersDir/FtIncludedInSV/*.formatted.sorted.bed] ne ""} { ; # Don't put {$FtIncludedInSVtext ne ""}: the user BED could have only 1 annotation column, and so $UserText can be equel to "" (without "\t")
	    append TextToWrite "\t$FtIncludedInSVtext"
	}
	
	####### "Breakpoints annotations"
	if {$g_AnnotSV(gcContentAnn)} {
	    append TextToWrite "\t$gcContentText"
	}
	if {$g_AnnotSV(repeatAnn)} {
	    append TextToWrite "\t$repeatText"
	}
	if {$g_AnnotSV(gapAnn)} {
	    append TextToWrite "\t$gapText"
	}
	if {$g_AnnotSV(segdupAnn)} {
	    append TextToWrite "\t$segdupText"
	}
	if {$g_AnnotSV(ENCODEblacklistAnn)} {
	    append TextToWrite "\t$ENCODEblacklistText"
	}	
	####### "Genes-based annotations"
	if {$g_AnnotSV(genesBasedAnn)} {
	    append TextToWrite "\t$genesBasedText"
	}
	####### "Exomiser annotation"
	if {$g_AnnotSV(hpo) ne ""} {
	    append TextToWrite "\t$exomiserText"
	}
	####### "Ranking annotations"  ( => done through "$TextToWrite" )
	if {$g_AnnotSV(ranking)} {
	    set rank [SVranking $TextToWrite]
	    # To select the SV of a user-defined specific class (from 1 to 5)
	    if {[lsearch -exact $g_AnnotSV(rankFiltering) $rank] eq -1} {
		continue
	    }
	    if {[lsearch -exact "$g_AnnotSV(outputColHeader)" "ranking decision criteria"] ne -1} {
		append TextToWrite "\t$rank\t$g_rankingExplanations($AnnotSV_ID)"
	    } else {
		append TextToWrite "\t$rank"
	    }
	}

	###### End of 1 annotation line
	lappend L_TextToWrite "$TextToWrite"	
    }
    close $f

    
    ################################################
    ################### Writing ####################
    ################################################
    puts "\n...writing of $outputFile ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    WriteTextInFile [join $L_TextToWrite "\n"] "$outputFile"


    
    ##############################################################
    ################### Delete temporary file ####################
    ##############################################################
    file delete -force $FullAndSplitBedFile
    file delete -force $g_AnnotSV(bedFile) ; # => ".formatted.bed" tmp file
    file delete -force $g_AnnotSV(fullAndSplitBedFile)
    regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".breakpoints.bed" tmpBreakpointsFile
    file delete -force $tmpBreakpointsFile
    foreach tmpBedFile [glob -nocomplain $g_AnnotSV(outputDir)/*_AnnotSV_inputSVfile.bed] {
	file delete -force $tmpBedFile ; # => a bedfile is present only if "-SVinputFile" is a VCF
    }
    if {[info exists headerFileToRemove]} {
	regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".header.tsv" BEDinputHeaderFile
	file delete -force $BEDinputHeaderFile
    }   

    ################################################
    ################### Display ####################
    ################################################
    puts "\n...output columns annotation ([clock format [clock seconds] -format "%B %d %Y - %H:%M"]):"
    regsub -all "\t" $headerOutput "; " t
    puts "\t$t\n"

}
