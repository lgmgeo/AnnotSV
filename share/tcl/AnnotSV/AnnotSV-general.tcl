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


proc DescendingSortOnElement1 {X Y {N 1}} {
    return [expr {[lindex $X $N]<[lindex $Y $N]}]
}
proc AscendingSortOnElement0 {X Y {N 0}} {
    return [expr {[lindex $X $N]>[lindex $Y $N]}]
}
proc AscendingSortOnElement1 {X Y {N 1}} {
    return [expr {[lindex $X $N]>[lindex $Y $N]}]
}
proc AscendingSortOnElement2 {X Y {N 1}} {
    return [expr {[lindex $X $N]>[lindex $Y $N]}]
}
proc AscendingSortOnElement4 {X Y {N 4}} {
    return [expr {[lindex $X $N]>[lindex $Y $N]}]
}
proc AscendingSortOnElement5 {X Y {N 5}} {
    return [expr {[lindex $X $N]>[lindex $Y $N]}]
}

##############################################################################
#                          WORKING WITH FILES
##############################################################################
proc FirstLineFromFile {{File ""}} {

    if {[regexp ".gz$" $File]} {
	set F [open "|gzip -cd $File"] 
    } else {
	set F [open "$File"]
    }

    set L_Lines {}
    set First 1

    while {[gets $F Line]>=0} {
	if {$First} {set First 0;break}
    }

    return $Line
}

proc ContentFromFile {{File ""}} {
    if {[string equal $File ""]} {return ""}
    set f     [open $File r]
    set Text [read -nonewline $f]
    close $f
    return $Text
}

proc LinesFromFile {{File ""}} {
    return [split [ContentFromFile $File] "\n"]
}

proc ContentFromGZFile {{File ""}} {
    if {[string equal $File ""]} {return ""}
    set f     [open "|gzip -cd $File" r]
    set Texte [read -nonewline $f]
    close $f
    return $Texte
}

proc LinesFromGZFile {{File ""}} {
    return [split [ContentFromGZFile $File] "\n"]
}

proc ReplaceTextInFile {NewText fichier} {
    set    fifi [open $fichier w]
    puts  $fifi $NewText
    close $fifi 
    return 1
}


proc WriteTextInFile {text fichier} {
    set    fifi [open $fichier a]
    puts  $fifi $text
    close $fifi 
    return 1
}



##############################################################################
#                          Settings of the annotSV_ID
##############################################################################

# INPUT: "${chrom}_${pos}_${end}_${svtype}" "$ref" "$alt"
proc settingOfTheAnnotSVID {deb ref alt} {
    global g_ID
    global g_deb

    if {![info exists g_ID($deb,$ref,$alt)]} {
	set i 1	
	while {[info exists g_deb(${deb}_$i)]} {incr i; continue}
	set g_deb(${deb}_$i) "1"
	set g_ID($deb,$ref,$alt) "${deb}_$i"
    }
    return $g_ID($deb,$ref,$alt)
}



##############################################################################
#                          WORKING WITH rsID
##############################################################################

proc isNotAnRS {rs} {
    if {[regexp "^rs\[0-9\]+$" $rs]} {return 0} else {return 1}
}

##############################################################################
#                     WORKING WITH ANNOTATION LIST
##############################################################################

# Example: 
# 1 annotation file (F1) with 2 annotation columns (C1 and C2) by gene
# (numberOfAnnCol = 2)
#
# Input (= a list):
# -----------------
# Annotations for the gene1: g1C1 g1C2
# Annotations for the gene2: g2C1 g2C2 
# Annotations for the gene3:  ""   ""  

# <=> Input = L_ann = {g1C1\tg1C2} {g2C1\tg2C2}} {\t}
#
# Output (= a tsv text):
# ----------------------
# Annotation desired for g1/g2/g3: 
# <=> Output = L_newAnn = "g1C1/g2C1/g3C1\tg1C2/g2C2/g3/C2"
# 
proc MergeAnnotation {L_ann numberOfannCol} {

    # There is only one gene annotation => no need to merge
    if {[llength $L_ann] eq 1} {
	return [join $L_ann]
    }

    # Merge genes annotations
    set g 0
    set L_ig {}
    foreach geneAnn $L_ann {          ;# geneAnn = g1C1\tg1C2
	lappend L_ig $g
	if {$geneAnn eq ""} {
	    set tann($g,0) "" ; # With 1 annotation column: An empty string can not be split and can not return a colAnn value => don't enter in the following foreach
	} else {
	    set c 0
	    foreach colAnn [split $geneAnn "\t"] {    ;# colAnn = g1C1
		set tann($g,$c) $colAnn               ;# tann($g,$c) = g1C1
		incr c
	    }
	}
	incr g
    }
    set L_newAnn {}
    set c 0
    while {$c < $numberOfannCol} {
	set annByCol ""
	foreach g $L_ig {
	    lappend annByCol $tann($g,$c)
	}
	incr c 
	set 1ann [join $annByCol "/"]
	if {[regexp "^( |/)+$" $1ann]} {set 1ann ""}
	lappend L_newAnn $1ann
    }
    set toreturn [join $L_newAnn "\t"]

    return $toreturn
}


##############################################################################
#                          WORKING WITH SV TYPE
##############################################################################

proc normalizeSVtype {SVtype} {

    # SVtype in which category: DUP? DEL? INV? INS? None?
    if {[regexp -nocase "del|loss|<CN0>" $SVtype]} {
	set SVtype "DEL"
    } elseif {[regexp -nocase "dup|gain|MCNV" $SVtype ]} {
	set SVtype "DUP"
    } elseif {[regexp -nocase "<CN(\[0-9\]+)>" $SVtype match i]} {
	if {$i>1} {set SVtype "DUP"}
    } elseif {[regexp -nocase "inv" $SVtype]} {
	set SVtype "INV"
    } elseif {[regexp -nocase "ins|MEI|alu|line|sva" $SVtype]} { ;# "DEL_ALU" is set to "DEL", OK!
	set SVtype "INS"
    } else {
	set SVtype "None"
    }

    return $SVtype
}

##############################################################################
#                          WORKING WITH bedFile
##############################################################################

# - Create a formatted file (*.formatted.bed)
# - If the bed file is located in the "Users" directory, create also a *.header.tsv file (if it doesn't exist yet)
proc checkBed {bedFile {bedDir ""} {mergeOverlap 0}} {

    global g_AnnotSV

    regsub -nocase ".bed$" $bedFile ".formatted.bed" newBed
    if {$bedDir eq ""} {
	set newBed "$g_AnnotSV(outputDir)/[file tail $newBed]"
    } else {
	set newBed "$bedDir/[file tail $newBed]"
    }

    if {[file exists $newBed]} {return 1}


    # Bedfile should be sorted, should not have "chr" in the first column, should not have 2 identical lines 
    ##############################################################################################################################################
    # If mergeOverlap is set to 1, should not have overlapping lines

    # Removing non-standard contigs (other than the standard 1-22,X,Y,MT) and remove "chr" if present
    # Check that the length is the same in each line of the bed file
    set f [open $bedFile]
    set test 0
    set L_length {}
    set header ""
    while {![eof $f]} {
        set L [gets $f]
        if {$L eq ""} {continue}
        if {[regexp "^#" $L]} {set header $L; continue}

        if {[regsub "chr" [lindex $L 0] "" chrom] ne 0} {
            lappend L_Text($chrom) "[string range $L 3 end]"
        } else {
            lappend L_Text($chrom) "$L"
        }
	lappend L_length [llength [split $L "\t"]]
    }
    close $f
    # With or without annotations?
    set L_length [lsort -unique $L_length]
    if {[llength $L_length] ne 1} {
	puts "$bedFile: Unexpected file format, not the same length for each line."
	puts "\t-> not using the associated annotations"
	set UseAnn 0
    } else {
	if {$L_length > 3} {set UseAnn 1} else {set UseAnn 0}
    }

    # Creation of the *.header.tsv file (if the bed file is located in the "Users" directory)
    if {[regexp "/Users/" $bedFile]} {
	regsub -nocase "(.formatted)?.bed$" $bedFile ".header.tsv" userHeaderFile
	if {![file exists $userHeaderFile]} {
	    if {[llength [split $header "\t"]] eq $L_length} {
		WriteTextInFile $header $userHeaderFile
	    } else {
		WriteTextInFile "[join [lrepeat $L_length ""] "\t"]" $userHeaderFile
	    }
	}
    }
    
    foreach chrom {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M MT} {
        if {![info exists L_Text($chrom)]} {continue}

        ## Remove identical lines (due to liftover_to_hg19 for example)
        set L_Text($chrom) [lsort -unique $L_Text($chrom)]

        ## Put in karyotypic order
        set L_Text($chrom) [lsort -command AscendingSortOnElement1 [lsort -command AscendingSortOnElement2 $L_Text($chrom)]]

	if {$mergeOverlap} {
	    ## Merge overlapping lines
	    set L_newText {}
	    set Ls [split [lindex $L_Text($chrom) 0] "\t"]
	    set previousStart [lindex $Ls 1]
	    set previousEnd [lindex $Ls 2]
	    if {$UseAnn} {
		set j 3
		while {$j < $L_length} {
		    set previousAnnotations($j) "[lindex $Ls $j]"
		    incr j
		}
	    }
	    foreach L [lrange $L_Text($chrom) 1 end] {
		set Ls [split $L "\t"]
		set start [lindex $Ls 1]
		set end [lindex $Ls 2]
		if {$UseAnn} {
		    set j 3
		    while {$j < $L_length} {
			set annotations($j) "[lindex $Ls $j]"
			incr j
		    }
		}
		if {$start < $previousEnd} {
                    if {$end > $previousEnd} {set previousEnd $end}
		    if {$UseAnn} {
			set j 3
			while {$j < $L_length} {
			    set previousAnnotations($j) "$previousAnnotations($j)/$annotations($j)"
			    incr j
			}
		    } 
		} else {               
		    if {$UseAnn} {
			set n "$chrom\t$previousStart\t$previousEnd"
			set j 3
			while {$j < $L_length} {
			    append n "\t$previousAnnotations($j)"
			    incr j
			}
			lappend L_newText "$n"
		    } else {
			lappend L_newText "$chrom\t$previousStart\t$previousEnd"
		    }
		    set previousStart $start
		    set previousEnd $end
		    set j 3
		    while {$j < $L_length} {
			set previousAnnotations($j) "$annotations($j)"
			incr j
		    }		    
		}
	    }
	    # Treatment of the last line(s)
	    if {$UseAnn} {
		set n "$chrom\t$previousStart\t$previousEnd"
		set j 3
		while {$j < $L_length} {
		    append n "\t$previousAnnotations($j)"
		    incr j
		}
		lappend L_newText "$n"
	    } else {
		lappend L_newText "$chrom\t$previousStart\t$previousEnd"
	    }
	}

        ## Write the formatted bedfile
	if {$mergeOverlap} {
	    WriteTextInFile [join $L_newText "\n"] $newBed
	} else {
	    WriteTextInFile [join $L_Text($chrom) "\n"] $newBed
	}
        unset L_Text($chrom)
    }
    file delete -force $bedFile


    return $newBed
}



# Check if a bed or VCF file is empty
proc isAnEmptyFile {bedOrVCFfile} {

    # Return 1 if the extension is not .bed or .vcf
    if {![regexp -nocase "\\.(bed|vcf(.gz)?)$" $bedOrVCFfile]} {return 1}

    # After filtering lines beginning with "#", at least 1 line should be still present
    if {[regexp ".gz$" $bedOrVCFfile]} {
	set f [open "|gzip -cd $bedOrVCFfile"] 
    } else {
	set f [open "$bedOrVCFfile"]
    }
    set test 0
    while {![eof $f]} {
        set L [gets $f]
        if {$L eq ""} {continue}
        if {[regexp "^#" $L]} {continue}
	incr test
	break
    }

    if {$test} {return 0} else {return 1}
}


# Check the -svtBEDcol option (=number of the SVtype column)
# - Should be > 3 (the 3 first columns are chrom, start, end)
# - Should be an existing column (<= total number of columns)
proc checksvtBEDcol {bedFile} {

    global g_AnnotSV

    # Look at the length of a line in the bed file
    set f [open $bedFile]
    while {![eof $f]} {
        set L [gets $f]
        if {$L eq ""} {continue}
        if {[regexp "^#" $L]} {continue}

	set thelength [llength [split $L "\t"]]
	break
    }
    close $f

    # Check the -svtBEDcol option, is set to -1 if the column doesn't exist
    if {$g_AnnotSV(svtBEDcol) ne -1} {
	if {$g_AnnotSV(svtBEDcol) > $thelength} {
	    puts "\nWARNING: -svtBEDcol = $g_AnnotSV(svtBEDcol)"
	    puts "This value should correspond to a column of $bedFile ($thelength columns)"
	    puts "-svtBEDcol is set to -1. NO SV RANKING WILL BE PERFORMED.\n"
	    set g_AnnotSV(svtBEDcol) "-1"
	    set g_AnnotSV(ranking) 0
	} elseif {$g_AnnotSV(svtBEDcol) < 4} {
	    puts "\nWARNING: -svtBEDcol = $g_AnnotSV(svtBEDcol)"
	    puts "This value should be > 3"
	    puts "-svtBEDcol is set to -1. NO SV RANKING WILL BE PERFORMED.\n"
	    set g_AnnotSV(svtBEDcol) "-1"
	    set g_AnnotSV(ranking) 0
	} else {
	    # 0 is the first element of an informatic list: -1
	    # AND 2 columns (AnnotSV ID + SV length) are added: +2
	    # Finally: +1
	    incr g_AnnotSV(svtBEDcol) 1
	}
    }
}
