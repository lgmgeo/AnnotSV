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



##################################################################
# Prepare 2 new annotation columns for annotated bedfile:
# - nb Hom in the SV
# - nb Htz in the SV
##################################################################


proc DichotomySearch {number OrderedNumberList} {
    # Search for the place of a number into a big ordered number list. 
    # Return the indice after which we can range the number
    # Return -1 if number < [lindex $ListeOrdonnee 0]
    set i 0
    set j [expr {[llength $OrderedNumberList]-1}]
    if {$number<[lindex $OrderedNumberList 0]} {return -1}
    if {$number>=[lindex $OrderedNumberList $j]} {return $j}
    set k [expr {($i+$j)/2}]

    if {$number<[lindex $OrderedNumberList $k]} {set j $k} else {set i $k}
    while {[expr {$j-$i}] > 1} {
        set k [expr {($i+$j)/2}]
	if {$number<[lindex $OrderedNumberList $k]} {set j $k} else {set i $k}
    }
    return $i
}

proc VCFannotation {SVchrom SVstart SVend} {

    global g_AnnotSV
    global lPos

    ## VCFs parsing is done only 1 time (After, g_AnnotSV(vcfParsing) is set to "done")
    if {![info exists g_AnnotSV(vcfParsing)]} {
	set g_AnnotSV(vcfParsing) "done"
	# parsing of $g_AnnotSV(snvIndelFiles): creation of lPos($SVchrom,htz,sample) and lPos($SVchrom,hom,sample)
	puts "\n\n...parsing of VCF file(s) for \"$g_AnnotSV(snvIndelSamples)\" ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])" 

	# "eval glob" accept regular expression ("*.vcf) as well as a list of files ("sample1.vcf sample2.vcf.gz"):
	foreach vcfF [eval glob -nocomplain $g_AnnotSV(snvIndelFiles)] {
	    puts "\t...parsing of $vcfF ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])" 
	    set iIntersect 0
	    set iSV 0
	    set iLoaded 0
	    set iNotPASS 0
	    set iNotGT 0
	    set chrManipulation 0

	    # Intersect each VCF with the bedfile
	    if {[regexp ".gz$" $vcfF]} {
		regsub ".gz$" [file tail $vcfF] ".tmp.gz" tmpVCFgz
		set tmpVCFgz "$g_AnnotSV(outputDir)/$tmpVCFgz"
		regsub ".gz$" $tmpVCFgz "" tmpVCF

		# VCF should not contain "chr" prefix
		catch {exec gunzip -c $vcfF | grep -c "^chr"} Message
		if {[string index $Message 0] ne 0} {
		    exec gunzip -c $vcfF | sed "s/^chr//" | gzip  > $g_AnnotSV(outputDir)/[file tail $vcfF].withoutchr.vcf.gz
		    set vcfF $g_AnnotSV(outputDir)/[file tail $vcfF].withoutchr.vcf.gz
		    set chrManipulation 1
		}

		set f [open "| gzip -cd $vcfF"]	 
		set L_tmpVCF_Header {} ; # Header is used to find "samples" columns names
		while {! [eof $f]} {
		    set L [gets $f]
		    if {[regexp "^#" $L]} {
			set L_tmpVCF_Header $L
		    } else {
			break
		    }
		}
		unset f
		
		WriteTextInFile "$L_tmpVCF_Header" $tmpVCF
		if {[catch {exec gzip $tmpVCF} Message]} {
		    puts "-- VCFannotation --"
		    puts "gzip $tmpVCF"
		    puts $Message
		}
		if {[catch {exec gunzip -c $vcfF | $g_AnnotSV(bedtools) intersect -a stdin -b $g_AnnotSV(bedFile) > $tmpVCFgz.bis} Message]} {
		    if {[regexp -nocase "error" $Message]} {
			WriteTextInFile "$Message" $vcfF.intersect.error
			puts "-- VCFannotation --"
			puts "gunzip -c $vcfF | $g_AnnotSV(bedtools) intersect -a stdin -b $g_AnnotSV(bedFile) > $tmpVCFgz.bis"
			puts "Error: look at $vcfF.intersect.error"
		    }
		}	
	    } else {
		set tmpVCFgz "$g_AnnotSV(outputDir)/[file tail $vcfF].tmp.gz"
		set tmpVCF "$g_AnnotSV(outputDir)/[file tail $vcfF].tmp"
		
		# VCF should not contain "chr" prefix
		catch {exec grep -c "^chr" $vcfF} Message
		if {[string index $Message 0] ne 0} {
		    exec sed "s/^chr//" $vcfF > $g_AnnotSV(outputDir)/[file tail $vcfF].withoutchr.vcf
		    set vcfF $g_AnnotSV(outputDir)/[file tail $vcfF].withoutchr.vcf
		    set chrManipulation 1
		}

		set f [open "$vcfF"]	 
		set L_tmpVCF_Header {} ; # Header is used to find "samples" columns names
		while {! [eof $f]} {
		    set L [gets $f]
		    if {[regexp "^#" $L]} {
			set L_tmpVCF_Header $L
		    } else {
			break
		    }
		}
		close $f

		WriteTextInFile "$L_tmpVCF_Header" $tmpVCF
		if {[catch {exec gzip $tmpVCF} Message]} {
		    puts "-- VCFannotation --"
		    puts "gzip $tmpVCF"
		    puts $Message
		}
		if {[catch {exec $g_AnnotSV(bedtools) intersect -a $vcfF -b $g_AnnotSV(bedFile) > $tmpVCFgz.bis} Message]} {
		    if {[regexp -nocase "error" $Message]} {
			WriteTextInFile "$Message" "$g_AnnotSV(outputDir)/[file tail $vcfF].intersect.error"
			puts "-- VCFannotation --"
			puts "$g_AnnotSV(bedtools) intersect -a $vcfF -b $g_AnnotSV(bedFile) > $tmpVCFgz.bis"
			puts "Error: look $g_AnnotSV(outputDir)/[file tail $vcfF].intersect.error"
		    }
		}
	    }

	    if {$chrManipulation} {file delete -force $vcfF}

	    if {[file size $tmpVCFgz.bis] eq 0} {
		# no intersection between the bed and the VCF, no SNV/Indel to load in memory
		puts "\t\t-> no variant in the intersection"
		file delete -force $tmpVCFgz
		file delete -force $tmpVCFgz.bis
		continue
	    } else {
		if {[catch {exec cat $tmpVCFgz.bis | gzip >> $tmpVCFgz} Message]} {
		    puts "-- VCFannotation --"
		    puts "cat $tmpVCFgz.bis | gzip >> $tmpVCFgz"
		    puts $Message
		}
		file delete -force $tmpVCFgz.bis
	    }

	    # load SNV/Indel in memory
	    set f [open "| gzip -cd $tmpVCFgz"]	  
	    set L_samples {}
	    while {![eof $f]} {
		set L [gets $f]

		if {[string range $L 0 5] eq "#CHROM"} {
		    foreach sample $g_AnnotSV(snvIndelSamples) {
			set i_sample($sample) [lsearch -exact [split $L "\t"] $sample]
			if {$i_sample($sample) ne -1} {lappend L_samples $sample}
		    }
		    continue
		}
		if {[string index $L 0] eq "#" || $L eq ""} {continue}

		incr iIntersect

		# set the variables
		set Ls [split $L "\t"] 
		set chrom  [lindex $Ls 0]
		regsub -nocase "chr" $chrom "" chrom
		set pos    [lindex $Ls 1]

		# Consider only the SNV/indel (not the SV in the VCF file)
		##########################################################
		# Example of SV: 
		# - Type1: ref="G" and alt="ACTGCTAACGATCCGTTTGCTGCTAACGATCTAACGATCGGGATTGCTAACGATCTCGGG"
		# - Type2: "<INS>", "<DEL>", ...
		# - Type3: complex rearrangements with breakends: "G]17:1584563]"
		set ref [lindex $Ls 3]
		set alt [lindex $Ls 4]
		if {[regexp "<|\\\[|\\\]" $alt]} {incr iSV; continue}; # it is an SV
		set variantLength [expr {[string length $ref]-[string length $alt]}]
		if {[expr {abs($variantLength)}]>$g_AnnotSV(SVminSize)} {incr iSV; continue}; # it is an SV

		set filter [lindex $Ls 6]
		set formatData [lindex $Ls 8]
		set j_GT [lsearch -exact [split $formatData ":"] "GT"]
		if {$j_GT eq -1} {incr iNotGT; continue}
	      
		# keep only variant with FILTER == PASS
		if {$g_AnnotSV(snvIndelPASS) && $filter ne "PASS"} {incr iNotPASS; continue}
		incr iLoaded
		
		foreach sample $L_samples {
		    set sampleData [lindex $Ls $i_sample($sample)]
		    # set the GT
		    set GTsample [lindex [split $sampleData ":"] $j_GT]
		    set GTsample [split $GTsample "/|\\|"]
		    set GTsampleA [lindex $GTsample 0]
		    set GTsampleB [lindex $GTsample 1]
		    if {$GTsampleA eq "" || $GTsampleB eq ""} {continue}
		    if {$GTsampleA eq "0" && $GTsampleB eq "0"} {continue}
		    if {[lindex $GTsample 0] ne [lindex $GTsample 1]} {set GT "htz"} else {set GT "hom"} 
		    # set the lPos($chrom,$GT,$sample)
		    lappend lPos($chrom,$GT,$sample) $pos
		}
	    }
	    close $f
	    file delete -force $tmpVCFgz
	    if {$iIntersect} {
		puts "\t\t-> $iIntersect variants located in the SV"
	    }
	    if {$iLoaded} {
		puts "\t\t-> $iLoaded variants loaded"
	    }
	    if {$iSV} {
		puts "\t\t-> $iSV SV excluded (considering only SNV/indel from the VCF)"
	    }
	    if {$iNotPASS} {
		puts "\t\t-> $iNotPASS variants excluded beacause of the FILTER value not equal to \"PASS\""
	    }
	    if {$iNotGT} {
		puts "\t\t-> $iNotGT variants excluded because of the absence of GT information"
	    }
	}
    }

    set textToReturn {}
    foreach sample $g_AnnotSV(snvIndelSamples) {
	# if $chrom not present in the VCF file:
	if {![info exists lPos($SVchrom,htz,$sample)] && ![info exists lPos($SVchrom,hom,$sample)]} {lappend textToReturn "0\t0"; continue}
	# Count for hom variants
	set countHom 0
	if {[info exists lPos($SVchrom,hom,$sample)]} {
	    set i_first [DichotomySearch $SVstart $lPos($SVchrom,hom,$sample)]
	    set i_last [DichotomySearch $SVend $lPos($SVchrom,hom,$sample)]
	    if {$SVstart ne [lindex $lPos($SVchrom,hom,$sample) $i_first]} {set i_first [expr {$i_first+1}]}
	    set L_posHom [lsort -unique [lrange $lPos($SVchrom,hom,$sample) $i_first $i_last]]
	    set countHom [llength $L_posHom]
	} 
	# Count for htz variants
	set countHtz 0
	if {[info exists lPos($SVchrom,htz,$sample)]} {
	    set i_first [DichotomySearch $SVstart $lPos($SVchrom,htz,$sample)]
	    set i_last [DichotomySearch $SVend $lPos($SVchrom,htz,$sample)]
	    if {$SVstart ne [lindex $lPos($SVchrom,htz,$sample) $i_first]} {set i_first [expr {$i_first+1}]}
	    set L_posHtz [lsort -unique [lrange $lPos($SVchrom,htz,$sample) $i_first $i_last]]
	    set countHtz [llength $L_posHtz]
	}
	lappend textToReturn "$countHom\t$countHtz"
    }

    return "[join $textToReturn "\t"]" 
}


proc VCFsToBED {SV_VCFfiles} {

    global g_AnnotSV
    global VCFheader
    global g_SVLEN

    set SV_BEDfile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_AnnotSV_inputSVfile.bed"
    file delete -force "$SV_BEDfile"
    set VCFheader "" 
    # Variants from the input file that are not annotated by AnnotSV are reported in $unannotatedOutputFile
    if {[regexp "annotated" $g_AnnotSV(outputDir)/$g_AnnotSV(outputFile)]} {
	regsub "annotated" $g_AnnotSV(outputDir)/$g_AnnotSV(outputFile) "unannotated" unannotatedOutputFile
    } else {
	regsub ".tsv$" $g_AnnotSV(outputDir)/$g_AnnotSV(outputFile) ".unannotated.tsv" unannotatedOutputFile
    }
    file delete -force $unannotatedOutputFile

    foreach VCFfile $SV_VCFfiles {
	set L_TextToWrite {}
	
	if {[regexp ".gz$" $VCFfile]} {
	    set f [open "| gzip -cd $VCFfile"]	 
	} else {		
	    set f [open "$VCFfile"]
	}	 
	
	set i 0
	set VCFheaderNotPresent 1
	while {![eof $f]} {
	    if {$i eq "500000"} {
		WriteTextInFile [join $L_TextToWrite "\n"] $SV_BEDfile
		set L_TextToWrite {}	    
		set i 0
	    }
	    set L [gets $f]
	    set Ls [split $L "\t"]

	    if {[string index $L 0] eq "#" || $L eq ""} {
		if {[regexp "^#CHROM" $L]} {
		    set VCFheaderNotPresent 0
		    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG00096
		    if {$g_AnnotSV(SVinputInfo)} {
			set VCFheader "SV type\t[join [lrange $Ls 2 end] "\t"]"
		    } else {
			set VCFheader "SV type\t[join [lrange $Ls 3 4] "\t"]\t[join [lrange $Ls 8 end] "\t"]"
		    }
		}
		continue
	    }
	    incr i
	    
	    if {$VCFheaderNotPresent} {
		set VCFheaderNotPresent 0
		puts "WARNING:\n$VCFfile: no VCF header line (prefixed with \"#CHROM\"). Check your VCF.\n"
	    }

	    # Consider only the SV (not the SNV/indel)
	    ##########################################
	    # Example of SV: 
	    # - Type1: ref="G" and alt="ACTGCTAACGATCCGTTTGCTGCTAACGATCTAACGATCGGGATTGCTAACGATCTCGGG" (length > 50bp)
	    # - Type2: alt="<INS>", "<DEL>", ...
	    # - Type3: complex rearrangements with breakends: alt="G]17:1584563]" or alt="G]chr17:1584563]"
	    set chrom [lindex $Ls 0]
	    regsub -nocase "chr" $chrom "" chrom
            set pos [lindex $Ls 1]
	    regsub "chr" [lindex $Ls 3] "" ref
	    regsub "chr" [lindex $Ls 4] "" alt 
	    regsub -all "\"" [lindex $Ls 7] "" INFOcol
            if {[regexp "SVLEN=(\[0-9-\]+)" $INFOcol match SVLEN]} {set svlen $SVLEN} else {set svlen ""}
	    if {[regexp ";END=(\[0-9\]+)" $INFOcol match END] || [regexp "^END=(\[0-9\]+)" $INFOcol match END]} {set end $END} else {set end ""}
	    if {[regexp "SVTYPE=(\[^;\]+)" $INFOcol match SVTYPE]} {set svtype $SVTYPE} else {set svtype ""}      
	    if {[regexp "^<" $alt]} {
		# Type2
		if {$end eq ""} {
		    # INS:ME (LINE1, ALU or SVA)
		    set end [expr {$pos+1}]
		} 
		if {[regexp "<TRA>" $alt]} {
		    # Example of a strange VCF format line (translocation):
		    # chr1    63705386   N    .     <TRA>   51      PASS    "PRECISE;CT=5to5;CIPOS=-10,10;CIEND=-10,10;SOMATIC;SVTYPE=TRA;CHR2=chrX;END=478444;SVLEN=0"        GT      ./.   
		    # POS = 63705386 --> sur chrom1
		    # END = 478444   --> sur chromX
		    ## => annotation only of the first breakpoint
		    set end [expr {$pos+1}]
		}
	    } elseif {[regexp "(\\\[|\\\])(\[^:\]+):(\[0-9\]+)" $alt]} {
		# Type3 	
		# Only one breakend is annotated with this kind of line
		set end [expr {$pos+1}]		
	    } elseif {[regexp -nocase "^\[ACGTN.*\]+$" $ref$alt]} {
		# Type1
		regsub -all "\[*.\]" $ref "" refbis ;# cf GRIDSS comment, just below
		regsub -all "\[*.\]" $alt "" altbis
		set variantLengthType1 [expr {[string length $altbis]-[string length $refbis]}]
		if {[expr {abs($variantLengthType1)}]<$g_AnnotSV(SVminSize)} {
		    # it is an indel
		    set AnnotSV_ID [settingOfTheAnnotSVID "${chrom}_${pos}_${end}_${svtype}" "$ref" "$alt"]
		    WriteTextInFile "$AnnotSV_ID: variantLength ([expr {abs($variantLengthType1)}]) < SVminSize ($g_AnnotSV(SVminSize))" $unannotatedOutputFile
		    continue
		}
		if {$variantLengthType1>0} {
		    # insertion
		    if {$end eq ""} {set end [expr {$pos+1}]} 
		    if {$svtype eq ""} {set svtype "INS"}
		} else {
		    # deletion
		    if {$end eq ""} {set end [expr $pos-$variantLengthType1]}
		    if {$svtype eq ""} {set svtype "DEL"}
		}
		if {[regexp "\\\." $ref] || [regexp "\\\." $alt]} {
		    # The GRIDSS author says that a . followed by bases refers to a single breakend where the reads cannot be uniquely mapped back. 
		    # e.g.: 2       39564894        gridss28_45b    T       .TTCTCTCATAACAAACCATGACATCCAGTCATTTAATACAATATGTCTGGGGTGGCTGGGCCCCTTTTTT 246.24  LOW_QUAL
		    # => AnnotSV can not determine the SV length
		    set variantLengthType1 ""
		} elseif {[string length $refbis] < 10000 && [string length $altbis] < 10000} { ;# else we can have "couldn't compile regular expression pattern: out of memory" in the next regexp!!!
		    if {![regexp -nocase "$refbis" $altbis] && ![regexp -nocase "$altbis" $refbis]} {
			# Complex SV: AGT>ATTGCATGGACCTGAGTCCCCAAAAAAAAAAATTTTTTTTTTGGGGGGGGGGCCCCCCCCCCGGGGGGGGGG 
			# => AnnotSV can not determine the SV length
			set variantLengthType1 ""
		    }
		} else {
		    # AnnotSV can not check the ref and alt values
		    # => AnnotSV can not determine the SV length
		    set variantLengthType1 ""		    
		}
		if {$svlen eq ""} {set svlen $variantLengthType1}
	    } else {
		set AnnotSV_ID [settingOfTheAnnotSVID "${chrom}_${pos}_${end}_${svtype}" "$ref" "$alt"]
		WriteTextInFile "$AnnotSV_ID: not a known SV format" $unannotatedOutputFile
		continue
	    }

	    if {$end < $pos} {set tutu $end; set end $pos; set pos $tutu}
	    if {$end eq $pos} {set end [expr {$pos+1}]}

	    set AnnotSV_ID [settingOfTheAnnotSVID "${chrom}_${pos}_${end}_${svtype}" "$ref" "$alt"]

	    # chrUn_KN707671v1_decoy 
	    if {[lsearch -exact "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M MT" $chrom] eq -1} {
		WriteTextInFile "$AnnotSV_ID: chromosome \"$chrom\" unknown" $unannotatedOutputFile
		continue
	    }

	    if {$end eq ""} {
		WriteTextInFile "$AnnotSV_ID: END of the SV not defined" $unannotatedOutputFile
		continue
	    }

	    if {$svtype eq "CNV" || $svtype eq ""} {set svtype $alt}
	
	    
	    if {$g_AnnotSV(SVinputInfo)} {
		lappend L_TextToWrite "$chrom\t$pos\t$end\t$svtype\t[join [lrange $Ls 2 end] "\t"]"
	    } else {
		lappend L_TextToWrite "$chrom\t$pos\t$end\t$svtype\t$ref\t$alt\t[join [lrange $Ls 8 end] "\t"]"
	    }
	    # Definition of g_SVLEN:
	    # (If not defined here, the variant length can be calculated in AnnotSV-write.tcl for some type of SV)
	    if {$svlen ne ""} {
		set g_SVLEN($AnnotSV_ID) $svlen
	    }
	}
	close $f
	
	WriteTextInFile [join $L_TextToWrite "\n"] $SV_BEDfile

	if {$VCFheader eq ""} { ; # No header in the VCF input file
	    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG00096
	    if {$g_AnnotSV(SVinputInfo)} {
		set length [llength [lrange $Ls 9 end]]
		set VCFheader "SV type\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t[join [lrepeat $length ""] "\t"]"
	    } else {
		set length [llength [lrange $Ls 8 end]]
		set VCFheader "SV type\tREF\tALT\t[join [lrepeat $length ""] "\t"]"
	    }
	}

    }

    # Bedfile should be sorted and should not have "chr" in the first column
    # -> This treatment will be done in the 'refGeneAnnotation' proc

    ## Bedfile: no SV to annotate if it is an empty file, or a file with only SNV/indel
    if {[isAnEmptyFile $SV_BEDfile]} {
	puts "############################################################################"
	puts "No SV to annotate in the SVinputFile - Exit without error."
	puts "############################################################################"
	exit 0
    }

    return "$SV_BEDfile"
}
