############################################################################################################
# AnnotSV 3.0.9                                                                                            #
#                                                                                                          #
# AnnotSV: An integrated tool for Structural Variations annotation and ranking                             #
#                                                                                                          #
# Copyright (C) 2017-2021 Veronique Geoffroy (veronique.geoffroy@inserm.fr)                                #
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
# Search for the place of a number into a big ordered number list. 
# Return the indice after which we can range the number
# Return -1 if number < [lindex $ListeOrdonnee 0]
##################################################################
proc DichotomySearch {number OrderedNumberList} {
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


##################################################################
# Prepare 3 new annotation columns:
# - Number of homozygous SNV/indel detected in the sample and located in the SV
# - Number of heterozygous SNV/indel detected in the sample and  located in the SV
# - Total number of calls detected in all samples and located in the SV
##################################################################
proc VCFannotation {SVchrom SVstart SVend SVtype} {

    global g_AnnotSV
    global lPos
    global L_allSamples

    ## VCFs parsing is done only 1 time (After, g_AnnotSV(vcfParsing) is set to "done")
    if {![info exists g_AnnotSV(vcfParsing)]} {

	set g_AnnotSV(vcfParsing) "done"
	
	# parsing of $g_AnnotSV(snvIndelFiles): creation of lPos($SVchrom,htz,sample), lPos($SVchrom,hom,sample) and lPos($SVchrom,tot,cohort)
	puts "\n\n...parsing of snvIndelFiles for \"$g_AnnotSV(snvIndelSamples)\" ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])" 

	# "eval glob" accept regular expression ("*.vcf) as well as a list of files ("sample1.vcf sample2.vcf.gz"):
	set L_allSamples {}
	foreach vcfF [eval glob -nocomplain $g_AnnotSV(snvIndelFiles)] {
	    
	    # Split multiallelic sites into multiple rows
	    puts "\t...split multiallelic sites into multiple rows for $vcfF ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
            set tmpVCF1 "$g_AnnotSV(outputDir)/[file tail $vcfF].[clock seconds].separate.vcf"
            set tmpVCFgz1 "$g_AnnotSV(outputDir)/[file tail $vcfF].[clock seconds].separate.vcf.gz"
	    catch {eval exec $g_AnnotSV(bcftools) norm -m -both $vcfF > $tmpVCF1} Message
	    if {[file size $tmpVCF1] eq 0} {
		# we continue AnnotSV without splitting the multiallelic sites !
		puts "\t   -- VCFannotation --"
		puts "\t   $g_AnnotSV(bcftools) norm -m -both $vcfF > $tmpVCF1"
		puts "\t   $tmpVCF1: file empty."
		foreach line [split $Message "\n"] {
		    if {[regexp -nocase "error|fail" $line]} {puts "\t   $line"}
		}
		puts "\t   => No multiallelic treatment done."
		if {[regexp ".gz$" $vcfF]} {
		    file copy -force $vcfF $tmpVCFgz1 
		} else {
		    catch {eval exec cat $vcfF | gzip > $tmpVCFgz1} Message 
		}
	    } elseif {[regexp -nocase "error|fail" $Message]} {
		# we continue AnnotSV without splitting the multiallelic sites !
		puts "\t   -- VCFannotation --"
		puts "\t   $g_AnnotSV(bcftools) norm -m -both $vcfF > $tmpVCF1"
		foreach line [split $Message "\n"] {
		    if {[regexp -nocase "error|fail" $line]} {puts "\t   $line"}
		}
		puts "\t   => No multiallelic treatment done."
		if {[regexp ".gz$" $vcfF]} {
		    file copy -force $vcfF $tmpVCFgz1 
		} else {
		    catch {eval exec cat $vcfF | gzip > $tmpVCFgz1} Message 
		}		
	    } else {
		catch {eval exec gzip $tmpVCF1} Message
	    }
	    file delete -force $tmpVCF1
	    
	    # Parsing of $vcfF
	    puts "\t...parsing of $vcfF ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])" 
	    set iIntersect 0
	    set iSV 0
	    set iLoaded 0
	    set iNotPASS 0
	    set iNotGT 0

	    # VCF should not contain "chr" prefix
	    regsub ".vcf.gz$" $tmpVCFgz1 ".withoutChr.vcf.gz" tmpVCFgz2
	    catch {exec gunzip -c $tmpVCFgz1 | grep -c "^chr"} Message
	    if {[string index $Message 0] ne 0} {
		exec gunzip -c $tmpVCFgz1 | sed "s/^chr//" | gzip > $tmpVCFgz2
	    } else {
		set tmpVCFgz2 "$tmpVCFgz1"
	    }

	    # Header for the intersection with the bedfile
	    regsub ".vcf.gz$" $tmpVCFgz2 ".intersect.vcf.gz" tmpVCFgz3
	    regsub ".gz$" $tmpVCFgz3 "" tmpVCF3
	    set f [open "| gzip -cd $tmpVCFgz2"]	 
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
	    lappend L_allSamples {*}[lrange $L_tmpVCF_Header 9 end]	    
	    WriteTextInFile "$L_tmpVCF_Header" $tmpVCF3
	    if {[catch {exec gzip $tmpVCF3} Message]} { ; #creation of $tmpVCFgz3
		puts "-- VCFannotation --"
		puts "gzip $tmpVCF3"
		puts $Message
	    }

	    # Intersect with the bedfile
	    if {[catch {exec gunzip -c $tmpVCFgz2 | $g_AnnotSV(bedtools) intersect -a stdin -b $g_AnnotSV(bedFile) > $tmpVCF3} Message]} {
		if {[regexp -nocase "error" $Message]} {
		    WriteTextInFile "$Message" $tmpVCF3.error
		    puts "-- VCFannotation --"
		    puts "gunzip -c $tmpVCFgz2 | $g_AnnotSV(bedtools) intersect -a stdin -b $g_AnnotSV(bedFile) > $tmpVCF3"
		    puts "Error: look at $tmpVCF3.error"
		}
	    }

	    # Nettoyage
	    file delete -force $tmpVCFgz1
	    file delete -force $tmpVCFgz2

	    if {[file size $tmpVCF3] eq 0} {
		# no intersection between the bed and the VCF, no SNV/Indel to load in memory
		puts "\t\t-> no variant in the intersection"
		file delete -force $tmpVCF3
		file delete -force $tmpVCFgz3
		continue
	    } else {
		if {[catch {exec cat $tmpVCF3 | gzip >> $tmpVCFgz3} Message]} {
		    puts "-- VCFannotation --"
		    puts "cat $tmpVCF3 | gzip >> $tmpVCFgz3"
		    puts $Message
		}
		file delete -force $tmpVCF3	
	    }
	    
	    # load SNV/Indel in memory
	    set f [open "| gzip -cd $tmpVCFgz3"]	  
	    set L_samples1VCF {}
	    while {![eof $f]} {
		set L [gets $f]
		set Ls [split $L "\t"]

		if {[string range $L 0 5] eq "#CHROM"} {
		    set L_samples1VCF [lrange $Ls 9 end]
		    set k 9
		    foreach sample [lrange $Ls 9 end] {
			set i_sample($sample) $k
			incr k
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
		
		foreach sample $L_samples1VCF {
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
		if {[info exists i_samples]} {
		    unset i_samples
		}
	    }
	    close $f
	    file delete -force $tmpVCFgz3
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
    set lPosTotInTheCohort {}
    foreach sample $L_allSamples {
	# Count for hom variants
	set countHom($sample) 0
	if {[info exists lPos($SVchrom,hom,$sample)]} {
	    set i_first [DichotomySearch $SVstart $lPos($SVchrom,hom,$sample)]
	    set i_last [DichotomySearch $SVend $lPos($SVchrom,hom,$sample)]
	    if {$SVstart ne [lindex $lPos($SVchrom,hom,$sample) $i_first]} {set i_first [expr {$i_first+1}]}
	    set L_posHom [lsort -unique [lrange $lPos($SVchrom,hom,$sample) $i_first $i_last]]
	    set countHom($sample) [llength $L_posHom]
	    lappend lPosTotInTheCohort {*}$L_posHom
	} 
	# Count for htz variants
	set countHtz($sample) 0
	if {[info exists lPos($SVchrom,htz,$sample)]} {
	    set i_first [DichotomySearch $SVstart $lPos($SVchrom,htz,$sample)]
	    set i_last [DichotomySearch $SVend $lPos($SVchrom,htz,$sample)]
	    if {$SVstart ne [lindex $lPos($SVchrom,htz,$sample) $i_first]} {set i_first [expr {$i_first+1}]}
	    set L_posHtz [lsort -unique [lrange $lPos($SVchrom,htz,$sample) $i_first $i_last]]
	    set countHtz($sample) [llength $L_posHtz]
	    lappend lPosTotInTheCohort {*}$L_posHtz
	}
    }
    set lPosTotInTheCohort [lsort -unique $lPosTotInTheCohort]
    set totalCount [llength $lPosTotInTheCohort]


    if {[regexp "DEL" [normalizeSVtype $SVtype]]} {
	if {$totalCount>50} {
	    set textToReturn ""
	    foreach sample $g_AnnotSV(snvIndelSamples) {
		# Count for hom variants (Hom + WT hom)
		set CountAllHom [expr {$totalCount-$countHtz($sample)}]
		# htz/AllHom
		if {$CountAllHom eq 0} {
		    set HtzHomRatio "NA"
		} else {
		    set HtzHomRatio [format "%.4f" [expr {$countHtz($sample)*1.0/$CountAllHom}]]
		}
		# htz(sample)/total(cohort)
		if {$totalCount eq 0} {
		    set HtzTotRatio "NA"
		} else {
		    set HtzTotRatio [format "%.4f" [expr {$countHtz($sample)*1.0/$totalCount}]]
		}
		append textToReturn "$countHom($sample)\t$countHtz($sample)\t$HtzHomRatio\t$HtzTotRatio\t"
	    }
	    append textToReturn "$totalCount"
	} else {
	    foreach sample $g_AnnotSV(snvIndelSamples) {
		append textToReturn "$countHom($sample)\t$countHtz($sample)\tNA\tNA\t"
	    }
	    append textToReturn "$totalCount"
	}
    } else {
	foreach sample $g_AnnotSV(snvIndelSamples) {
	    append textToReturn "NA\tNA\tNA\tNA\t"
	}
	append textToReturn "NA"
    }

    return "$textToReturn" 
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
	
	# Check if the GT feature is present for at least 1 variant
	set GTabsent 1
	
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
			if {[lsearch -exact $g_AnnotSV(outputColHeader) "Samples_ID"] eq "-1"} {
			    set VCFheader "SV_type\t[join [lrange $Ls 2 end] "\t"]"
			} else {
			    set VCFheader "SV_type\tSamples_ID\t[join [lrange $Ls 2 end] "\t"]"
			}
		    } else {
			if {[lsearch -exact $g_AnnotSV(outputColHeader) "Samples_ID"] eq "-1"} {
			    set VCFheader "SV_type\t[join [lrange $Ls 3 4] "\t"]\t[join [lrange $Ls 8 end] "\t"]"
			} else {
			    set VCFheader "SV_type\tSamples_ID\t[join [lrange $Ls 3 4] "\t"]\t[join [lrange $Ls 8 end] "\t"]"
			}
		    }
		    set L_allSamples [lrange $Ls 9 end]
		}
		continue
	    }
	    incr i
	    
	    if {$VCFheaderNotPresent} {
		set VCFheaderNotPresent 0
		puts "ERROR:\n$VCFfile: no VCF header line (prefixed with \"#CHROM\"). Check your VCF."
		puts "Exit with error\n"
		exit 2
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
	    regsub "chr" [lindex $Ls 3] "" ref ; # for alt like: "alt="G]chr17:1584563]"
	    regsub "chr" [lindex $Ls 4] "" alt 
	    regsub -all "\"" [lindex $Ls 7] "" INFOcol
            if {[regexp "SVLEN=(\[0-9-\]+)" $INFOcol match SVLEN]} {set svlen $SVLEN} else {set svlen ""}
	    if {[regexp ";END=(\[0-9\]+)" $INFOcol match END] || [regexp "^END=(\[0-9\]+)" $INFOcol match END]} {set end $END} else {set end ""}
	    if {[regexp "SVTYPE=(\[^;\]+)" $INFOcol match SVTYPE]} {set svtype $SVTYPE} else {set svtype ""}      
	    if {[regexp "CIPOS=(\[^;\]+)" $INFOcol match CIPOS]} {set cipos $CIPOS} else {set cipos ""}      
	    if {[regexp "CIEND=(\[^;\]+)" $INFOcol match CIEND]} {set ciend $CIEND} else {set ciend ""}
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

	    if {$g_AnnotSV(includeCI)} {
		# Correction of the "start" / "end" SV positions by using the confidence interval around the boundaries:
		##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
		##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
		# e.g. CIPOS=-10,10;CIEND=-10,10;  or  CIPOS=-6,5;CIEND=-6,5;  or  CIPOS=-410,410;CIEND=-410,410
		set ciposleft [lindex [split $cipos ","] 0]
		set ciendright [lindex [split $ciend ","] end]
		if {$ciposleft < 0} {
		    catch {incr pos $ciposleft}
		}
		if {$ciendright > 0} {
		    catch {incr end $ciendright}
		}
	    }

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
	
	    # Set up for the Samples_ID value
	    set i_gt [lsearch -exact [split [lindex $Ls 8] ":"] "GT"]
	    if {$i_gt eq -1} {
		# If the GT is not given, AnnotSV considers that the SV is present in all the samples
		set L_samplesid "$L_allSamples"
	    } else {
		set GTabsent 0
		set L_samplesid {}
		set isample 0
		foreach sampleValue [lrange $Ls 9 end] {
		    set gt [lindex [split $sampleValue ":"] $i_gt]
		    if {[regexp "1|2|3|4|5|6|7|8|9" $gt]} {lappend L_samplesid [lindex $L_allSamples $isample]}
		    incr isample
		}
	    }

	    # Text to write
	    if {$g_AnnotSV(SVinputInfo)} {
		if {[lsearch -exact $g_AnnotSV(outputColHeader) "Samples_ID"] eq "-1"} {
		    lappend L_TextToWrite "$chrom\t$pos\t$end\t$svtype\t[join [lrange $Ls 2 end] "\t"]"
		} else {
		    lappend L_TextToWrite "$chrom\t$pos\t$end\t$svtype\t[join $L_samplesid ","]\t[join [lrange $Ls 2 end] "\t"]"
		}
	    } else {
		if {[lsearch -exact $g_AnnotSV(outputColHeader) "Samples_ID"] eq "-1"} {
		    lappend L_TextToWrite "$chrom\t$pos\t$end\t$svtype\t$ref\t$alt\t[join [lrange $Ls 8 end] "\t"]"
		} else {
		    lappend L_TextToWrite "$chrom\t$pos\t$end\t$svtype\t[join $L_samplesid ","]\t$ref\t$alt\t[join [lrange $Ls 8 end] "\t"]"
		}
	    }
	    
	    # Definition of g_SVLEN:
	    # (If not defined here, the variant length can be calculated in AnnotSV-write.tcl for some type of SV)
	    if {$svlen ne ""} {
		set g_SVLEN($AnnotSV_ID) $svlen
	    }
	}
	if {![regexp ".gz$" $VCFfile]} {close $f}
	
	WriteTextInFile [join $L_TextToWrite "\n"] $SV_BEDfile

	if {$g_AnnotSV(candidateSnvIndelFiles) ne "" && $GTabsent} {
	    puts "...parsing of $SV_VCFfiles"
	    puts "\t...WARNING: The SV genotype is not indicated in the FORMAT column under the “GT” field"
	    puts "\t            => Compound heterozygosity analysis won't be processed!\n"
	    set g_AnnotSV(candidateSnvIndelFiles)   ""
	    set g_AnnotSV(candidateSnvIndelSamples) ""
	}

	if {$VCFheader eq ""} { ; # No header in the VCF input file
	    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG00096
	    if {$g_AnnotSV(SVinputInfo)} {
		set length [llength [lrange $Ls 9 end]]
		set VCFheader "SV_type\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t[join [lrepeat $length ""] "\t"]"
	    } else {
		set length [llength [lrange $Ls 8 end]]
		set VCFheader "SV_type\tREF\tALT\t[join [lrepeat $length ""] "\t"]"
	    }
	}

    }

    # Bedfile should be sorted and should not have "chr" in the first column
    # -> This treatment will be done in the 'genesAnnotation' proc

    ## Bedfile: no SV to annotate if it is an empty file, or a file with only SNV/indel
    if {[isAnEmptyFile $SV_BEDfile]} {
	puts "############################################################################"
	puts "No SV to annotate in the SVinputFile - Exit without error."
	puts "############################################################################"
	exit 0
    }

    return "$SV_BEDfile"
}
