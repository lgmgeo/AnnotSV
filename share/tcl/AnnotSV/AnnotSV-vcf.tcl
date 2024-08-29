############################################################################################################
# AnnotSV 3.4.3                                                                                            #
#                                                                                                          #
# AnnotSV: An integrated tool for Structural Variations annotation and ranking                             #
#                                                                                                          #
# Copyright (C) 2017-2024 Veronique Geoffroy (veronique.geoffroy@inserm.fr)                                #
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
# Prepare 5 new annotation columns:
# - Number of homozygous SNV/indel detected in the sample and located in the SV --- $countHom($sample)
# - Number of heterozygous SNV/indel detected in the sample and  located in the SV --- $countHtz($sample)
# - Ratio of htz / Hom ---  $HtzHomRatio
# - Ratio of htz / Total --- $HtzTotRatio
# - Total number of calls detected in all samples and located in the SV --- $totalCount
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
        
        # "eval glob" accept regular expression ("*.vcf") as well as a list of files ("sample1.vcf sample2.vcf.gz"):
        set L_allSamples {}
        foreach vcfF [eval glob -nocomplain $g_AnnotSV(snvIndelFiles)] {
            
            # Split multiallelic sites into multiple rows
            puts "\t...split multiallelic sites into multiple rows for $vcfF ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
            set tmpVCF1 "$g_AnnotSV(outputDir)/[file tail $vcfF].[clock seconds].separate.vcf"
            set tmpVCFgz1 "$g_AnnotSV(outputDir)/[file tail $vcfF].[clock seconds].separate.vcf.gz"
            catch {eval exec $g_AnnotSV(bcftools) norm -m -both $vcfF > $tmpVCF1} Message
            if {[file size $tmpVCF1] eq 0} {
                # we continue AnnotSV without splitting the multiallelic sites!
                puts "\t   -- VCFannotation --"
                puts "\t   $g_AnnotSV(bcftools) norm -m -both $vcfF > $tmpVCF1"
                puts "\t   $tmpVCF1: file empty."
                foreach line [split $Message "\n"] {
                    if {[regexp -nocase "error|fail" $line]} {puts "\t   $line"}
                }
                puts "\t   => No multiallelic treatment done."
                if {[regexp -nocase ".gz$" $vcfF]} {
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
                if {[regexp -nocase ".gz$" $vcfF]} {
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
            
            # Header for the intersection with the SV bedfile
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
                puts "\t\t-> $iIntersect SNV/indel located in all the SV (potentially with redundancy depending on SV overlap)"
            }
            if {$iLoaded} {
                puts "\t\t-> $iLoaded SNV/indel loaded"
            }
            if {$iSV} {
                puts "\t\t-> $iSV SV excluded (considering only SNV/indel from the VCF)"
            }
            if {$iNotPASS} {
                puts "\t\t-> $iNotPASS SNV/indel excluded beacause of the FILTER value not equal to \"PASS\""
            }
            if {$iNotGT} {
                puts "\t\t-> $iNotGT SNV/indel excluded because of the absence of GT information"
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
    global g_ID
    global g_deb
    
    
    ## BED: 0-based, non-inclusive-end [startFromBED, endFromBED[
    ## VCF: 1-based, inclusive-end [startFromBED+1, endFromBED]
    
    ## It is to notice in the VCF specification:
    ## - For imprecise structural variants (i.e. symbolic allele with bracketed notation; e.g. <DUP>):
    ##   END = POS + length of REF allele
    ## - For precise structural variants (e.g. ref="G" and alt="ACTGCTA):
    ##   END = POS + length of REF allele - 1
    
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
        puts "...VCF to BED ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
        set L_TextToWrite {}
        # TextToWrite_rescue($SV_ID)
        set L_squBrack_SV_ID_Written {}
        
        if {[regexp -nocase ".gz$" $VCFfile]} {
            set f [open "| gzip -cd $VCFfile"]
        } else {
            set f [open "$VCFfile"]
        }
        
        set i 0
        set VCFheaderNotPresent 1
        
        # Check if the GT feature is present for at least 1 variant
        #
        # GT info:
        # - 0/1 and 1/0 functionally mean the same thing (that the individual is a heterozygote) - since the genotype is unphased, the alleles aren't ordered.
        #   The / symbol tells you the genotype is unphased.
        #   The convention is to write the GT field in ascending order, so 0/1 rather than 1/0
        # - 0|1 and 1|0 do mean something different, since the pipe symbol "|" tell us the order of the alleles matters
        # - ”.” must be specified for each missing allele in the GT field (for example ./.) (=> if a call cannot be made for a sample at a given locus)
        #
        # In the code:
        # - Sample ID with GT "./." or ".|." (unknown GT) are reported in the "Samples_ID" output column
        set GTabsent 1; #To check if the “GT” field is indicated in the FORMAT column
        
        set VCFlineNumber 0
        set nUnknownGT 0; # Number of Samples with unknown GT for an SV
        
        while {![eof $f]} {
            
            if {$i eq "500000"} {
                WriteTextInFile [join $L_TextToWrite "\n"] $SV_BEDfile
                set L_TextToWrite {}
                set i 0
            }
            
            set L [gets $f]
            set Ls [split $L "\t"]
            incr VCFlineNumber
            if {[string index $L 0] eq "#" || $L eq ""} {
                if {[regexp "^#CHROM" $L]} {
                    # Check if this header is complete
                    foreach classicValue {#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT} {
                        if {[lsearch -exact $Ls $classicValue] eq -1} {
                            puts "VCF header file badly formatted:"
                            puts $L
                            puts "($classicValue not present)"
                            puts "Exit with error"
                            exit 2
                        }
                    }
                    set VCFheaderNotPresent 0
                    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG00096
                    if {$g_AnnotSV(SVinputInfo)} {
                        set VCFheader "SV_type\tSamples_ID\t[join [lrange $Ls 2 end] "\t"]"
                    } else {
                        set VCFheader "SV_type\tSamples_ID\t[join [lrange $Ls 3 4] "\t"]\t[join [lrange $Ls 8 end] "\t"]"
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
            # - Type3: squared-bracketed SV notation: alt="G]17:1584563]" or alt="G]chr17:1584563]" (length > 50bp)
            #          (developed with the assistance and guidance of Rodrigo Martin, BSC, Spain)
            set BNDrescue 0 ;# The ALT of a rescue BND is not square-bracketed anymore
            set TRArescue 0
            set SquareBracketedSV 0
            set chrom [lindex $Ls 0]
            regsub -nocase "chr" $chrom "" chrom
            set posVCF [lindex $Ls 1]
            set pos [expr {$posVCF-1}]; # VCF: 1-based ==> BED: 0-based
            regsub "chr" [lindex $Ls 3] "" ref ; # for alt like: "alt="G]chr17:1584563]"
            regsub "chr" [lindex $Ls 4] "" alt
            set altVCF $alt ;# To use with the proc "settingOfTheAnnotSVID" (because alt is modified if it is a square-bracketed notation)
            set refVCF $ref
            regsub "CHR2=chr" $Ls "CHR2=" Ls
            regsub -all "\"" [lindex $Ls 7] "" INFOcol
            set svlen ""; set end ""; set svtype ""; set cipos ""; set ciend ""; set chr2 ""; set cipos95 ""; set ciend95 ""
            if {[regexp "SVLEN=(\[0-9-\]+)" $INFOcol match SVLEN]} {set svlen $SVLEN}
            if {[regexp ";END=(\[0-9\]+)" $INFOcol match END] || [regexp "^END=(\[0-9\]+)" $INFOcol match END]} {set end $END}
            if {[regexp "SVTYPE=(\[^;\]+)" $INFOcol match SVTYPE]} {set svtype $SVTYPE}
            if {[regexp "CIPOS=(\[^;\]+)" $INFOcol match CIPOS]} {set cipos $CIPOS}
            if {[regexp "CIEND=(\[^;\]+)" $INFOcol match CIEND]} {set ciend $CIEND}
            if {[regexp "CHR2=(\[^;\]+)" $INFOcol match CHR2]} {set chr2 $CHR2} ;# Translocation
            if {[regexp "CIPOS95=(\[^;\]+)" $INFOcol match CIPOS95]} {set cipos95 $CIPOS95}
            if {[regexp "CIEND95=(\[^;\]+)" $INFOcol match CIEND95]} {set ciend95 $CIEND95}
            if {[regexp "^<" $alt]} {
                # Type2: angle-bracketed notation <...>
                
                if {$end eq ""} {
                    # INS:ME (LINE1, ALU or SVA)
                    set end $posVCF
                }
                
                if {[regexp "<TRA>" $alt]} {
                    # Example of a translocation VCF format line:
                    # chr1   63705386   N    .     <TRA>   51      PASS    "PRECISE;CT=5to5;CIPOS=-10,10;CIEND=-10,10;SOMATIC;SVTYPE=TRA;CHR2=chrX;END=478444;SVLEN=0"        GT      ./.
                    # POS = 63705386 --> sur chrom1
                    # END = 478444   --> sur chromX
                    ## => annotation should consider only this breakpoint on chr1 (not the BND on chrX because the SVINFO features are specific to this BND on chr1, not to the other one on chrX)
                    
                    #2       321682  tra_a   T       <TRA>   .       PASS    SVTYPE=BND;CHR2=13;END=123456  GT      0/0     0/1
                    set end $posVCF
                    set svtype "TRA"
                    set svlen 0
                } else {
                    # First, we choose the SVtype in the ALT column.
                    # Second, we choose the SVTYPE in the INFO column.
                    if {[regexp "^<(.+)>$" $alt match svtype2]} {
                        set svtype $svtype2
                    }
                    set svtype [normalizeSVtype $svtype] ;# DEL or DUP or INS or INV or None
                }
                
                if {$svlen eq "" && $end ne ""} {
                    if {$svtype eq "DEL"} {
                        set svlen [expr {$posVCF-$end}]
                    } elseif {[regexp "DUP|INV" $svtype]} {
                        set svlen [expr {$end-$posVCF}]
                    }
                }
                
                # Check the SVminSize for DUP, DEL and INS
                if {$svlen ne "" && ($svtype eq "DUP" || $svtype eq "DEL" || $svtype eq "INS")} {
                    if {[expr {abs($svlen)}]<$g_AnnotSV(SVminSize)} {
                        # it is a small variant
                        WriteTextInFile "${chrom}_${posVCF}_${end}_${svtype}_${ref}_${altVCF}: variantLength ([expr {abs($svlen)}]) < SVminSize ($g_AnnotSV(SVminSize)) (line $VCFlineNumber)" $unannotatedOutputFile
                        continue
                    }
                }
                
            } elseif {[regexp "(\[acgtnACGTN\]*)(\\\[|\\\])(\[^:\]+):(\[0-9\]+)(\\\[|\\\])(\[acgtnACGTN\]*)" $alt match baseLeft bracketLeft inBracketChrom inBracketStart bracketRight baseRight]} {
                # Type3: square-bracketed notation ]...]
                # => Developed with the assistance and guidance of Rodrigo Martin, BSC, Spain
                
                # With this bracketed notation, we have 1 line in the VCF for each breakend (=> 2 lines)
                # => It corresponds to only 1 line in the BED file, except for translocation (where both breakends are annotated)
                #
                # Rules applied here:
                # - DUP, DEL, INS and INV: We look first at the lines whith pos_POS < pos_ALT.
                #                          The other lines (pos_POS > pos_ALT) are only used for BND rescue
                # - TRA: We look at the 2 breakpoints lines.
                
                # The ALT of a recovered BND is not square-bracketed anymore but angle-bracketed
                # The REF is set to "N"
                set SquareBracketedSV 1
                
                if {$chrom != $inBracketChrom} {
                    # TRA (chrom_#CHROM != chrom_ALT)
                    #################################
                    # Example:
                    # 17      198982  trn_no_mateid_a A       A]2:321681]	=> Line analysed by AnnotSV
                    # 2       321681  trn_no_mateid_b G       G]17:198982]	=> Line analysed by AnnotSV
                    # (INFO: Can not have a translocation on the same chrom. A translocation on the same chrom is always represented with another thing: DUP, INV or DEL)
                    #        Can not rescue the orientation of the square-brackets from one BND to the paired BND ("]" or "[") => for the rescue: ALT = <TRA>
                    set ref "N"
                    set alt "<TRA>"
                    set end $posVCF
                    set svtype "TRA"
                    set svlen 0
                    
                    set TRArescue 1
                    set altVCF2 "<TRA>"
                    set chrom2 $inBracketChrom
                    set pos2 [expr {$inBracketStart-1}]; # VCF: 1-based ==> BED: 0-based
                    set posVCF2 $inBracketStart
                    set end2 $inBracketStart
                    set INFOcol2 "$INFOcol;BNDrescue"
                    
                    # #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample1  sample2
                    set L2 "$chrom2\t$inBracketStart\t[lindex $Ls 2]\t$ref\t$altVCF2\t[join [lrange $Ls 5 6] "\t"]\t$INFOcol2\t[join [lrange $Ls 8 end] "\t"]"
                    set L2s [split $L2 "\t"]
                    regsub -all "\\\\" $L2s "" L2s
                    
                } else {
                    if {$posVCF > $inBracketStart} {
                        # Lines only used for BND rescue
                        
                        # DUP, DEL, INS and INV: We first look at the lines whith pos_POS < pos_ALT
                        # => Else, search for the reciprocal breakend in case of rescue needed
                        # With one breakend, you can always infer the other. Indeed, the only thing that could be different from the mate breakend is
                        # its CIPOS INFO field (but it should be provided in the CIEND field of the other breakend).
                        # Regarding GT, both breakends must have the same GT because they represent the same thing.
                        # The CIEND/CIPOS relationship is that you can use the CIEND info from the breakend you already have to set the CIPOS field of the new "inferred" breakend.
                        
                        set BNDrescue 1; # The other BND from this BND pair can be recovered here
                        
                        set remember $posVCF
                        set posVCF $inBracketStart
                        set inBracketStart $remember
                        set pos [expr {$posVCF-1}]; # VCF: 1-based ==> BED: 0-based
                        
                        # Modif in $INFOcol: Correction of the CIPOS, CIEND...
                        #                    Add "BNDrescue"
                        set remember $cipos
                        set cipos $ciend
                        set ciend $remember
                        regsub "CIPOS=\[^;\]+" $INFOcol "CIPOS=$cipos" INFOcol
                        regsub "CIEND=\[^;\]+" $INFOcol "CIEND=$ciend" INFOcol
                        set remember $cipos95
                        set cipos95 $ciend95
                        set ciend95 $remember
                        regsub "CIPOS95=\[^;\]+" $INFOcol "CIPOS95=$cipos95" INFOcol
                        regsub "CIEND95=\[^;\]+" $INFOcol "CIEND95=$ciend95" INFOcol
                        append INFOcol ";BNDrescue"
                        
                        set ref "N"
                        if {[string length $baseLeft] > 1} {
                            set baseLeft "[string range $baseLeft 1 end]N"
                        } elseif {[string length $baseRight] > 1} {
                            set baseRight "N[string range $baseRight 0 end-1]"
                        }
                        set remember $baseLeft
                        set baseLeft $baseRight
                        set baseRight $remember
                        if {$bracketLeft == "\["} {
                            set bracketLeft "\]"
                            set bracketRight "\]"
                        } else {
                            set bracketLeft "\["
                            set bracketRight "\["
                        }
                        set alt "$baseLeft$bracketLeft${inBracketChrom}:$inBracketStart$bracketRight$baseRight"
                        set altVCF $alt
                        
                        # #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample1  sample2
                        set L "$chrom\t$posVCF\t[lindex $Ls 2]\t$ref\t$altVCF\t[join [lrange $Ls 5 6] "\t"]\t$INFOcol\t[join [lrange $Ls 8 end] "\t"]"
                        set Ls [split $L "\t"]
                    }
                    
                    if {[string length $baseLeft] > 1 || [string length $baseRight] > 1} {
                        # INS ("first mapped base" is followed by the inserted sequence)
                        ################################################################
                        # Example:
                        # 13      53040041        ins_by_gridss   T       TATATATATACACAC[13:53040042[	=> Line analysed by AnnotSV
                        # 13      53040042        ins_by_gridss   A       ]13:53040041]ATATATATACACACA	=> Line NOT reported by AnnotSV (if the previous BND has already been analysed)
                        set svlen [expr {[string length $baseLeft]-1}]
                        set svtype "INS"
                        set end $posVCF
                        set alt "<INS>"
                        if {$svlen<$g_AnnotSV(SVminSize)} {
                            # it is an indel
                            WriteTextInFile "${chrom}_${posVCF}_${end}_${svtype}_${ref}_${altVCF}: variantLength ($svlen) < SVminSize ($g_AnnotSV(SVminSize)) (line $VCFlineNumber)" $unannotatedOutputFile
                            continue
                        }
                    } elseif {$bracketRight eq "]" && $baseRight ne ""} {
                        # DUP ("first mapped base" is NOT contained in the bracket: "N[" or "]N"; REF is after the brackets
                        ###################################################################################################
                        # Example:
                        # 2       3000    breakend_dup_a  T       ]2:5000]T  	=> Line analysed by AnnotSV
                        # 2       5000    breakend_dup_b  T       T[2:3000[ 	=> Line NOT reported by AnnotSV (if the previous BND has already been analysed)
                        set svlen [expr {$inBracketStart-$posVCF}]
                        set svtype "DUP"
                        set end $inBracketStart
                        set alt "<DUP>"
                        if {$svlen<$g_AnnotSV(SVminSize)} {
                            # it is an indel
                            WriteTextInFile "${chrom}_${posVCF}_${end}_${svtype}_${ref}_${altVCF}: variantLength ($svlen) < SVminSize ($g_AnnotSV(SVminSize)) (line $VCFlineNumber)" $unannotatedOutputFile
                            continue
                        }
                    } elseif {($bracketRight eq "]" && $baseLeft ne "") || ($bracketRight eq "\[" && $baseRight ne "")} {
                        # INV ("first mapped base" is contained in the bracket: "N]" or "[N")
                        #####################################################################
                        # Example 1:
                        # 3       2999    breakend_inv_1_a        T       T]3:5000]  	=> Line analysed by AnnotSV
                        # 3       5000    breakend_inv_1_b        T       [3:2999[T  	=> Line NOT reported by AnnotSV (if the previous BND has already been analysed)
                        # Example 2:
                        # 3       3000    breakend_inv_2_a        T       [3:5001[T  	=> Line analysed by AnnotSV
                        # 3       5001    breakend_inv_2_b        T       T]3:3000] 	=> Line NOT reported by AnnotSV (if the previous BND has already been analysed)
                        # Example 3:
                        # 3       3000    breakend_inv_3_a        T       [3:5001[T     => Line analysed by AnnotSV
                        # 3       5001    breakend_inv_3_b        T       [3:3000[T     => Line NOT reported by AnnotSV (if the previous BND has already been analysed)
                        set svlen [expr {$inBracketStart-$posVCF}]
                        set svtype "INV"
                        set end $inBracketStart
                        set alt "<INV>"
                    } elseif {$bracketRight eq "\[" && $baseLeft ne ""} {
                        # DEL ("first mapped base" is NOT contained in the bracket: "N[" or "]N"; "first mapped base" is before the brackets
                        ####################################################################################################################
                        # Example:
                        # 12      3000    breakend_del_1_a        T       T[12:5000[ 	=> Line analysed by AnnotSV
                        # 12      5000    breakend_del_1_b        T       ]12:3000]T  	=> Line NOT reported by AnnotSV (if the previous BND has already been analysed)
                        set svlen [expr {$posVCF-$inBracketStart}];# negative value for del
                        set svtype "DEL"
                        set end $inBracketStart
                        set alt "<DEL>"
                        if {[expr {abs($svlen)}]<$g_AnnotSV(SVminSize)} {
                            # it is an indel
                            WriteTextInFile "${chrom}_${posVCF}_${end}_${svtype}_${ref}_${altVCF}: variantLength ([expr {abs($svlen)}]) < SVminSize ($g_AnnotSV(SVminSize)) (line $VCFlineNumber)" $unannotatedOutputFile
                            continue
                        }
                    }
                }
            } elseif {[regexp -nocase "^\[acgtnACGTN.*\]+$" $ref$alt]} {
                # Type1
                regsub -all "\[*.\]" $ref "" refbis ;# cf GRIDSS comment, just below
                regsub -all "\[*.\]" $alt "" altbis
                
                set variantLengthType1 [expr {[string length $altbis]-[string length $refbis]}]
                if {[expr {abs($variantLengthType1)}] < $g_AnnotSV(SVminSize)} {
                    # it is a small insertion or deletion (indel)
                    WriteTextInFile "${chrom}_${posVCF}_${end}_${svtype}_${ref}_${altVCF}: variantLength ([expr {abs($variantLengthType1)}]) < SVminSize ($g_AnnotSV(SVminSize)) (line $VCFlineNumber)" $unannotatedOutputFile
                    continue
                }
                
                if {$variantLengthType1>0} {
                    # insertion
                    if {$end eq ""} {set end $posVCF}
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
                WriteTextInFile "${chrom}_${posVCF}_${end}_${svtype}_${ref}_${altVCF}: not a known SV format (line $VCFlineNumber)" $unannotatedOutputFile
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
                    if {$pos <=0} {set pos 1}
                }
                if {$ciendright > 0} {
                    catch {incr end $ciendright}
                }
            }
            
            # chrUn_KN707671v1_decoy
            set L_orderedChrom {1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M MT}
            if {[lsearch -exact "$L_orderedChrom" $chrom] eq -1} {
                WriteTextInFile "${chrom}_${posVCF}_${end}_${svtype}_${ref}_${altVCF}: chromosome \"$chrom\" unknown (line $VCFlineNumber)" $unannotatedOutputFile
                continue
            }
            # END of the SV not defined
            if {$end eq ""} {
                WriteTextInFile "${chrom}_${posVCF}_${end}_${svtype}_${ref}_${altVCF}: END of the SV not defined (line $VCFlineNumber)" $unannotatedOutputFile
                continue
            }
            
            # ALT with square-bracketed notation (For duplication, inversion, deletion and insertion. Not for translocation)
            # If both BND (in a SV pair) are presents and ordered, we don't use the line with the 2d BND of the pair.
            if {$BNDrescue} {
                set idWithN "${chrom}_${posVCF}_${end}_${svtype}_${ref}_[replaceREFwithNinALT ${altVCF}]"
                regsub -all "\\\]|\\\[" $idWithN ";" idWithN
                if {[lsearch -exact $L_squBrack_SV_ID_Written "$idWithN"] ne -1} {
                    WriteTextInFile "$idWithN: reciprocal breakend (line $VCFlineNumber)" $unannotatedOutputFile
                    continue
                }
                if {$svtype eq "INV"} {
                    set idWithN "${chrom}_${posVCF}_${end}_${svtype}_${ref}_[replaceREFwithNinALT ${altVCF} "inv"]"
                    regsub -all "\\\]|\\\[" $idWithN ";" idWithN
                    if {[lsearch -exact $L_squBrack_SV_ID_Written "$idWithN"] ne -1} {
                        WriteTextInFile "$idWithN: reciprocal breakend (line $VCFlineNumber)" $unannotatedOutputFile
                        continue
                    }
                }
            }
            
            
            if {$SquareBracketedSV} {
                if {$svtype eq "TRA"} {
                    # If both BND (in a TRA pair) are presents: the rescue of the 2d BND is already listed in $L_squBrack_SV_ID_Written
                    set idTRArescue "${chrom2}_${posVCF2}_${end2}_${svtype}_${ref}_<TRA>"
                    if {[lsearch -exact $L_squBrack_SV_ID_Written "$idTRArescue"] ne -1} {
                        set TRArescue 0
                    }
                }
            }
            
            # Set up for the Samples_ID value
            # nUnknownGT = Number of Samples with unknown GT for an SV
            set i_gt [lsearch -exact [split [lindex $Ls 8] ":"] "GT"]
            if {$i_gt eq -1} {
                # If the GT is not given, AnnotSV considers that the SV is present in all the samples
                set L_samplesid "$L_allSamples"
            } else {
                set GTabsent 0; # => the “GT” field is indicated in the FORMAT column
                set L_samplesid {}
                set isample 0
                foreach sampleValue [lrange $Ls 9 end] {
                    set gt [lindex [split $sampleValue ":"] $i_gt]
                    if {$gt ne "0/0" && $gt ne "0\|0"} {lappend L_samplesid [lindex $L_allSamples $isample]}; # AnnotSV reports the sample_id with unknown GT ("./." and ".|.")
                    if {$gt eq "./." || $gt eq ".\|."} {incr nUnknownGT}
                    incr isample
                }
            }
            
            # Text to write + Text to possibly rescue
            #
            # AnnotSV_ID is required to set the "g_SVLEN($AnnotSV_ID)" variable
            # - For the setting of AnnotSV_ID, we use the "ref" and "alt" written in the BED file!
            # - At the end, we unset the global variables "g_ID" and "g_deb"
            
            if {$g_AnnotSV(SVinputInfo)} {
                if {$BNDrescue} {
                    # With the BND recovery, the REF and ALT fields have an undefined NT ("N")
                    # #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample1  sample2
                    set idWithN "${chrom}_${posVCF}_${end}_${svtype}_${ref}_[replaceREFwithNinALT ${altVCF}]"
                    regsub -all "\\\]|\\\[" "$idWithN" ";" idWithN
                    set TextToWrite_rescue($idWithN) "$chrom\t$pos\t$end\t$svtype\t[join $L_samplesid ","]\t[lindex $Ls 2]\t$ref\t$alt\t[join [lrange $Ls 5 end] "\t"]"
                    set VCFline_rescue($idWithN) $VCFlineNumber
                    set AnnotSV_ID [settingOfTheAnnotSVID "${chrom}_${pos}_${end}_${svtype}" "$ref" "$alt"]
                } else {
                    lappend L_TextToWrite "$chrom\t$pos\t$end\t$svtype\t[join $L_samplesid ","]\t[join [lrange $Ls 2 end] "\t"]"
                    if {$SquareBracketedSV} {
                        if {$svtype eq "TRA"} {
                            set idTRA "${chrom}_${posVCF}_${end}_${svtype}_${ref}_<TRA>"
                            lappend L_squBrack_SV_ID_Written "$idTRA"
                        } else {
                            # lappend add "\" before "]" => we remove all []
                            regsub -all "\\\]|\\\[" "${chrom}_${posVCF}_${end}_${svtype}_${refVCF}_${altVCF}" ";" idWithoutN
                            lappend L_squBrack_SV_ID_Written "$idWithoutN"
                            regsub -all "\\\]|\\\[" "${chrom}_${posVCF}_${end}_${svtype}_N_[replaceREFwithNinALT ${altVCF}]" ";" idWithN
                            lappend L_squBrack_SV_ID_Written "$idWithN"
                            if {$svtype eq "INV"} {
                                regsub -all "\\\]|\\\[" "${chrom}_${posVCF}_${end}_${svtype}_N_[replaceREFwithNinALT ${altVCF} "inv"]" ";" idWithN
                                lappend L_squBrack_SV_ID_Written "$idWithN"
                            }
                        }
                    }
                    set AnnotSV_ID [settingOfTheAnnotSVID "${chrom}_${pos}_${end}_${svtype}" "$refVCF" "$altVCF"]
                }
                if {$TRArescue} {
                    set idTRA "${chrom2}_${posVCF2}_${end2}_${svtype}_${ref}_<TRA>"
                    set TextToWrite_rescue($idTRA) "$chrom2\t$pos2\t$end2\t$svtype\t[join $L_samplesid ","]\t[lindex $L2s 2]\t$ref\t$alt\t[join [lrange $L2s 5 end] "\t"]"
                    set VCFline_rescue($idTRA) $VCFlineNumber
                    set AnnotSV_ID2 [settingOfTheAnnotSVID "${chrom2}_${pos2}_${end2}_${svtype}" "$ref" "$alt"]
                }
            } else {
                # SVinputInfo = 0
                if {$BNDrescue} {
                    set idWithN "${chrom}_${posVCF}_${end}_${svtype}_${ref}_[replaceREFwithNinALT ${altVCF}]"
                    regsub -all "\\\]|\\\[" "$idWithN" ";" idWithN
                    set TextToWrite_rescue($idWithN) "$chrom\t$pos\t$end\t$svtype\t[join $L_samplesid ","]\t$ref\t$alt\t[join [lrange $Ls 8 end] "\t"]"
                    set VCFline_rescue($idWithN) $VCFlineNumber
                    set AnnotSV_ID [settingOfTheAnnotSVID "${chrom}_${pos}_${end}_${svtype}" "$ref" "$altVCF"]
                } else {
                    lappend L_TextToWrite "$chrom\t$pos\t$end\t$svtype\t[join $L_samplesid ","]\t$refVCF\t$altVCF\t[join [lrange $Ls 8 end] "\t"]"
                    if {$SquareBracketedSV} {
                        if {$svtype eq "TRA"} {
                            set idTRA "${chrom}_${posVCF}_${end}_${svtype}_${ref}_<TRA>"
                            lappend L_squBrack_SV_ID_Written "$idTRA"
                        } else {
                            # lappend add "\" before "]" => we remove all []
                            regsub -all "\\\]|\\\[" "${chrom}_${posVCF}_${end}_${svtype}_${refVCF}_${altVCF}" ";" idWithoutN
                            lappend L_squBrack_SV_ID_Written "$idWithoutN"
                            regsub -all "\\\]|\\\[" "${chrom}_${posVCF}_${end}_${svtype}_N_[replaceREFwithNinALT ${altVCF}]" ";" idWithN
                            lappend L_squBrack_SV_ID_Written "$idWithN"
                            if {$svtype eq "INV"} {
                                regsub -all "\\\]|\\\[" "${chrom}_${posVCF}_${end}_${svtype}_N_[replaceREFwithNinALT ${altVCF} "inv"]" ";" idWithN
                                lappend L_squBrack_SV_ID_Written "$idWithN"
                            }
                        }
                    }
                    set AnnotSV_ID [settingOfTheAnnotSVID "${chrom}_${pos}_${end}_${svtype}" "$refVCF" "$altVCF"]
                }
                if {$TRArescue} {
                    set idTRA "${chrom2}_${posVCF2}_${end2}_${svtype}_${ref}_<TRA>"
                    set TextToWrite_rescue($idTRA) "$chrom2\t$pos2\t$end2\t$svtype\t[join $L_samplesid ","]\t$ref\t$alt\t[join [lrange $L2s 8 end] "\t"]"
                    set VCFline_rescue($idTRA) $VCFlineNumber
                    set AnnotSV_ID2 [settingOfTheAnnotSVID "${chrom2}_${pos2}_${end2}_${svtype}" "$ref" "$alt"]
                }
            }
            
            # Definition of g_SVLEN:
            # (If not defined here, the variant length can be calculated in AnnotSV-write.tcl for some type of SV)
            if {$svlen ne ""} {
                set g_SVLEN($AnnotSV_ID) $svlen
                if {$TRArescue} {
                    set g_SVLEN($AnnotSV_ID2) 0
                }
            }
        }
        if {![regexp -nocase ".gz$" $VCFfile]} {close $f}
        
        
        # Warning for the unknown GT
        ############################
        if {$nUnknownGT > 1} {
            puts "\t...WARNING: $nUnknownGT sample IDs with missing alleles in the GT field (./. or .\|.) have been reported in the \"Samples_ID\" output field\n"
        } elseif {$nUnknownGT eq 1} {
            puts "\t...WARNING: 1 sample ID with missing alleles in the GT field (./. or .\|.) has been reported in the \"Samples_ID\" output field\n"
        } else {puts "\n"}
        
        # Writing of the BED file
        #########################
        WriteTextInFile [join $L_TextToWrite "\n"] $SV_BEDfile
        
        # Writing of the BND rescue lines (if needed)
        #############################################
        set L_TextToWrite {}
        foreach SV_ID [array names TextToWrite_rescue] {
            if {[lsearch -exact $L_squBrack_SV_ID_Written "$SV_ID"] eq -1} {
                lappend L_TextToWrite $TextToWrite_rescue($SV_ID)
            } else {
                WriteTextInFile "$SV_ID: reciprocal breakend (line $VCFline_rescue($SV_ID))" $unannotatedOutputFile
            }
        }
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
    
    # Unset the global variables "g_ID" and "g_deb" before to set the AnnotSV_ID variables in AnnotSV-write.tcl
    unset g_ID
    unset g_deb
    
    return "$SV_BEDfile"
}

