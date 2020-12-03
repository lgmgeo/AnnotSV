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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                            #
# GNU General Public License for more details.                                                             #
#                                                                                                          #
# You should have received a copy of the GNU General Public License                                        #
# along with this program; If not, see <http://www.gnu.org/licenses/>.                                     #
############################################################################################################


# Creation / update (if some data source files are presents) of:
# - $benignLossFile_Sorted
# - $benignGainFile_Sorted
# - $benignInsFile_Sorted
# - $benignInvFile_Sorted
proc checkBenignFiles {} {

    global g_AnnotSV

    foreach genomeBuild {GRCh37 GRCh38} {
	set benignDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/SVincludedInFt/BenignSV/$genomeBuild"

	# Files to create / update
	set benignLossFile_Sorted "$benignDir/benign_Loss_SV_$genomeBuild.sorted.bed"
	set benignGainFile_Sorted "$benignDir/benign_Gain_SV_$genomeBuild.sorted.bed"
	set benignInsFile_Sorted "$benignDir/benign_Ins_SV_$genomeBuild.sorted.bed"
	set benignInvFile_Sorted "$benignDir/benign_Inv_SV_$genomeBuild.sorted.bed"
	
	# Text to write only if an update of the sources is requested (in the check* proc).
	set g_AnnotSV(benignText) "\t...creation/update of the \"benign_*_SV_$genomeBuild.sorted.bed\" files ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])\n\t   (done only once)"	    

	# Creation of benign*File_Tmp (not formatted and sorted)
	set benignLossFile_Tmp "$benignDir/benign_Loss_SV_$genomeBuild.tmp.bed"
	set benignGainFile_Tmp "$benignDir/benign_Gain_SV_$genomeBuild.tmp.bed"
	set benignInsFile_Tmp "$benignDir/benign_Ins_SV_$genomeBuild.tmp.bed"
	set benignInvFile_Tmp "$benignDir/benign_Inv_SV_$genomeBuild.tmp.bed"
	file delete -force $benignLossFile_Tmp
	file delete -force $benignGainFile_Tmp
	file delete -force $benignInsFile_Tmp
	file delete -force $benignInvFile_Tmp	
	checkGnomAD_benignFile $genomeBuild
	checkDGV_benignFile $genomeBuild
	checkDDD_benignFile $genomeBuild
	check1000g_benignFile $genomeBuild
	checkClinGenHITS_benignFile $genomeBuild
	checkClinVar_benignFile $genomeBuild
	checkIMH_benignFile $genomeBuild
	catch {unset g_AnnotSV(benignText)}

	# Creation of *.tmp.formatted.bed
	foreach SVtype {Loss Gain Ins Inv} {
	    set benignFile_Sorted "$benignDir/benign_${SVtype}_SV_$genomeBuild.sorted.bed"
	    set benignFile_Tmp "$benignDir/benign_${SVtype}_SV_$genomeBuild.tmp.bed"
	    if {[file exists $benignFile_Tmp]} {
		# Add the SV from $benign*File_Sorted 
		if {[file exists $benignFile_Sorted]} {
		    set id_file1 [open "$benignFile_Sorted" r]
		    set id_file2 [open "$benignFile_Tmp" a]
		    while { [gets $id_file1 line] >= 0 } {
			puts $id_file2 $line
		    }
		    close $id_file1
		    close $id_file2
		}
		# sort
		catch {checkBed $benignFile_Tmp $benignDir}
		file delete -force $benignFile_Tmp
	    }
	}
	
	# Creation of *.sorted.bed:
	# Intersection with very large files can cause trouble with excessive memory usage.
	# A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm. 
	foreach SVtype {Loss Gain Ins Inv} {
	    set benignFile_TmpFormatted "$benignDir/benign_${SVtype}_SV_$genomeBuild.tmp.formatted.bed"
	    if {![file exist $benignFile_TmpFormatted]} {continue}
	    set benignFile_Sorted "$benignDir/benign_${SVtype}_SV_$genomeBuild.sorted.bed"
	    set sortTmpFile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_sort.tmp.bash"
	    ReplaceTextInFile "#!/bin/bash" $sortTmpFile
	    WriteTextInFile "# The locale specified by the environment can affects the traditional sort order. We need to use native byte values." $sortTmpFile
	    WriteTextInFile "export LC_ALL=C" $sortTmpFile
	    WriteTextInFile "sort -k1,1 -k2,2n $benignFile_TmpFormatted > $benignFile_Sorted" $sortTmpFile
	    file attributes $sortTmpFile -permissions 0755
	    if {[catch {eval exec bash $sortTmpFile} Message]} {
		puts "-- checkBenignFiles --"
		puts "sort -k1,1 -k2,2n $benignFile_TmpFormatted > $benignFile_Sorted"
		puts "$Message"
		puts "Exit with error"
		exit 2
	    }
	    file delete -force $sortTmpFile 
	    file delete -force $benignFile_TmpFormatted
	}
    }
    
    # Remove some columns from g_AnnotSV(outputColHeader) if the corresponding annotation file doesn't exist
    foreach svtype {"Gain" "Loss" "Ins" "Inv"} {
	set benignFile_Sorted "$benignDir/benign_${svtype}_SV_$genomeBuild.sorted.bed"
	if {![file exists $benignFile_Sorted]} {
	    set newList {}
	    foreach e "$g_AnnotSV(outputColHeader)" {
		if {[regexp "^B_[string tolower ${svtype}]_" $e]} {continue}
		lappend newList $e
	    }
	    set g_AnnotSV(outputColHeader) $newList
	}
    }

    return
}




################################################################################################
## In all the next "check*File" proc, create and then complete the following file:
##   - Benign_Gain_SV_$genomeBuild.tmp.bed
##   - Benign_Loss_SV_$genomeBuild.tmp.bed
##   - Benign_Ins_SV_$genomeBuild.tmp.bed
##   - Benign_Inv_SV_$genomeBuild.tmp.bed
################################################################################################

proc checkClinVar_benignFile {genomeBuild} {
    
    global g_AnnotSV
    
    ## Check if ClinVar file has been downloaded 
    ############################################
    set benignDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/SVincludedInFt/BenignSV/$genomeBuild"
    set ClinVarFileDownloaded [lindex [glob -nocomplain "$benignDir/clinvar*.vcf.gz"] end]

    if {$ClinVarFileDownloaded ne ""} {
	# We have some ClinVar annotations to add in $benign*File
	if {[info exists g_AnnotSV(benignText)]} {
	    puts $g_AnnotSV(benignText)
	    unset g_AnnotSV(benignText)
	}
	puts "\t   >>> ClinVar parsing ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
		
	set benignLossFile_Tmp "$benignDir/benign_Loss_SV_$genomeBuild.tmp.bed"
	set benignGainFile_Tmp "$benignDir/benign_Gain_SV_$genomeBuild.tmp.bed"
	set benignInsFile_Tmp "$benignDir/benign_Ins_SV_$genomeBuild.tmp.bed"
	set benignInvFile_Tmp "$benignDir/benign_Inv_SV_$genomeBuild.tmp.bed"
	set L_toWriteLoss {}
	set L_toWriteGain {}
	set L_toWriteIns {}
	set L_toWriteInv {}
	set f [open "| gzip -cd $ClinVarFileDownloaded"]
	while {! [eof $f]} {
	    set L [gets $f]
	    set Ls [split $L "\t"]
	    if {[regexp "^#" $L]} {continue}
	    set infos [lindex $Ls 7]
	    
	    # Selection of the benign variants to keep:
	    ###########################################
	    # Criteria:
	    # - “benign” or “benign/likely benign” clinical significance (CLNSIG)
	    # - “criteria_provided”, “_multiple_submitters” or “reviewed_by_expert_panel” SV review status (CLNREVSTAT)
	    # - “Deletion” or “Duplication” SV type (CLNVC)
	    # - ≥ $g_AnnotSV(SVminSize) bp in size (default 50)	    
	    if {![regexp "ALLELEID=(\[^;\]+)?;.*CLNREVSTAT=(\[^;\]+)?;.*CLNSIG=Benign.*?CLNVC=(\[^;\]+)?" $infos match ALLELEID CLNREVSTAT CLNVC]} {continue}
	    if {![regexp "criteria_provided|_multiple_submitters|reviewed_by_expert_panel" $CLNREVSTAT]} {continue}
	    if {[regexp "no_assertion_criteria_provided" $CLNREVSTAT]} {continue}
	    
	    set ref [lindex $Ls 3]
	    set alt [lindex $Ls 4]
	    set SVlength [expr {abs([string length $ref]-[string length $alt])}]
	    if {$SVlength < $g_AnnotSV(SVminSize)} {continue}
	    set chrom [lindex $Ls 0]
	    set start [lindex $Ls 1]
	    set end [expr {$start+$SVlength}]
	    set coord "${chrom}:${start}-$end"
	    
	    if {$CLNVC eq "Deletion"} {
		lappend L_toWriteLoss "$chrom\t$start\t$end\tCLN:$ALLELEID\t$coord"
	    } elseif {$CLNVC eq "Duplication"} {
		lappend L_toWriteGain "$chrom\t$start\t$end\tCLN:$ALLELEID\t$coord"
	    } elseif {$CLNVC eq "Insertion"} {
		lappend L_toWriteIns "$chrom\t$start\t$end\tCLN:$ALLELEID\t$coord"
	    } elseif {$CLNVC eq "Inversion"} {
		lappend L_toWriteInv "$chrom\t$start\t$end\tCLN:$ALLELEID\t$coord"
	    }
	    
	}
	close $f
	
	# Writing:
	##########
	puts "\t       ([llength $L_toWriteLoss] SV Loss + [llength $L_toWriteGain] SV Gain)"
 	if {$L_toWriteLoss ne {}} {
	    WriteTextInFile [join $L_toWriteLoss "\n"] $benignLossFile_Tmp
	}
	if {$L_toWriteGain ne {}} {
	    WriteTextInFile [join $L_toWriteGain "\n"] $benignGainFile_Tmp
	}
	if {$L_toWriteIns ne {}} {
	    WriteTextInFile [join $L_toWriteIns "\n"]  $benignInsFile_Tmp
	}
	if {$L_toWriteInv ne {}} {
	    WriteTextInFile [join $L_toWriteInv "\n"]  $benignInvFile_Tmp
	}
	
	# Clean:
	########
	file delete -force $ClinVarFileDownloaded
    }

    return
}


proc checkClinGenHITS_benignFile {genomeBuild} {
    
    global g_AnnotSV
    
    ## Check if ClinGen files have been downloaded 
    ##############################################
    set benignDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/SVincludedInFt/BenignSV/$genomeBuild"
    set ClinGenFileDownloaded1 [glob -nocomplain "$benignDir/ClinGen_gene_curation_list_$genomeBuild.tsv"] 
    set ClinGenFileDownloaded2 [glob -nocomplain "$benignDir/ClinGen_region_curation_list_$genomeBuild.tsv"] 

    
    if {$ClinGenFileDownloaded1 ne "" || $ClinGenFileDownloaded2 ne ""} {
	# We have some ClinGen annotations to add in $benign*File
	if {[info exists g_AnnotSV(benignText)]} {
	    puts $g_AnnotSV(benignText)
	    unset g_AnnotSV(benignText)
	}
	puts "\t   >>> ClinGen parsing ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    } else {return}

    set L_files "$ClinGenFileDownloaded1"
    lappend L_files "$ClinGenFileDownloaded2"

    set L_toWriteLoss {}
    set L_toWriteGain {}
    set benignLossFile_Tmp "$benignDir/benign_Loss_SV_$genomeBuild.tmp.bed"
    set benignGainFile_Tmp "$benignDir/benign_Gain_SV_$genomeBuild.tmp.bed"

    foreach ClinGenFileDownloaded "$L_files" {
	if {$ClinGenFileDownloaded ne ""} {	
	    set f [open "$ClinGenFileDownloaded"]
	    while {! [eof $f]} {
		set L [gets $f]
		set Ls [split $L "\t"]		
		# Header 1 and 2:
		#################
		#Gene Symbol    Gene ID cytoBand        Genomic Location        Haploinsufficiency Score        Haploinsufficiency Description  Haploinsufficiency PMID1        Haploinsufficiency PMID2        Haploinsufficiency PMID3 Triplosensitivity Score  Triplosensitivity Description   Triplosensitivity PMID1 Triplosensitivity PMID2 Triplosensitivity PMID3 Date Last Evaluated     Loss phenotype OMIM ID  Triplosensitive phenotype OMIM ID
		
		#ISCA ID        ISCA Region Name        cytoBand        Genomic Location        Haploinsufficiency Score        Haploinsufficiency Description  Haploinsufficiency PMID1        Haploinsufficiency PMID2        Haploinsufficiency PMID3  Triplosensitivity Score Triplosensitivity Description   Triplosensitivity PMID1 Triplosensitivity PMID2 Triplosensitivity PMID3 Date Last Evaluated     Loss phenotype OMIM ID  Triplosensitive phenotype OMIM ID
		
		# Selection of the benign variants to keep:
		###########################################
		if {[regexp "(^#Gene Symbol)|(^#ISCA ID)" $L]} {
		    set i_id 0 ;# variantaccession
		    set i_coord [lsearch -exact $Ls "Genomic Location"]; if {$i_coord == -1} {puts "$ClinGenFileDownloaded"; puts "Bad header line syntax. Genomic Location column not found - Exit with error"; exit 2}
		    set i_hi [lsearch -exact $Ls "Haploinsufficiency Score"]; if {$i_hi == -1} {puts "$ClinGenFileDownloaded"; puts "Bad header line syntax. Haploinsufficiency Score column not found - Exit with error"; exit 2}
		    set i_ts [lsearch -exact $Ls "Triplosensitivity Score"]; if {$i_ts == -1} {puts "$ClinGenFileDownloaded"; puts "Bad header line syntax. Triplosensitivity column not found - Exit with error"; exit 2}
		    continue
		}
		if {[regexp "^#" $L]} {continue}
		
		set ID [lindex $Ls $i_id]
		set coord [lindex $Ls $i_coord]
		regsub "(chr)| " $coord "" coord
		if { ![regexp "(\[0-9XYMT\]+):(\[0-9\]+)-(\[0-9\]+)" $coord match chrom start end] } {continue} ;# In GRCh38, some coordinates = "tbd" (to be determined)
		set hi [lindex $Ls $i_hi]
		set ts [lindex $Ls $i_ts]
		if {$hi eq "40"} {
		    lappend L_toWriteLoss "$chrom\t$start\t$end\tHI40:$ID\t$coord"
		}
		if {$ts eq "40"} {
		    lappend L_toWriteGain "$chrom\t$start\t$end\tTS40:$ID\t$coord"
		} 
	    }
	    close $f
	}
    }
    
    # Writing:
    ##########
    puts "\t       ([llength $L_toWriteLoss] SV Loss + [llength $L_toWriteGain] SV Gain)"
    if {$L_toWriteLoss ne {}} {
	WriteTextInFile [join $L_toWriteLoss "\n"] $benignLossFile_Tmp
    }
    if {$L_toWriteGain ne {}} {
	WriteTextInFile [join $L_toWriteGain "\n"] $benignGainFile_Tmp
    }
    
    # Clean:
    ########
    file delete -force $ClinGenFileDownloaded1
    file delete -force $ClinGenFileDownloaded2

    return
}




proc checkDGV_benignFile {genomeBuild} {

    global g_AnnotSV

    ## Check if DGV files have been downloaded
    ##########################################
    set benignDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/SVincludedInFt/BenignSV/$genomeBuild"
    set DGVfileDownloaded [glob -nocomplain "$benignDir/GRCh3*_hg*_variants_*.txt"]

    if {$DGVfileDownloaded ne ""} {
	# We have some DGV annotations to add in $benign*File
	if {[info exists g_AnnotSV(benignText)]} {
	    puts $g_AnnotSV(benignText)
	    unset g_AnnotSV(benignText)
	}
	puts "\t   >>> DGV parsing ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	
	set benignLossFile_Tmp "$benignDir/benign_Loss_SV_$genomeBuild.tmp.bed"
	set benignGainFile_Tmp "$benignDir/benign_Gain_SV_$genomeBuild.tmp.bed"
	set benignInsFile_Tmp "$benignDir/benign_Ins_SV_$genomeBuild.tmp.bed"
	set benignInvFile_Tmp "$benignDir/benign_Inv_SV_$genomeBuild.tmp.bed"
	set L_toWriteLoss {}
	set L_toWriteGain {}
	set L_toWriteIns {}
	set L_toWriteInv {}
   
 	set f [open "$DGVfileDownloaded"]
	foreach L [LinesFromFile $DGVfileDownloaded] {
	    set Ls [split $L "\t"]

	    # Header:
	    #########
	    # variantaccession   chr   start   end   varianttype   variantsubtype   reference   pubmedid   method   platform   mergedvariants   supportingvariants   mergedorsample   frequency   samplesize   observedgains   observedlosses   cohortdescription   genes   samples
	    # WARNING: the "frequency" field is empty
	    if {[regexp "^variantaccession" $L]} {
		set i_id 0 ;# variantaccession
		set i_chr [lsearch -exact $Ls "chr"];     if {$i_chr == -1} {puts "$DGVfileDownloaded"; puts "Bad header line syntax. chr column not found - Exit with error"; exit 2}
		set i_start [lsearch -exact $Ls "start"]; if {$i_start == -1} {puts "$DGVfileDownloaded"; puts "Bad header line syntax. start column not found - Exit with error"; exit 2}
		set i_end [lsearch -exact $Ls "end"];     if {$i_end == -1} {puts "$DGVfileDownloaded"; puts "Bad header line syntax. end column not found - Exit with error"; exit 2}
		set i_samplesize [lsearch -exact $Ls "samplesize"];  if {$i_samplesize == -1} {puts "$DGVfileDownloaded"; puts "Bad header line syntax. samplesize column not found - Exit with error"; exit 2}
		set i_obsgains [lsearch -exact $Ls "observedgains"]; if {$i_obsgains == -1} {puts "$DGVfileDownloaded"; puts "Bad header line syntax. observedgains column not found - Exit with error"; exit 2}
		set i_obsloss [lsearch -exact $Ls "observedlosses"]; if {$i_obsloss == -1} {puts "$DGVfileDownloaded"; puts "Bad header line syntax. observedlosses column not found - Exit with error"; exit 2}
		continue
	    }

	    # Selection of the benign variants to keep:
	    ###########################################
	    # ≥ 500 individuals tested 
	    set samplesize [lindex $Ls $i_samplesize]
	    if {$samplesize <= $g_AnnotSV(minTotalNumber)} {continue}
	    # allele frequency > 1%
	    set obsgains [lindex $Ls $i_obsgains]
	    set obsloss [lindex $Ls $i_obsloss]
	    set freqgain [expr {$obsgains*100/(2*$samplesize)}]
	    set freqloss [expr {$obsloss*100/(2*$samplesize)}]
	    set id [lindex $Ls $i_id]
	    set chr [lindex $Ls $i_chr]
	    set start [lindex $Ls $i_start]
	    set end [lindex $Ls $i_end]
	    set coord "${chr}:${start}-$end"
	    set infos "$chr\t$start\t$end\t$id\t$coord"
	    if {$freqgain > 1} {
		lappend L_toWriteGain "$infos"
	    }
	    if {$freqloss > 1} {
		lappend L_toWriteLoss "$infos"
	    }
	}

	# Writing:
	##########
	puts "\t       ([llength $L_toWriteLoss] SV Loss + [llength $L_toWriteGain] SV Gain)"
	if {$L_toWriteLoss ne {}} {
	    WriteTextInFile [join $L_toWriteLoss "\n"] $benignLossFile_Tmp
	}
	if {$L_toWriteGain ne {}} {
	    WriteTextInFile [join $L_toWriteGain "\n"] $benignGainFile_Tmp
	}
	if {$L_toWriteIns ne {}} {
	    WriteTextInFile [join $L_toWriteIns "\n"]  $benignInsFile_Tmp
	}
	if {$L_toWriteInv ne {}} {
	    WriteTextInFile [join $L_toWriteInv "\n"]  $benignInvFile_Tmp
	}
	
	# Clean:
	########
	file delete -force $DGVfileDownloaded

    }
    
    return
}



proc checkGnomAD_benignFile {genomeBuild} {
    
    global g_AnnotSV
    
    ## Check if GnomAD file has been downloaded 
    ###########################################
    set benignDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/SVincludedInFt/BenignSV/$genomeBuild"
    set GnomADfileDownloaded [glob -nocomplain "$benignDir/gnomad_v2.1_sv.sites.bed.gz"] ;# GRCh37 (not yet available in GRCh38)

    if {$GnomADfileDownloaded ne ""} {
	# We have some GnomAD annotations to add in $benign*File
	if {[info exists g_AnnotSV(benignText)]} {
	    puts $g_AnnotSV(benignText)
	    unset g_AnnotSV(benignText)
	}
	puts "\t   >>> GnomAD parsing ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
	
	set benignLossFile_Tmp "$benignDir/benign_Loss_SV_$genomeBuild.tmp.bed"
	set benignGainFile_Tmp "$benignDir/benign_Gain_SV_$genomeBuild.tmp.bed"
	set benignInsFile_Tmp "$benignDir/benign_Ins_SV_$genomeBuild.tmp.bed"
	set benignInvFile_Tmp "$benignDir/benign_Inv_SV_$genomeBuild.tmp.bed"
	set L_toWriteLoss {}
	set L_toWriteGain {}
	set L_toWriteIns {}
	set L_toWriteInv {}

	set f [open "| gzip -cd $GnomADfileDownloaded"]
	while {![eof $f]} {
	    set L [gets $f]
	    if {$L eq ""} {continue}
	    set Ls [split $L "\t"] 

	    # Header:
	    #########
	    #chrom  start   end     name    svtype  ALGORITHMS      BOTHSIDES_SUPPORT       CHR2    CPX_INTERVALS   CPX_TYPE        END2    END     EVIDENCE        HIGH_SR_BACKGROUND      PCRPLUS_DEPLETED        PESR_GT_OVERDISPERSION    POS2    PROTEIN_CODING__COPY_GAIN       PROTEIN_CODING__DUP_LOF PROTEIN_CODING__DUP_PARTIAL     PROTEIN_CODING__INTERGENIC      PROTEIN_CODING__INTRONIC        PROTEIN_CODING__INV_SPAN        PROTEIN_CODING__LOF       PROTEIN_CODING__MSV_EXON_OVR    PROTEIN_CODING__NEAREST_TSS     PROTEIN_CODING__PROMOTER        PROTEIN_CODING__UTR     SOURCE  STRANDS SVLEN   SVTYPE  UNRESOLVED_TYPE UNSTABLE_AF_PCRPLUS     VARIABLE_ACROSS_BATCHES   AN      AC      AF      N_BI_GENOS      N_HOMREF        N_HET   N_HOMALT        FREQ_HOMREF     FREQ_HET        FREQ_HOMALT     MALE_AN MALE_AC MALE_AF MALE_N_BI_GENOS MALE_N_HOMREF   MALE_N_HET      MALE_N_HOMALT     MALE_FREQ_HOMREF        MALE_FREQ_HET   MALE_FREQ_HOMALT        MALE_N_HEMIREF  MALE_N_HEMIALT  MALE_FREQ_HEMIREF       MALE_FREQ_HEMIALT       PAR     FEMALE_AN       FEMALE_AC       FEMALE_AF       FEMALE_N_BI_GENOS FEMALE_N_HOMREF FEMALE_N_HET    FEMALE_N_HOMALT FEMALE_FREQ_HOMREF      FEMALE_FREQ_HET FEMALE_FREQ_HOMALT      POPMAX_AF       AFR_AN  AFR_AC  AFR_AF  AFR_N_BI_GENOS  AFR_N_HOMREF    AFR_N_HET       AFR_N_HOMALT      AFR_FREQ_HOMREF AFR_FREQ_HET    AFR_FREQ_HOMALT AFR_MALE_AN     AFR_MALE_AC     AFR_MALE_AF     AFR_MALE_N_BI_GENOS     AFR_MALE_N_HOMREF       AFR_MALE_N_HET  AFR_MALE_N_HOMALT       AFR_MALE_FREQ_HOMREF    AFR_MALE_FREQ_HET AFR_MALE_FREQ_HOMALT    AFR_MALE_N_HEMIREF      AFR_MALE_N_HEMIALT      AFR_MALE_FREQ_HEMIREF   AFR_MALE_FREQ_HEMIALT   AFR_FEMALE_AN   AFR_FEMALE_AC   AFR_FEMALE_AF   AFR_FEMALE_N_BI_GENOS   AFR_FEMALE_N_HOMREF       AFR_FEMALE_N_HET        AFR_FEMALE_N_HOMALT     AFR_FEMALE_FREQ_HOMREF  AFR_FEMALE_FREQ_HET     AFR_FEMALE_FREQ_HOMALT  AMR_AN  AMR_AC  AMR_AF  AMR_N_BI_GENOS  AMR_N_HOMREF    AMR_N_HET       AMR_N_HOMALT      AMR_FREQ_HOMREF AMR_FREQ_HET    AMR_FREQ_HOMALT AMR_MALE_AN     AMR_MALE_AC     AMR_MALE_AF     AMR_MALE_N_BI_GENOS     AMR_MALE_N_HOMREF       AMR_MALE_N_HET  AMR_MALE_N_HOMALT       AMR_MALE_FREQ_HOMREF    AMR_MALE_FREQ_HET AMR_MALE_FREQ_HOMALT    AMR_MALE_N_HEMIREF      AMR_MALE_N_HEMIALT      AMR_MALE_FREQ_HEMIREF   AMR_MALE_FREQ_HEMIALT   AMR_FEMALE_AN   AMR_FEMALE_AC   AMR_FEMALE_AF   AMR_FEMALE_N_BI_GENOS   AMR_FEMALE_N_HOMREF       AMR_FEMALE_N_HET        AMR_FEMALE_N_HOMALT     AMR_FEMALE_FREQ_HOMREF  AMR_FEMALE_FREQ_HET     AMR_FEMALE_FREQ_HOMALT  EAS_AN  EAS_AC  EAS_AF  EAS_N_BI_GENOS  EAS_N_HOMREF    EAS_N_HET       EAS_N_HOMALT      EAS_FREQ_HOMREF EAS_FREQ_HET    EAS_FREQ_HOMALT EAS_MALE_AN     EAS_MALE_AC     EAS_MALE_AF     EAS_MALE_N_BI_GENOS     EAS_MALE_N_HOMREF       EAS_MALE_N_HET  EAS_MALE_N_HOMALT       EAS_MALE_FREQ_HOMREF    EAS_MALE_FREQ_HET EAS_MALE_FREQ_HOMALT    EAS_MALE_N_HEMIREF      EAS_MALE_N_HEMIALT      EAS_MALE_FREQ_HEMIREF   EAS_MALE_FREQ_HEMIALT   EAS_FEMALE_AN   EAS_FEMALE_AC   EAS_FEMALE_AF   EAS_FEMALE_N_BI_GENOS   EAS_FEMALE_N_HOMREF       EAS_FEMALE_N_HET        EAS_FEMALE_N_HOMALT     EAS_FEMALE_FREQ_HOMREF  EAS_FEMALE_FREQ_HET     EAS_FEMALE_FREQ_HOMALT  EUR_AN  EUR_AC  EUR_AF  EUR_N_BI_GENOS  EUR_N_HOMREF    EUR_N_HET       EUR_N_HOMALT      EUR_FREQ_HOMREF EUR_FREQ_HET    EUR_FREQ_HOMALT EUR_MALE_AN     EUR_MALE_AC     EUR_MALE_AF     EUR_MALE_N_BI_GENOS     EUR_MALE_N_HOMREF       EUR_MALE_N_HET  EUR_MALE_N_HOMALT       EUR_MALE_FREQ_HOMREF    EUR_MALE_FREQ_HET EUR_MALE_FREQ_HOMALT    EUR_MALE_N_HEMIREF      EUR_MALE_N_HEMIALT      EUR_MALE_FREQ_HEMIREF   EUR_MALE_FREQ_HEMIALT   EUR_FEMALE_AN   EUR_FEMALE_AC   EUR_FEMALE_AF   EUR_FEMALE_N_BI_GENOS   EUR_FEMALE_N_HOMREF       EUR_FEMALE_N_HET        EUR_FEMALE_N_HOMALT     EUR_FEMALE_FREQ_HOMREF  EUR_FEMALE_FREQ_HET     EUR_FEMALE_FREQ_HOMALT  OTH_AN  OTH_AC  OTH_AF  OTH_N_BI_GENOS  OTH_N_HOMREF    OTH_N_HET       OTH_N_HOMALT      OTH_FREQ_HOMREF OTH_FREQ_HET    OTH_FREQ_HOMALT OTH_MALE_AN     OTH_MALE_AC     OTH_MALE_AF     OTH_MALE_N_BI_GENOS     OTH_MALE_N_HOMREF       OTH_MALE_N_HET  OTH_MALE_N_HOMALT       OTH_MALE_FREQ_HOMREF    OTH_MALE_FREQ_HET OTH_MALE_FREQ_HOMALT    OTH_MALE_N_HEMIREF      OTH_MALE_N_HEMIALT      OTH_MALE_FREQ_HEMIREF   OTH_MALE_FREQ_HEMIALT   OTH_FEMALE_AN   OTH_FEMALE_AC   OTH_FEMALE_AF   OTH_FEMALE_N_BI_GENOS   OTH_FEMALE_N_HOMREF       OTH_FEMALE_N_HET        OTH_FEMALE_N_HOMALT     OTH_FEMALE_FREQ_HOMREF  OTH_FEMALE_FREQ_HET     OTH_FEMALE_FREQ_HOMALT  FILTER
	    if {[regexp "^#chrom" $L]} {
		set i_chrom      [lsearch -exact $Ls "#chrom"];       if {$i_chrom == -1} {puts "Bad header line syntax. #CHROM column not found - Exit with error"; exit 2}
		set i_start      [lsearch -exact $Ls "start"];        if {$i_start == -1} {puts "Bad header line syntax. START column not found - Exit with error"; exit 2}
		set i_end        [lsearch -exact $Ls "end"];          if {$i_end == -1} {puts "Bad header line syntax. END column not found - Exit with error"; exit 2}
		set i_svid       [lsearch -exact $Ls "name"];         if {$i_svid == -1} {puts "Bad header line syntax. NAME column not found - Exit with error"; exit 2}
		set i_svtype     [lsearch -exact $Ls "svtype"];       if {$i_svtype == -1} {puts "Bad header line syntax. SVTYPE column not found - Exit with error"; exit 2}
		set i_nhomalt    [lsearch -exact $Ls "N_HOMALT"];     if {$i_nhomalt == -1} {puts "Bad header line syntax. N_HOMALT column not found - Exit with error"; exit 2}
		set i_popmaxaf   [lsearch -exact $Ls "POPMAX_AF"];    if {$i_popmaxaf == -1} {puts "Bad header line syntax. POPMAX_AF column not found - Exit with error"; exit 2}
		continue
	    }

	    # Selection of the benign variants to keep:
	    ###########################################
	    set SVTYPE [lindex $Ls $i_svtype]
	    regsub ":.+" $SVTYPE "" SVTYPE	    
	    # WARNING:
	    # - MCNV have multiple coma-separated-values for nhet, nhomalt, af and popmaxaf... 
	    # - "gnomad_v2_sv.sites.bed"      => 445858 lines of which only 1148 are MCNV
	    #   "gnomad_v2.1_sv.sites.bed.gz" => 387478 lines of which only 1108 are MCNV (<=> "CN=0")
	    # => These SV are not selected.
	    if {[lsearch -exact {DEL DUP INS INV} "$SVTYPE"] eq -1} {continue}
	    
	    # At least one population allele frequency > 1% OR at least 5 homozygous individuals
	    set nhomalt    [lindex $Ls $i_nhomalt]
	    set popmaxaf   [lindex $Ls $i_popmaxaf]
	    if {$popmaxaf < 0.01 && $nhomalt < 5} {continue}

	    set chrom      [lindex $Ls $i_chrom]
	    set start      [lindex $Ls $i_start]
	    set end        [lindex $Ls $i_end]
	    set svid       [lindex $Ls $i_svid]
	    set coord "$chrom:${start}-$end"
	    set infos "$chrom\t$start\t$end\t$svid\t$coord"
	    
	    if {$SVTYPE eq "DEL"} {
		lappend L_toWriteLoss "$infos"
	    }
	    if {$SVTYPE eq "DUP"} {
		lappend L_toWriteGain "$infos"
	    }
	    if {$SVTYPE eq "INS"} {
		lappend L_toWriteIns "$infos"
	    }
	    if {$SVTYPE eq "INV"} {
		lappend L_toWriteInv "$infos"
	    }
	}
	close $f
	
	# Writing:
	##########
	puts "\t       ([llength $L_toWriteLoss] SV Loss + [llength $L_toWriteGain] SV Gain)"
 	if {$L_toWriteLoss ne {}} {
	    WriteTextInFile [join $L_toWriteLoss "\n"] $benignLossFile_Tmp
	}
	if {$L_toWriteGain ne {}} {
	    WriteTextInFile [join $L_toWriteGain "\n"] $benignGainFile_Tmp
	}
	if {$L_toWriteIns ne {}} {
	    WriteTextInFile [join $L_toWriteIns "\n"]  $benignInsFile_Tmp
	}
	if {$L_toWriteInv ne {}} {
	    WriteTextInFile [join $L_toWriteInv "\n"]  $benignInvFile_Tmp
	}  
	# Clean
	file delete -force $GnomADfileDownloaded
    }

    return
}


proc checkDDD_benignFile {genomeBuild} {

    global g_AnnotSV

    ## Check if DDD file has been downloaded 
    ########################################
    set benignDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/SVincludedInFt/BenignSV/$genomeBuild"
    set DDDfileDownloaded [glob -nocomplain "$benignDir/population_cnv.txt.gz"]

    if {$DDDfileDownloaded ne ""} {
	# We have some DDD annotations to add in $benign*File
	if {[info exists g_AnnotSV(benignText)]} {
	    puts $g_AnnotSV(benignText)
	    unset g_AnnotSV(benignText)
	}
	puts "\t   >>> DDD parsing ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"

	set benignLossFile_Tmp "$benignDir/benign_Loss_SV_$genomeBuild.tmp.bed"
	set benignGainFile_Tmp "$benignDir/benign_Gain_SV_$genomeBuild.tmp.bed"
	set L_toWriteLoss {}
	set L_toWriteGain {}
	
	set f [open "| gzip -cd $DDDfileDownloaded"]
	while {! [eof $f]} {
	    set L [gets $f]
	    set Ls [split $L "\t"]

	    # Header:
	    #########
	    # #population_cnv_id      chr     start   end     deletion_observations   deletion_frequency      deletion_standard_error duplication_observations        duplication_frequency   duplication_standard_error      observations     frequency       standard_error  type    sample_size     study

	    if {[regexp  "^#population_cnv_id" $L]} {		
		set i_id 0 ;# #population_cnv_id
		set i_chr      [lsearch -regexp $Ls "chr"];     if {$i_chr == -1} {puts "Bad syntax into $DDDfileDownloaded.\nchr field not found - Exit with error"; exit 2}
		set i_start    [lsearch -regexp $Ls "start"];   if {$i_start == -1} {puts "Bad syntax into $DDDfileDownloaded.\nstart field not found - Exit with error"; exit 2}
		set i_end      [lsearch -regexp $Ls "end"];     if {$i_end == -1} {puts "Bad syntax into $DDDfileDownloaded.\nend field not found - Exit with error"; exit 2}
		set i_delfreq  [lsearch -regexp $Ls "deletion_frequency"];    if {$i_delfreq == -1} {puts "Bad syntax into $DDDfileDownloaded.\ndeletion_frequency field not found - Exit with error"; exit 2}
		set i_dupfreq  [lsearch -regexp $Ls "duplication_frequency"]; if {$i_dupfreq == -1} {puts "Bad syntax into $DDDfileDownloaded.\nduplication_frequency field not found - Exit with error"; exit 2}
		set i_sample_size [lsearch -regexp $Ls "sample_size"];        if {$i_sample_size == -1} {puts "Bad syntax into $DDDfileDownloaded.\nsample_size field not found - Exit with error"; exit 2}
		continue
	    }

	    # Selection of the benign variants to keep:
	    ###########################################
	    set Ls [split $L "\t"]

	    set sample_size [lindex $Ls $i_sample_size]
	    # ≥ 500 individuals tested 
	    if {$sample_size < $g_AnnotSV(minTotalNumber)} {continue}
	    # allele frequency > 1%
	    set id "DDD:[lindex $Ls $i_id]"
	    set chr [lindex $Ls $i_chr]
	    set start [lindex $Ls $i_start]
	    set end [lindex $Ls $i_end]
	    set delfreq [lindex $Ls $i_delfreq]
	    set dupfreq [lindex $Ls $i_dupfreq]
	    set coord "$chr:${start}-$end"
	    set infos "$chr\t$start\t$end\t$id\t$coord"
	    if {$delfreq > 0.01} {
		lappend L_toWriteLoss "$infos"
	    }
	    if {$dupfreq > 0.01} {
		lappend L_toWriteGain "$infos"
	    }
	}
	close $f
	
	# Writing:
	##########
	puts "\t       ([llength $L_toWriteLoss] SV Loss + [llength $L_toWriteGain] SV Gain)"
 	if {$L_toWriteLoss ne {}} {
	    WriteTextInFile [join $L_toWriteLoss "\n"] $benignLossFile_Tmp
	}
	if {$L_toWriteGain ne {}} {
	    WriteTextInFile [join $L_toWriteGain "\n"] $benignGainFile_Tmp
	}
	
	# Clean:
	########
	file delete -force $DDDfileDownloaded
    }

    return
}

proc check1000g_benignFile {genomeBuild} {

    global g_AnnotSV

    ## Check if 1000g file has been downloaded 
    ##########################################
    set benignDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/SVincludedInFt/BenignSV/$genomeBuild"
    set 1000gFileDownloaded [glob -nocomplain "$benignDir/ALL.wgs.mergedSV*.vcf.gz"]

    if {$1000gFileDownloaded ne ""} {
	# We have some 1000g annotations to add in $benign*File
	if {[info exists g_AnnotSV(benignText)]} {
	    puts $g_AnnotSV(benignText)
	    unset g_AnnotSV(benignText)
	}
	puts "\t   >>> 1000g parsing ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"

	set benignLossFile_Tmp "$benignDir/benign_Loss_SV_$genomeBuild.tmp.bed"
	set benignGainFile_Tmp "$benignDir/benign_Gain_SV_$genomeBuild.tmp.bed"
	set benignInsFile_Tmp "$benignDir/benign_Ins_SV_$genomeBuild.tmp.bed"
	set benignInvFile_Tmp "$benignDir/benign_Inv_SV_$genomeBuild.tmp.bed"
	set L_toWriteLoss {}
	set L_toWriteGain {}
	set L_toWriteIns {}
	set L_toWriteInv {}

	# Split multiallelic sites into multiple rows:
	##############################################
	regsub ".vcf.gz" $1000gFileDownloaded ".tmp.vcf" 1000gFileTmp 	
	catch {eval exec $g_AnnotSV(bcftools) norm -m -both $1000gFileDownloaded > $1000gFileTmp} Message
	if {[file size $1000gFileTmp] eq 0} {
	    # we continue AnnotSV without splitting the multiallelic sites !
	    puts "\t   -- check1000g --"
	    puts "\t   $g_AnnotSV(bcftools) norm -m -both $1000gFileDownloaded > $1000gFileTmp"
	    puts "\t   $1000gFileTmp: file empty."
	    puts "\t   No multiallelic treatment done."
	    file delete -force $1000gFileTmp	   
	} 
	
	# Selection of the benign variants to keep:
	###########################################
	if {[file exists $1000gFileTmp]} {
	    set f [open "$1000gFileTmp"]
	} else {
	    set f [open "| gzip -cd $1000gFileDownloaded"]
	}
	while {![eof $f]} {
	    set L [gets $f]
	    if {[string index $L 0] eq "#" || $L eq ""} {continue}
	    # Line example:
	    # 1  645710  ALU_umary_ALU_2 A  <INS:ME:ALU>  .  .  AC=35;AF=0.00698882;AFR_AF=0;AMR_AF=0.0072;AN=5008;CS=ALU_umary;EAS_AF=0.0069;EUR_AF=0.0189;MEINFO=AluYa4_5,1,223,-;NS=2504;SAS_AF=0.0041;SITEPOST=0.9998;SVLEN=222;SVTYPE=ALU;TSD=null  GT  0|0  ...
	    set Ls [split $L "\t"]

	    set chrom [lindex $Ls 0]
	    set pos [lindex $Ls 1]
 	    set L_infos [split [lindex $Ls 7] ";"]
	    set ref [lindex $Ls 3]
	    set alt [lindex $Ls 4]
	    # Consider only the SV (not the SNV/indel in the VCF file)
	    ##########################################################
	    # Example of SV:
	    # - Type1: ref="G" and alt="ACTGCTAACGATCCGTTTGCTGCTAACGATCTAACGATCGGGATTGCTAACGATCTCGGG"
	    # - Type2: "<INS>", "<DEL>", ...
	    # - Type3: complex rearrangements with breakends: "G]17:1584563]"
	    if {![regexp "<|\\\[|\\\]" $alt]} {
		set variantLength [expr {[string length $ref]-[string length $alt]}]
		if {[expr {abs($variantLength)}]<$g_AnnotSV(SVminSize)} {continue}; # it is not an SV
	    }
	    set max -1
	    catch {unset AF}
	    catch {unset END}
	    set SVTYPE "CNV"
	    foreach inf $L_infos {
		set inf [split $inf "="]
		set val [lindex $inf 0]
		set $val [lindex $inf 1]
		if {[regexp "AF$" $val]} {
		    if {[set $val]>$max} {set max [set $val]}
		}
	    }
	    # allele frequency > 1%
	    if {$max < 0.01} {continue}
	    
	    if {![info exists END]} {
		# INS:ME (LINE1, ALU or SVA)
		set END [expr {$pos+1}]
	    }
	    # SVTYPE is defined in $L_infos
	    # In case of a single value for "alt", SVTYPE is correctly defined (e.g. SVTYPE=DUP).
	    # In case of multiple alt ("<CN0>,<CN2>"), SVTYPE has a unique value (often: SVTYPE=CNV)
	    #   => We replace CNV with a more precise value <CN0> or <CNV2>
	    if {$SVTYPE eq "CNV"} {set SVTYPE "$alt"}
	    # Sometimes, we have "SVTYPE=DUP" and alt=<CN0>
	    # or "SVTYPE=DEL" and alt=<CN2>
	    # It's due to multiallelic sites. The good value is in "alt"
	    if {[regexp "^<CN" $alt]} {
		set SVTYPE "$alt"
	    }
	    set coord "$chrom:${pos}-$END"
	    set infos "$chrom\t$pos\t$END\t1000g\t$coord"
	    
	    if {[regexp "<CN0>|DEL" $SVTYPE]} {
		lappend L_toWriteLoss "$infos"
	    }
	    if {[regexp "<CN\[1-9\]|ALU|LINE1|SVA" $SVTYPE]} {
		lappend L_toWriteGain "$infos"
	    }
	    if {$SVTYPE eq "INS"} {
		lappend L_toWriteIns "$infos"
	    }
	    if {$SVTYPE eq "INV"} {
		lappend L_toWriteInv "$infos"
	    }
	}
    
	close $f
	
	# Writing:
	##########
	puts "\t       ([llength $L_toWriteLoss] SV Loss + [llength $L_toWriteGain] SV Gain)"
 	if {$L_toWriteLoss ne {}} {
	    WriteTextInFile [join $L_toWriteLoss "\n"] $benignLossFile_Tmp
	}
	if {$L_toWriteGain ne {}} {
	    WriteTextInFile [join $L_toWriteGain "\n"] $benignGainFile_Tmp
	}
	if {$L_toWriteIns ne {}} {
	    WriteTextInFile [join $L_toWriteIns "\n"]  $benignInsFile_Tmp
	}
	if {$L_toWriteInv ne {}} {
	    WriteTextInFile [join $L_toWriteInv "\n"]  $benignInvFile_Tmp
	}
	
	# Clean:
	########
	file delete -force $1000gFileDownloaded
	file delete -force $1000gFileTmp	   
    }

    return
}

proc checkIMH_benignFile {genomeBuild} {

    global g_AnnotSV

    ## Check if IMH file has been downloaded 
    ########################################
    set benignDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/SVincludedInFt/BenignSV/$genomeBuild"
    set IMHfileDownloaded [glob -nocomplain "$benignDir/*callset.public.bedpe.gz"]

    if {$IMHfileDownloaded ne ""} {
	# We have some IMH annotations to add in $benign*File
	if {[info exists g_AnnotSV(benignText)]} {
	    puts $g_AnnotSV(benignText)
	    unset g_AnnotSV(benignText)
	}
	puts "\t   >>> IMH parsing ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"

	set benignLossFile_Tmp "$benignDir/benign_Loss_SV_$genomeBuild.tmp.bed"
	set benignGainFile_Tmp "$benignDir/benign_Gain_SV_$genomeBuild.tmp.bed"
	set benignInsFile_Tmp "$benignDir/benign_Ins_SV_$genomeBuild.tmp.bed"
	set benignInvFile_Tmp "$benignDir/benign_Inv_SV_$genomeBuild.tmp.bed"
	set L_toWriteLoss {}
	set L_toWriteGain {}
	set L_toWriteIns {}
	set L_toWriteInv {}

	set f [open "| gzip -cd $IMHfileDownloaded"]
	while {![eof $f]} {
	    set L [gets $f]
	    # Header(bedpe):
	    ################
	    #CHROM_A        START_A END_A   CHROM_B START_B END_B   ID      QUAL    STRAND_A        STRAND_B        TYPE    FILTER  NAME_A  REF_A   ALT_A   NAME_B  REF_B   ALT_B   INFO_A  INFO_B
	    #1       10361   10590   15      102520301       102520447       2838    47911.3 -       -       BND     LOW     2838_1  N       [15:102520406[N 2838_2  N       [1:10567[N      SVTYPE=BND;POS=10567;STRANDS=--:97;IMPRECISE;CIPOS=-205,24;CIEND=-104,42;CIPOS95=-20,20;CIEND95=-6,6;MATEID=2838_2;EVENT=2838;SU=97;PE=97;SR=0;ALG=PROD;AF=0.07568;NSAMP=1274;MSQ=37.02;AN=16834;AC=1274;NS=8438;AC_Hom=0;AC_Het=1275;AC_Hemi=0;MAF=0.0755511;HWE=4.03332e-23;NFAM=1263     SVTYPE=BND;POS=102520406;STRANDS=--:97;IMPRECISE;CIPOS=-104,42;CIEND=-205,24;CIPOS95=-6,6;CIEND95=-20,20;MATEID=2838_1;EVENT=2838;SECONDARY;SU=97;PE=97;SR=0;ALG=PROD;AF=0.07568;NSAMP=1274;MSQ=37.02;AN=16834;AC=1274;NS=8438;AC_Hom=0;AC_Het=1275;AC_Hemi=0;MAF=0.0755511;HWE=4.03332e-23;NFAM=1263
	    set Ls [split $L "\t"] 
	    if {[regexp  "^#CHROM_A" $L]} {		
		set i_chrA 0 ;# #CHROM_A
		set i_chrB [lsearch -regexp $Ls "CHROM_B"]; if {$i_chrB == -1} {puts "Bad syntax into $IMHfileDownloaded.\nCHROM_B field not found - Exit with error"; exit 2}
		set i_start [lsearch -regexp $Ls "START_A"]; if {$i_start == -1} {puts "Bad syntax into $IMHfileDownloaded.\nSTART_A field not found - Exit with error"; exit 2}
		set i_end [lsearch -regexp $Ls "END_B"]; if {$i_end == -1} {puts "Bad syntax into $IMHfileDownloaded.\nEND_B field not found - Exit with error"; exit 2}
		set i_svtype [lsearch -regexp $Ls "TYPE"]; if {$i_svtype == -1} {puts "Bad syntax into $IMHfileDownloaded.\nTYPE field not found - Exit with error"; exit 2}
		set i_annotations [lsearch -regexp $Ls "INFO_A"]; if {$i_annotations == -1} {puts "Bad syntax into $IMHfileDownloaded.\nINFO_A field not found - Exit with error"; exit 2}
	    }
	    if {[string index $L 0] eq "#" || $L eq ""} {continue}

	    # Selection of the benign variants to keep:
	    ###########################################	    
	    regsub "chr" [lindex $Ls $i_chrA] "" chromA
	    regsub "chr" [lindex $Ls $i_chrB] "" chromB 
	    if {$chromA ne $chromB} {continue}	    
	    set Annotations [lindex $Ls $i_annotations]
	    if {![regexp ";AF=(.+?);" $Annotations match AF]} {continue}
	    if {$AF < 0.01} {continue}
	    set start  [lindex $Ls $i_start]
	    set end    [lindex $Ls $i_end]
	    set SVTYPE [normalizeSVtype [lindex $Ls $i_svtype]]
	    set coord "$chromA:${start}-$end"
	    set infos "$chromA\t$start\t$end\tIMH\t$coord"
	    
	    if {$SVTYPE eq "DEL"} {
		lappend L_toWriteLoss "$infos"
	    }
	    if {$SVTYPE eq "DUP"} {
		lappend L_toWriteGain "$infos"
	    }
	    if {$SVTYPE eq "INS"} {
		lappend L_toWriteIns "$infos"
	    }
	    if {$SVTYPE eq "INV"} {
		lappend L_toWriteInv "$infos"
	    }
	}
	close $f

	# Writing:
	##########
	puts "\t       ([llength $L_toWriteLoss] SV Loss + [llength $L_toWriteGain] SV Gain)"

 	if {$L_toWriteLoss ne {}} {
	    WriteTextInFile [join $L_toWriteLoss "\n"] $benignLossFile_Tmp
	}
	if {$L_toWriteGain ne {}} {
	    WriteTextInFile [join $L_toWriteGain "\n"] $benignGainFile_Tmp
	}
	if {$L_toWriteIns ne {}} {
	    WriteTextInFile [join $L_toWriteIns "\n"]  $benignInsFile_Tmp
	}
	if {$L_toWriteInv ne {}} {
	    WriteTextInFile [join $L_toWriteInv "\n"]  $benignInvFile_Tmp
	}
 	# Clean:
	########
	file delete -force $IMHfileDownloaded
	file delete -force [glob -nocomplain "Build*.public.v2.vcf.gz"]
	file delete -force [glob -nocomplain "Supplementary_File_*.zip"]
    }
    
    return
}

# 2020-11-27
############
# benign_Gain_SV_GRCh37.sorted.bed:  23 148 SV (  310 161 514 bp)
# benign_Loss_SV_GRCh37.sorted.bed:  60 442 SV (  433 579 309 bp)
# benign_Ins_SV_GRCh37.sorted.bed:   21 254 SV (    1 638 472 bp)
# benign_Inv_SV_GRCh37.sorted.bed:    1 693 SV (5 381 728 519 bp)




# Return the overlapped SV annotation 
proc benignSVannotation {SVchrom SVstart SVend} {

    global g_AnnotSV
    global benignText

    set benignDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/SVincludedInFt/BenignSV/$g_AnnotSV(genomeBuild)"
    
    if {![info exists benignText(DONE)]} {
	# headerOutput "B_gain_source B_gain_coord B_loss_source B_loss_coord"
	# + (if user selected) "B_ins_source B_ins_coord B_inv_source B_inv_coord" 

	set L_benignText(Empty) {}
	foreach svtype {"gain" "loss" "ins" "inv"} {
	    # Keep only the user requested columns (defined in the configfile)
	    if {[lsearch -regexp "$g_AnnotSV(outputColHeader)" "^B_${svtype}_"] eq -1} { continue }
	    lappend L_benignText(Empty) {*}{"" ""}
	}
	set benignText(Empty) "[join $L_benignText(Empty) "\t"]"
	
	foreach svtype {"Gain" "Loss" "Ins" "Inv"} {
	    
	    # Keep only the user requested columns (defined in the configfile)
	    if {[lsearch -regexp "$g_AnnotSV(outputColHeader)" "^B_[string tolower ${svtype}]_"] eq -1} { continue }
	    
	    # Intersect
	    set benignBEDfile [glob -nocomplain "$benignDir/benign_${svtype}_SV_$g_AnnotSV(genomeBuild).sorted.bed"]
	    regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".intersect.benign-$svtype" tmpFile
	    set tmpFile "$g_AnnotSV(outputDir)/[file tail $tmpFile]"
	    file delete -force $tmpFile
	    # -f 1.0 : the feature in B overlaps at least 100% of the A feature.
	    if {[catch {exec $g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(fullAndSplitBedFile) -b $benignBEDfile -f 1.0 -wa -wb > $tmpFile} Message]} {
		puts "-- benignSVAnnotation, $svtype --"
		puts "$g_AnnotSV(bedtools) intersect -sorted -a $g_AnnotSV(fullAndSplitBedFile) -b $benignBEDfile -f 1.0 -wa -wb > $tmpFile"
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

		set benign_coord [lindex $Ls end]
		set benign_source [lindex $Ls end-1]
		
		set SVtoAnn "$SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end"
		lappend L_allSVtoAnn $SVtoAnn
		lappend L_benign_coord($SVtoAnn,$svtype) $benign_coord
		lappend L_benign_source($SVtoAnn,$svtype) $benign_source
	    }
	    file delete -force $tmpFile
	}
	    
	# Loading benign final annotation for each SV
	if {[info exists L_allSVtoAnn]} {
	    foreach SVtoAnn [lsort -unique $L_allSVtoAnn] {
		
		foreach svtype {"Gain" "Loss" "Ins" "Inv"} {
		    
		    # Keep only the user requested columns (defined in the configfile)
		    if {[lsearch -regexp "$g_AnnotSV(outputColHeader)" "^B_[string tolower ${svtype}]_"] eq -1} { continue }
			
		    if {[info exists L_benign_coord($SVtoAnn,$svtype)]} {
			lappend L_benignText($SVtoAnn) "[join $L_benign_source($SVtoAnn,$svtype) " / "]"
			lappend L_benignText($SVtoAnn) "[join $L_benign_coord($SVtoAnn,$svtype) " / "]"
		    } else {
			lappend L_benignText($SVtoAnn) ""
			lappend L_benignText($SVtoAnn) ""
		    }
		}
		set benignText($SVtoAnn) [join $L_benignText($SVtoAnn) "\t"]
	    }
	}
	set benignText(DONE) 1	
    }
    
    if {[info exist benignText($SVchrom,$SVstart,$SVend)]} {
	return $benignText($SVchrom,$SVstart,$SVend)
    } else {
	return $benignText(Empty)
    }
}
