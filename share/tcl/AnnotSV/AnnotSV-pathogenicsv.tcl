############################################################################################################
# AnnotSV 3.0.8                                                                                            #
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


# Creation/update (if some data source files are presents) of:
# - $pathogenicLossFile_Sorted
# - $pathogenicGainFile_Sorted
# - $pathogenicInsFile_Sorted
# - $pathogenicInvFile_Sorted
proc checkPathogenicFiles {} {

    global g_AnnotSV

    foreach genomeBuild {GRCh37 GRCh38} {
	set pathogenicDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/FtIncludedInSV/PathogenicSV/$genomeBuild"

	# Files to create/update
	set pathogenicLossFile_Sorted "$pathogenicDir/pathogenic_Loss_SV_$genomeBuild.sorted.bed"
	set pathogenicGainFile_Sorted "$pathogenicDir/pathogenic_Gain_SV_$genomeBuild.sorted.bed"
	set pathogenicInsFile_Sorted "$pathogenicDir/pathogenic_Ins_SV_$genomeBuild.sorted.bed"
	set pathogenicInvFile_Sorted "$pathogenicDir/pathogenic_Inv_SV_$genomeBuild.sorted.bed"
	
	# Text to write only if an update of the sources is requested (in the check* proc).
	set g_AnnotSV(pathogenicText) "\t...creation/update of the \"pathogenic_*_SV_$genomeBuild.sorted.bed\" files ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])\n\t   (done only once)"
	
	# Creation of pathogenic*File_Tmp (not formatted and sorted)
	set pathogenicLossFile_Tmp "$pathogenicDir/pathogenic_Loss_SV_$genomeBuild.tmp.bed"
	set pathogenicGainFile_Tmp "$pathogenicDir/pathogenic_Gain_SV_$genomeBuild.tmp.bed"
	set pathogenicInsFile_Tmp "$pathogenicDir/pathogenic_Ins_SV_$genomeBuild.tmp.bed"
	set pathogenicInvFile_Tmp "$pathogenicDir/pathogenic_Inv_SV_$genomeBuild.tmp.bed"
	file delete -force $pathogenicLossFile_Tmp
	file delete -force $pathogenicGainFile_Tmp
	file delete -force $pathogenicInsFile_Tmp
	file delete -force $pathogenicInvFile_Tmp	
	checkClinVar_PathogenicFile $genomeBuild
	checkClinGenHITS_pathogenicFile $genomeBuild
	checkdbVar_PathogenicFile $genomeBuild
	checkOMIM_PathogenicFile $genomeBuild
	catch {unset g_AnnotSV(pathogenicText)}
    
	# Creation of *.formatted.bed
	foreach SVtype {Loss Gain Ins Inv} {
	    set pathogenicFile_Sorted "$pathogenicDir/pathogenic_${SVtype}_SV_$genomeBuild.sorted.bed"
	    set pathogenicFile_Tmp "$pathogenicDir/pathogenic_${SVtype}_SV_$genomeBuild.tmp.bed"
	    if {[file exists $pathogenicFile_Tmp]} {
		# Add the SV from $pathogenic*File_Sorted 
		if {[file exists $pathogenicFile_Sorted]} {
		    set id_file1 [open "$pathogenicFile_Sorted" r]
		    set id_file2 [open "$pathogenicFile_Tmp" a]
		    while { [gets $id_file1 line] >= 0 } {
			puts $id_file2 $line
		    }
		    close $id_file1
		    close $id_file2
		}
		# sort
		if {[catch {checkBed $pathogenicFile_Tmp $pathogenicDir} Message]} {
		    puts "-- checkPathogenicFiles --"
		    puts $Message
		}

		file delete -force $pathogenicFile_Tmp
	    }
	}
	
	# Creation of *.sorted.bed:
	# Intersection with very large files can cause trouble with excessive memory usage.
	# A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm. 
	foreach SVtype {Loss Gain Ins Inv} {
	    set pathogenicFile_TmpFormatted "$pathogenicDir/pathogenic_${SVtype}_SV_$genomeBuild.tmp.formatted.bed"
	    if {![file exist $pathogenicFile_TmpFormatted]} {continue}
	    set pathogenicFile_Sorted "$pathogenicDir/pathogenic_${SVtype}_SV_$genomeBuild.sorted.bed"
	    set sortTmpFile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_sort.tmp.bash"
	    ReplaceTextInFile "#!/bin/bash" $sortTmpFile
	    WriteTextInFile "# The locale specified by the environment can affects the traditional sort order. We need to use native byte values." $sortTmpFile
	    WriteTextInFile "export LC_ALL=C" $sortTmpFile
	    WriteTextInFile "sort -k1,1 -k2,2n $pathogenicFile_TmpFormatted > $pathogenicFile_Sorted" $sortTmpFile
	    file attributes $sortTmpFile -permissions 0755
	    if {[catch {eval exec bash $sortTmpFile} Message]} {
		puts "-- checkPathogenicFiles --"
		puts "sort -k1,1 -k2,2n $pathogenicFile_TmpFormatted > $pathogenicFile_Sorted"
		puts "$Message"
		puts "Exit with error"
		exit 2
	    }
	    file delete -force $sortTmpFile 
	    file delete -force $pathogenicFile_TmpFormatted
	}
    }
    
    # Remove some columns from g_AnnotSV(outputColHeader) if the corresponding annotation file doesn't exist
    foreach svtype {"Gain" "Loss" "Ins" "Inv"} {
	set pathogenicFile_Sorted "$pathogenicDir/pathogenic_${svtype}_SV_$genomeBuild.sorted.bed"
	if {![file exists $pathogenicFile_Sorted]} {
	    set newList {}
	    foreach e "$g_AnnotSV(outputColHeader)" {
		if {[regexp "^P_[string tolower ${svtype}]_" $e]} {continue}
		lappend newList $e
	    }
	    set g_AnnotSV(outputColHeader) $newList
	}
    }

    return
}




################################################################################################
## In all the next "check*File" proc, create and then complete the following file:
##   - Pathogenic_Gain_SV_$genomeBuild.tmp.bed
##   - Pathogenic_Loss_SV_$genomeBuild.tmp.bed
##   - Pathogenic_Ins_SV_$genomeBuild.tmp.bed
##   - Pathogenic_Inv_SV_$genomeBuild.tmp.bed
################################################################################################

proc checkClinVar_PathogenicFile {genomeBuild} {
    
    global g_AnnotSV
    
    ## Check if ClinVar file has been downloaded 
    ############################################
    set pathogenicDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/FtIncludedInSV/PathogenicSV/$genomeBuild"
    set ClinVarFileDownloaded [lindex [glob -nocomplain "$pathogenicDir/clinvar*.vcf.gz"] end]

    if {$ClinVarFileDownloaded ne ""} {
	# We have some ClinVar annotations to add in $pathogenic*File
	if {[info exists g_AnnotSV(pathogenicText)]} {
	    puts $g_AnnotSV(pathogenicText)
	    unset g_AnnotSV(pathogenicText)
	}
	puts "\t   >>> ClinVar parsing ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
		
	set pathogenicLossFile_Tmp "$pathogenicDir/pathogenic_Loss_SV_$genomeBuild.tmp.bed"
	set pathogenicGainFile_Tmp "$pathogenicDir/pathogenic_Gain_SV_$genomeBuild.tmp.bed"
	set pathogenicInsFile_Tmp "$pathogenicDir/pathogenic_Ins_SV_$genomeBuild.tmp.bed"
	set pathogenicInvFile_Tmp "$pathogenicDir/pathogenic_Inv_SV_$genomeBuild.tmp.bed"
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
	    
	    # Selection of the pathogenic variants to keep:
	    ###########################################
	    # Criteria:
	    # - “pathogenic” or “pathogenic/likely pathogenic” clinical significance (CLNSIG)
	    # - “criteria_provided”, “_multiple_submitters” or “reviewed_by_expert_panel” SV review status (CLNREVSTAT)
	    # - "Deletion" or "Duplication" or "Insertion" or "Inversion" SV type (CLNVC)
	    # - ≥ $g_AnnotSV(SVminSize) bp in size (default 50)	    
	    if {![regexp "ALLELEID=(\[^;\]+)?;.*CLNDISDB=(\[^;\]+)?;.*CLNDN\[INCL\]*=(\[^;\]+)?;.*CLNREVSTAT=(\[^;\]+)?;.*CLNSIG=Pathogenic.*?CLNVC=(\[^;\]+)?" $infos match ALLELEID CLNDISDB CLNDN CLNREVSTAT CLNVC]} {continue}
	    if {![regexp "criteria_provided|_multiple_submitters|reviewed_by_expert_panel" $CLNREVSTAT]} {continue}
	    if {[regexp "no_assertion_criteria_provided" $CLNREVSTAT]} {continue}
	    
	    set ref [lindex $Ls 3]
	    set alt [lindex $Ls 4]
	    set SVlength [expr {abs([string length $ref]-[string length $alt])+1}]
	    if {$SVlength < $g_AnnotSV(SVminSize)} {continue}
	    set chrom [lindex $Ls 0]
	    set start [lindex $Ls 1]
	    set end [expr {$start+$SVlength}]
	    set coord "${chrom}:${start}-$end"

	    # CLNDISDB example :
	    # Human_Phenotype_Ontology:HP:0000556,Human_Phenotype_Ontology:HP:0007736,Human_Phenotype_Ontology:HP:0007910,Human_Phenotype_Ontology:HP:0007974,Human_Phenotype_Ontology:HP:0007982,MONDO:MONDO:0019118,MeSH:D058499,MedGen:C0854723,Orphanet:ORPHA71862,SNOMED_CT:314407005|MONDO:MONDO:0012056,MedGen:C1837873,OMIM:608553|MedGen:CN517202
	    set L_HPO {}
	    set i 0
	    while {[regexp -start $i -indices "Human_Phenotype_Ontology:HP:\[0-9\]+" $CLNDISDB indices]} {
		set HPO [string range $CLNDISDB [lindex $indices 0] [lindex $indices 1]]
		regsub "^Human_Phenotype_Ontology:" $HPO "" HPO
		lappend L_HPO $HPO
		set i [lindex $indices 1]
	    }
	    regsub "(\\\|)?not_provided" $CLNDN "" CLNDN
	    set toWrite "$chrom\t$start\t$end\t$CLNDN\t$L_HPO\tCLN:$ALLELEID\t$coord"
	    if {$CLNVC eq "Deletion"} {
		lappend L_toWriteLoss "$toWrite"
	    } elseif {$CLNVC eq "Duplication"} {
		lappend L_toWriteGain "$toWrite"
	    } elseif {$CLNVC eq "Insertion"} {
		lappend L_toWriteIns "$toWrite"
	    } elseif {$CLNVC eq "Inversion"} {
		lappend L_toWriteInv "$toWrite"
	    }
	    
	}
	close $f
	
	# Writing:
	##########
	puts "\t       ([llength $L_toWriteLoss] SV Loss + [llength $L_toWriteGain] SV Gain)"
 	if {$L_toWriteLoss ne {}} {
	    WriteTextInFile [join $L_toWriteLoss "\n"] $pathogenicLossFile_Tmp
	}
	if {$L_toWriteGain ne {}} {
	    WriteTextInFile [join $L_toWriteGain "\n"] $pathogenicGainFile_Tmp
	}
	if {$L_toWriteIns ne {}} {
	    WriteTextInFile [join $L_toWriteIns "\n"]  $pathogenicInsFile_Tmp
	}
	if {$L_toWriteInv ne {}} {
	    WriteTextInFile [join $L_toWriteInv "\n"]  $pathogenicInvFile_Tmp
	}
	
	# Clean:
	########
	file delete -force $ClinVarFileDownloaded
    }

    return
}


proc checkClinGenHITS_pathogenicFile {genomeBuild} {
    
    global g_AnnotSV
    
    ## Check if ClinGen files have been downloaded 
    ##############################################
    set pathogenicDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/FtIncludedInSV/PathogenicSV/$genomeBuild"
    set ClinGenFileDownloaded1 [glob -nocomplain "$pathogenicDir/ClinGen_gene_curation_list_$genomeBuild.tsv"] 
    set ClinGenFileDownloaded2 [glob -nocomplain "$pathogenicDir/ClinGen_region_curation_list_$genomeBuild.tsv"] 
    
    if {$ClinGenFileDownloaded1 ne "" || $ClinGenFileDownloaded2 ne ""} {
	# We have some ClinGen annotations to add in $pathogenic*File
	if {[info exists g_AnnotSV(pathogenicText)]} {
	    puts $g_AnnotSV(pathogenicText)
	    unset g_AnnotSV(pathogenicText)
	}
	puts "\t   >>> ClinGen parsing ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    } else {return}

    set L_files "$ClinGenFileDownloaded1"
    lappend L_files "$ClinGenFileDownloaded2"

    set L_toWriteLoss {}
    set L_toWriteGain {}
    set pathogenicLossFile_Tmp "$pathogenicDir/pathogenic_Loss_SV_$genomeBuild.tmp.bed"
    set pathogenicGainFile_Tmp "$pathogenicDir/pathogenic_Gain_SV_$genomeBuild.tmp.bed"

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
		
		# Selection of the pathogenic variants to keep:
		###########################################
		if {[regexp "(^#Gene Symbol)|(^#ISCA ID)" $L]} {
		    set i_id 0 ;# variantaccession
		    set i_coord [lsearch -exact $Ls "Genomic Location"]; if {$i_coord == -1} {puts "$ClinGenFileDownloaded"; puts "Bad header line syntax. Genomic Location column not found - Exit with error"; exit 2}
		    set i_hi [lsearch -exact $Ls "Haploinsufficiency Score"]; if {$i_hi == -1} {puts "$ClinGenFileDownloaded"; puts "Bad header line syntax. Haploinsufficiency Score column not found - Exit with error"; exit 2}
		    set i_ts [lsearch -exact $Ls "Triplosensitivity Score"]; if {$i_ts == -1} {puts "$ClinGenFileDownloaded"; puts "Bad header line syntax. Triplosensitivity Score column not found - Exit with error"; exit 2}
		    set i_omimloss [lsearch -exact $Ls "Loss phenotype OMIM ID"]; if {$i_omimloss == -1} {puts "$ClinGenFileDownloaded"; puts "Bad header line syntax. Loss phenotype OMIM ID column not found - Exit with error"; exit 2}
		    set i_omimgain [lsearch -exact $Ls "Triplosensitive phenotype OMIM ID"]; if {$i_omimgain == -1} {puts "$ClinGenFileDownloaded"; puts "Bad header line syntax. Triplosensitive phenotype OMIM ID column not found - Exit with error"; exit 2}
		    continue
		}
		if {[regexp "^#" $L]} {continue}
		
		set ID [lindex $Ls $i_id]
		set coord [lindex $Ls $i_coord]
		regsub -all "(chr)| " $coord "" coord
		if { ![regexp "(\[0-9XYMT\]+):(\[0-9\]+)-(\[0-9\]+)" $coord match chrom start end] } {continue} ;# In GRCh38, some coordinates = "tbd" (to be determined)
		set hi [lindex $Ls $i_hi]
		set ts [lindex $Ls $i_ts]
		set hpo ""
		if {$hi eq "3"} {
		    set OMIM [lindex $Ls $i_omimloss]
		    set phenotype [fromOMIMtoPhenotype $OMIM]
		    lappend L_toWriteLoss "$chrom\t$start\t$end\t$phenotype\t$hpo\tHI3:$ID\t$coord"
		}
		if {$ts eq "3"} {
		    set OMIM [lindex $Ls $i_omimgain]
		    set phenotype [fromOMIMtoPhenotype $OMIM]
		    lappend L_toWriteGain "$chrom\t$start\t$end\t$phenotype\t$hpo\tTS3:$ID\t$coord"
		} 
	    }
	    close $f
	}
    }
    
    # Writing:
    ##########
    puts "\t       ([llength $L_toWriteLoss] SV Loss + [llength $L_toWriteGain] SV Gain)"
    if {$L_toWriteLoss ne {}} {
	WriteTextInFile [join $L_toWriteLoss "\n"] $pathogenicLossFile_Tmp
    }
    if {$L_toWriteGain ne {}} {
	WriteTextInFile [join $L_toWriteGain "\n"] $pathogenicGainFile_Tmp
    }
    
    # Clean:
    ########
    file delete -force $ClinGenFileDownloaded1
    file delete -force $ClinGenFileDownloaded2

    return
}


proc checkdbVar_PathogenicFile {genomeBuild} {

    global g_AnnotSV

    ## Check if dbVar files have been downloaded
    ############################€###############
    set pathogenicDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/FtIncludedInSV/PathogenicSV/$genomeBuild"
    set dbVarFileDownloaded1 [glob -nocomplain "$pathogenicDir/$genomeBuild.nr_deletions.pathogenic.tsv.gz"]    ; # Loss
    set dbVarFileDownloaded2 [glob -nocomplain "$pathogenicDir/$genomeBuild.nr_duplications.pathogenic.tsv.gz"] ; # Gain
    set dbVarFileDownloaded3 [glob -nocomplain "$pathogenicDir/$genomeBuild.nr_insertions.pathogenic.tsv.gz"]   ; # Ins

    if {$dbVarFileDownloaded1 ne "" || $dbVarFileDownloaded2 ne "" || $dbVarFileDownloaded3 ne ""} {
	# We have some dbVar annotations to add in $pathogenic*File
	if {[info exists g_AnnotSV(pathogenicText)]} {
	    puts $g_AnnotSV(pathogenicText)
	    unset g_AnnotSV(pathogenicText)
	}
	puts "\t   >>> dbVar parsing ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    } else {return}

    set L_files "$dbVarFileDownloaded1"
    lappend L_files "$dbVarFileDownloaded2"
    lappend L_files "$dbVarFileDownloaded3"

    set L_toWriteLoss {}
    set L_toWriteGain {}
    set L_toWriteIns {}
    set pathogenicLossFile_Tmp "$pathogenicDir/pathogenic_Loss_SV_$genomeBuild.tmp.bed"
    set pathogenicGainFile_Tmp "$pathogenicDir/pathogenic_Gain_SV_$genomeBuild.tmp.bed"
    set pathogenicInsFile_Tmp "$pathogenicDir/pathogenic_Ins_SV_$genomeBuild.tmp.bed"

    foreach dbVarFileDownloaded "$L_files" {

	if {[regexp "deletions" $dbVarFileDownloaded]} {
	    set svtype "Loss"
	} elseif {[regexp "duplications" $dbVarFileDownloaded]} {
	    set svtype "Gain"
	} elseif {[regexp "insertions" $dbVarFileDownloaded]} {
	    set svtype "Ins"
	}

	if {$dbVarFileDownloaded ne ""} {
	    set f [open "| gzip -cd $dbVarFileDownloaded"]
	    while {! [eof $f]} {
		set L [gets $f]
		set Ls [split $L "\t"]		
		# Header 1 and 2:
		#################
		# #chr    outermost_start outermost_stop  variant_count   variant_type    method  analysis        platform        study   variant clinical_assertion      clinvar_accession       bin_size

		# Selection of the pathogenic variants to keep:
		###############################################
		if {[regexp "^#chr" $L]} {
		    set i_chr 0
		    set i_start [lsearch -exact $Ls "outermost_start"]; if {$i_start == -1} {puts "$dbVarFileDownloaded"; puts "Bad header line syntax. outermost_start column not found - Exit with error"; exit 2}
		    set i_stop [lsearch -exact $Ls "outermost_stop"]; if {$i_stop == -1} {puts "$dbVarFileDownloaded"; puts "Bad header line syntax. outermost_stop column not found - Exit with error"; exit 2}
		    set i_source [lsearch -exact $Ls "variant"]; if {$i_source == -1} {puts "$dbVarFileDownloaded"; puts "Bad header line syntax. variant column not found - Exit with error"; exit 2}
		    
		    continue
		}
		if {[string index $L 0] eq "#" || $L eq ""} {continue}
		set chrom [lindex $Ls $i_chr]
		set start [lindex $Ls $i_start]
		set end [lindex $Ls $i_stop]
		set Source [lindex $Ls $i_source]
		set coord "$chrom:${start}-$end"
		set phenotype ""
		set hpo ""

		lappend L_toWrite$svtype "$chrom\t$start\t$end\t$phenotype\t$hpo\tdbVar:$Source\t$coord"
	    }
	    close $f
	}
    }
    
    # Writing:
    ##########
    puts "\t       ([llength $L_toWriteLoss] SV Loss + [llength $L_toWriteGain] SV Gain)"
    if {$L_toWriteLoss ne {}} {
	WriteTextInFile [join $L_toWriteLoss "\n"] $pathogenicLossFile_Tmp
    }
    if {$L_toWriteGain ne {}} {
	WriteTextInFile [join $L_toWriteGain "\n"] $pathogenicGainFile_Tmp
    }
    if {$L_toWriteIns ne {}} {
	WriteTextInFile [join $L_toWriteIns "\n"] $pathogenicInsFile_Tmp
    }
    
    # Clean:
    ########
    file delete -force $dbVarFileDownloaded1
    file delete -force $dbVarFileDownloaded2
    file delete -force $dbVarFileDownloaded3  

    return
}


proc checkOMIM_PathogenicFile {genomeBuild} {

    global g_AnnotSV

    ## Check if OMIM file has been downloaded 
    #########################################
    set pathogenicDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/FtIncludedInSV/PathogenicSV/$genomeBuild"
    set OMIMfileDownloaded [glob -nocomplain "$pathogenicDir/*_morbid.tsv.gz"]


    if {$OMIMfileDownloaded ne ""} {
	# We have some OMIM annotations to add in $pathogenicLossFile
	if {[info exists g_AnnotSV(pathogenicText)]} {
	    puts $g_AnnotSV(pathogenicText)
	    unset g_AnnotSV(pathogenicText)
	}
	puts "\t   >>> OMIM morbid genes parsing ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"

	# Genes coordinates
	set genesDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Genes/$g_AnnotSV(genomeBuild)"
	set genesFileFormatted "[glob -nocomplain $genesDir/genes.$g_AnnotSV(tx).sorted.bed]"
	foreach L [LinesFromFile $genesFileFormatted] {
	    set Ls [split $L "\t"]
	    set geneName [lindex $Ls 4]
	    set chrom($geneName) [lindex $Ls 0]
	    set start($geneName) [lindex $Ls 1]
	    set end($geneName)   [lindex $Ls 2]
	}

	# Morbid genes 
	set L_toWriteLoss {}
	set f [open "| gzip -cd $OMIMfileDownloaded"]
	set L_HPO ""
	while {! [eof $f]} {
	    set L [gets $f]
	    set g [lindex $L 0]
	    if {$g eq ""} {continue}
	    set pheno [fromOMIMgeneToPhenotype $g]
	    if {[info exists chrom($g)]} {
		set coord "$chrom($g):$start($g)"
		append coord "-$end($g)"
		regsub -all "{|}" $pheno "" pheno
		set toWrite "$chrom($g)\t$start($g)\t$end($g)\t$pheno\t$L_HPO\tmorbid:$g\t$coord"
		lappend L_toWriteLoss "$toWrite"
	    }
	}
	close $f
	
	# Writing:
	##########
	puts "\t       ([llength $L_toWriteLoss] SV Loss)"
	set pathogenicLossFile_Tmp "$pathogenicDir/pathogenic_Loss_SV_$genomeBuild.tmp.bed"
 	if {$L_toWriteLoss ne {}} {
	    WriteTextInFile [join $L_toWriteLoss "\n"] $pathogenicLossFile_Tmp
	}
	
	# Clean:
	########
	file delete -force $OMIMfileDownloaded
    }

   
    return
}




# Return the overlapped SV annotation 
proc pathogenicSVannotation {SVchrom SVstart SVend} {

    global g_AnnotSV
    global pathogenicText

    set pathogenicDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/FtIncludedInSV/PathogenicSV/$g_AnnotSV(genomeBuild)"
    
    if {![info exists pathogenicText(DONE)]} {
	# headerOutput:         P_gain_phen  P_gain_hpo  P_gain_source  P_gain_coord 
	#                       P_loss_phen  P_loss_hpo  P_loss_source  P_loss_coord
	# + (if user selected): P_ins_phen   P_ins_hpo   P_ins_source   P_ins_coord
	#                       P_inv_phen   P_inv_hpo   P_inv_source   P_inv_coord

	set L_pathogenicText(Empty) {}
	foreach svtype {"gain" "loss" "ins" "inv"} {
	    # Keep only the user requested columns (defined in the configfile)
	    if {[lsearch -regexp "$g_AnnotSV(outputColHeader)" "^P_${svtype}_"] eq -1} { continue }
	    lappend L_pathogenicText(Empty) {*}{"" "" "" ""}
	}
	set pathogenicText(Empty) "[join $L_pathogenicText(Empty) "\t"]"
	
	foreach svtype {"Gain" "Loss" "Ins" "Inv"} {
	    
	    # Keep only the user requested columns (defined in the configfile)
	    if {[lsearch -regexp "$g_AnnotSV(outputColHeader)" "^P_[string tolower ${svtype}]_"] eq -1} { continue }
	    
	    # Intersect
	    set pathogenicBEDfile [glob -nocomplain "$pathogenicDir/pathogenic_${svtype}_SV_$g_AnnotSV(genomeBuild).sorted.bed"]
	    regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".intersect.pathogenic-$svtype" tmpFile
	    set tmpFile "$g_AnnotSV(outputDir)/[file tail $tmpFile]"
	    file delete -force $tmpFile
	    # -f 1.0 : the feature in B overlaps at least 100% of the A feature.
	    if {[catch {exec $g_AnnotSV(bedtools) intersect -sorted -a $pathogenicBEDfile -b $g_AnnotSV(fullAndSplitBedFile) -f 1.0 -wa -wb > $tmpFile} Message]} {
		if {[catch {exec $g_AnnotSV(bedtools) intersect -a $pathogenicBEDfile -b $g_AnnotSV(fullAndSplitBedFile) -f 1.0 -wa -wb > $tmpFile} Message]} {
		    puts "-- pathogenicSVAnnotation, $svtype --"
		    puts "$g_AnnotSV(bedtools) intersect -sorted -a $pathogenicBEDfile -b $g_AnnotSV(fullAndSplitBedFile) -f 1.0 -wa -wb > $tmpFile"
		    puts "$Message"
		    puts "Exit with error"
		    exit 2
		}
	    }
	    # Parse
	    set f [open $tmpFile]
	    while {![eof $f]} {
		set L [gets $f]
		if {$L eq ""} {continue}
		set Ls [split $L "\t"]

		set pathogenic_phen [lindex $Ls 3]
		set pathogenic_hpo [lindex $Ls 4]
		set pathogenic_source [lindex $Ls 5]
		set pathogenic_coord [lindex $Ls 6]

		set SVtoAnn_chrom [lindex $Ls 7]
		set SVtoAnn_start [lindex $Ls 8]
		set SVtoAnn_end   [lindex $Ls 9]
		
		set SVtoAnn "$SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end"
		lappend L_allSVtoAnn $SVtoAnn
		if {$pathogenic_phen ne ""} {
		    lappend L_pathogenic_phen($SVtoAnn,$svtype) $pathogenic_phen
		}
		if {$pathogenic_hpo ne ""} {
		    lappend L_pathogenic_hpo($SVtoAnn,$svtype) {*}[split $pathogenic_hpo " "]
		}
		lappend L_pathogenic_source($SVtoAnn,$svtype) $pathogenic_source
		lappend L_pathogenic_coord($SVtoAnn,$svtype) $pathogenic_coord
	    }
	    file delete -force $tmpFile
	}
	    
	# Loading pathogenic final annotation for each SV
	if {[info exists L_allSVtoAnn]} {
	    foreach SVtoAnn [lsort -unique $L_allSVtoAnn] {
		
		foreach svtype {"Gain" "Loss" "Ins" "Inv"} {
		    
		    # Keep only the user requested columns (defined in the configfile)
		    if {[lsearch -regexp "$g_AnnotSV(outputColHeader)" "^P_[string tolower ${svtype}]_"] eq -1} { continue }
			
		    if {[info exists L_pathogenic_coord($SVtoAnn,$svtype)]} {
			if {[info exists L_pathogenic_phen($SVtoAnn,$svtype)]} {
			    lappend L_pathogenicText($SVtoAnn) "[join [lsort -unique $L_pathogenic_phen($SVtoAnn,$svtype)] ";"]"
			} else {
			    lappend L_pathogenicText($SVtoAnn) ""
			}
			if {[info exists L_pathogenic_hpo($SVtoAnn,$svtype)]} {
			    lappend L_pathogenicText($SVtoAnn) "[join [lsort -unique $L_pathogenic_hpo($SVtoAnn,$svtype)] ";"]"
			} else {
			    lappend L_pathogenicText($SVtoAnn) ""
			}
			lappend L_pathogenicText($SVtoAnn) "[join [lsort -unique $L_pathogenic_source($SVtoAnn,$svtype)] ";"]"
			lappend L_pathogenicText($SVtoAnn) "[join [lsort -unique $L_pathogenic_coord($SVtoAnn,$svtype)] ";"]"
		    } else {
			lappend L_pathogenicText($SVtoAnn) {*}{"" "" "" ""}
		    }
		}
		set pathogenicText($SVtoAnn) [join $L_pathogenicText($SVtoAnn) "\t"]
	    }
	}
	set pathogenicText(DONE) 1	
    }
    
    if {[info exist pathogenicText($SVchrom,$SVstart,$SVend)]} {
	return $pathogenicText($SVchrom,$SVstart,$SVend)
    } else {
	return $pathogenicText(Empty)
    }
}
