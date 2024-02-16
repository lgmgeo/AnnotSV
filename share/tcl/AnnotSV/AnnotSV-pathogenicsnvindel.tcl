############################################################################################################
# AnnotSV 3.4                                                                                              #
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


# Creation/update (if some data source files are presents) of:
# - $pathoSNVindelBEDfile
proc checkPathoSNVindelFile {} {
    
    global g_AnnotSV
    
    foreach genomeBuild {GRCh37 GRCh38} {
        ## Check if ClinVar file has been downloaded
        ############################################
        set pathoSNVindelDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/FtIncludedInSV/PathogenicSNVindel/$genomeBuild"
        set ClinVarFileDownloaded [lindex [glob -nocomplain "$pathoSNVindelDir/clinvar*.vcf.gz"] end]
        
        if {$ClinVarFileDownloaded ne ""} {
            # We have some ClinVar annotations to parse
            puts "\t   ...$genomeBuild pathogenic SNV/indels configuration ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
            
            set pathoSNVindelBEDfile "$pathoSNVindelDir/pathogenic_SNVindel_$genomeBuild.bed"
            file delete -force $pathoSNVindelBEDfile
            set L_toWrite {}
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
                # - < $g_AnnotSV(SVminSize) bp in size (default 50)
                if {![regexp "ALLELEID=(\[^;\]+)?;.*CLNDISDB=(\[^;\]+)?;.*CLNDN\[INCL\]*=(\[^;\]+)?;.*CLNREVSTAT=(\[^;\]+)?;.*CLNSIG=Pathogenic" $infos match ALLELEID CLNDISDB CLNDN CLNREVSTAT]} {continue}
                if {![regexp "criteria_provided|_multiple_submitters|reviewed_by_expert_panel" $CLNREVSTAT]} {continue}
                if {[regexp "no_assertion_criteria_provided" $CLNREVSTAT]} {continue}
                
                set ref [lindex $Ls 3]
                set alt [lindex $Ls 4]
                set SVlength [expr {abs([string length $ref]-[string length $alt])+1}]
                if {$SVlength >= $g_AnnotSV(SVminSize)} {continue}
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
                lappend L_toWrite "$toWrite"
            }
            close $f
            
            # Writing:
            ##########
            puts "\t       ([llength $L_toWrite] pathogenic SNV/indel)"
            if {$L_toWrite ne {}} {
                WriteTextInFile [join $L_toWrite "\n"] $pathoSNVindelBEDfile
            }
            
            # Clean:
            ########
            file delete -force $ClinVarFileDownloaded
            
            # Creation of *.formatted.bed
            #############################
            if {[file exists $pathoSNVindelBEDfile]} {
                regsub ".bed$" $pathoSNVindelBEDfile ".formatted.bed" pathoSNVindelBEDfile_Formatted
                if {[catch {checkBed $pathoSNVindelBEDfile $pathoSNVindelDir} Message]} {
                    puts "-- checkPathoSNVindelFile --"
                    puts $Message
                }
                file delete -force $pathoSNVindelBEDfile
            }
            
            # Creation of *.sorted.bed
            ##########################
            if {[file exists $pathoSNVindelBEDfile_Formatted]} {
                regsub ".bed$" $pathoSNVindelBEDfile ".sorted.bed" pathoSNVindelBEDfile_Sorted
                # Intersection with very large files can cause trouble with excessive memory usage.
                # A presort of the bed files by chromosome and then by start position combined with the use of the -sorted option will invoke a memory-efficient algorithm.
                set sortTmpFile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_sort.tmp.bash"
                ReplaceTextInFile "#!/bin/bash" $sortTmpFile
                WriteTextInFile "# The locale specified by the environment can affects the traditional sort order. We need to use native byte values." $sortTmpFile
                WriteTextInFile "export LC_ALL=C" $sortTmpFile
                WriteTextInFile "sort -k1,1 -k2,2n $pathoSNVindelBEDfile_Formatted > $pathoSNVindelBEDfile_Sorted" $sortTmpFile
                file attributes $sortTmpFile -permissions 0755
                if {[catch {eval exec bash $sortTmpFile} Message]} {
                    puts "-- checkPathoSNVindelFile --"
                    puts "sort -k1,1 -k2,2n $pathoSNVindelBEDfile_Formatted > $pathoSNVindelBEDfile_Sorted"
                    puts "$Message"
                    puts "Exit with error"
                    exit 2
                }
                file delete -force $sortTmpFile
                file delete -force $pathoSNVindelBEDfile_Formatted
            }
        }
    }
    
    return
}



# Return the overlapped SV annotation
proc pathoSNVindelAnnotation {SVchrom SVstart SVend} {
    
    global g_AnnotSV
    global pathoSNVindelText
    
    set pathoSNVindelDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/FtIncludedInSV/PathogenicSNVindel/$g_AnnotSV(genomeBuild)"
    set pathoSNVindelBEDfile "$pathoSNVindelDir/pathogenic_SNVindel_$g_AnnotSV(genomeBuild).sorted.bed"
    
    if {![info exists pathoSNVindelText(DONE)]} {
        # headerOutput:  P_snvindel_nb P_snvindel_phen
        set pathoSNVindelText(Empty) "\t"
        
        # Intersect
        if {[file exists $pathoSNVindelBEDfile]} {
            regsub -nocase "(.formatted)?.bed$" $g_AnnotSV(bedFile) ".intersect.pathogenic-snvindel" tmpFile
            set tmpFile "$g_AnnotSV(outputDir)/[file tail $tmpFile]"
            file delete -force $tmpFile
            # -f 1.0 : the feature in B overlaps at least 100% of the A feature.
            if {[catch {exec $g_AnnotSV(bedtools) intersect -sorted -a $pathoSNVindelBEDfile -b $g_AnnotSV(fullAndSplitBedFile) -f 1.0 -wa -wb > $tmpFile} Message]} {
                if {[catch {exec $g_AnnotSV(bedtools) intersect -a $pathoSNVindelBEDfile -b $g_AnnotSV(fullAndSplitBedFile) -f 1.0 -wa -wb > $tmpFile} Message]} {
                    puts "-- pathoSNVindelAnnotation --"
                    puts "$g_AnnotSV(bedtools) intersect -sorted -a $pathoSNVindelBEDfile -b $g_AnnotSV(fullAndSplitBedFile) -f 1.0 -wa -wb > $tmpFile"
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
                #set pathogenic_hpo [lindex $Ls 4]
                #set pathogenic_source [lindex $Ls 5]
                #set pathogenic_coord [lindex $Ls 6]
                
                set SVtoAnn_chrom [lindex $Ls 7]
                set SVtoAnn_start [lindex $Ls 8]
                set SVtoAnn_end   [lindex $Ls 9]
                
                set SVtoAnn "$SVtoAnn_chrom,$SVtoAnn_start,$SVtoAnn_end"
                lappend L_allSVtoAnn $SVtoAnn
                foreach p [split $pathogenic_phen "|"] {
                    if {$p ne ""} {lappend L_pathogenic_phen($SVtoAnn) "$p"}
                }
                if {![info exists L_pathogenic_phen($SVtoAnn)]} {set L_pathogenic_phen($SVtoAnn) {}}
                #lappend L_pathogenic_hpo($SVtoAnn) {*}[split $pathogenic_hpo "|"]
                #lappend L_pathogenic_source($SVtoAnn) $pathogenic_source
                #lappend L_pathogenic_coord($SVtoAnn) $pathogenic_coord
                if {[info exists nb_snvindel($SVtoAnn)]} {
                    incr nb_snvindel($SVtoAnn)
                } else {
                    set nb_snvindel($SVtoAnn) 1
                }
            }
            file delete -force $tmpFile
            
            # Loading pathogenic snvindel final annotation for each SV
            if {[info exists L_allSVtoAnn]} {
                foreach SVtoAnn [lsort -unique $L_allSVtoAnn] {
                    lappend L_pathoSNVindelText($SVtoAnn) $nb_snvindel($SVtoAnn)
                    lappend L_pathoSNVindelText($SVtoAnn) "[join [lsort -unique $L_pathogenic_phen($SVtoAnn)] ";"]"
                    #lappend L_pathoSNVindelText($SVtoAnn) "[join $L_pathogenic_hpo($SVtoAnn) ";"]"
                    #lappend L_pathoSNVindelText($SVtoAnn) "[join $L_pathogenic_source($SVtoAnn) ";"]"
                    #lappend L_pathoSNVindelText($SVtoAnn) "[join $L_pathogenic_coord($SVtoAnn) ";"]"
                    set pathoSNVindelText($SVtoAnn) [join $L_pathoSNVindelText($SVtoAnn) "\t"]
                }
            }
        }
        set pathoSNVindelText(DONE) 1
    }
    
    if {[info exist pathoSNVindelText($SVchrom,$SVstart,$SVend)]} {
        return $pathoSNVindelText($SVchrom,$SVstart,$SVend)
    } else {
        return $pathoSNVindelText(Empty)
    }
}

