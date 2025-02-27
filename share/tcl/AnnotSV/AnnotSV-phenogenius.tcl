############################################################################################################
# AnnotSV 3.4.4                                                                                            #
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



## - Check if the PhenoGeniusCli installation is ok
##   Install of PhenoGeniusCli is done here (if needed)
proc checkPhenoGeniusCli {} {
    
    global g_AnnotSV

	if {$g_AnnotSV(hpo) eq ""} {
		set g_AnnotSV(PhenoGeniusCli) "0" ;# No HPO terms in INPUT
		return
	}
	set checkResult [catch {exec $g_AnnotSV(bashDir)/checkPhenoGeniusCliInstall.sh} Message]
    if {$checkResult eq 0} {
        puts "\tINFO: AnnotSV takes use of PhenoGenius (Yauy et al., 2023) for the phenotype-driven analysis."
        set g_AnnotSV(PhenoGeniusCli) "1"
    } else {
        puts "\nWARNING: No PhenoGenius installation available for the phenotype-driven analysis."
		puts "\ncf $g_AnnotSV(installDir)/share/python3/phenogeniuscli/*.log\n"
        set g_AnnotSV(PhenoGeniusCli) "0"
    }
 
    return
}



# Creation of the g_PhenoGeniusCli variable:
# g_PhenoGeniusCli($geneName) = PhenoGeniusCli_score\tPhenoGeniusCli_phenotype\tPhenoGeniusCli_specificity
# g_PhenoGeniusCli($NCBIid) = PhenoGeniusCli_score\tPhenoGeniusCli_phenotype\tPhenoGeniusCli_specificity
# default = "-1.0\t\t"
# 
# INPUTS:
# L_NCBI_ID
# L_genes
# L_HPO:   e.g. "HP:0001156,HP:0001363,HP:0011304,HP:0010055"
#
# INFO: This proc is run from "AnnotSV-regulatoryelements.tcl"
proc runPhenoGeniusCli {L_Genes L_NCBI_ID L_HPO} {
    
    global g_AnnotSV
    global g_PhenoGeniusCli
    
    puts "...running PhenoGeniusCli for [llength $L_Genes] gene names ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])\n"
    
    set tmpCommandFile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_run_phenogeniuscli.tmp.bash"
    set tmpResultsFile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_phenogeniuscli_results.tmp.tsv"
    set tmpResultsLogFile "$g_AnnotSV(outputDir)/[clock format [clock seconds] -format "%Y%m%d-%H%M%S"]_phenogeniuscli_results.tmp.log"

	set codeToWrite "#!/bin/bash\n\n"
	append codeToWrite "cd $g_AnnotSV(installDir)/share/python3/phenogeniuscli/PhenoGeniusCli\n"
    append codeToWrite "poetry install\n"
    append codeToWrite "poetry run python3 phenogenius_cli.py --hpo_list $L_HPO --result_file $tmpResultsFile &> $tmpResultsLogFile \n"

    WriteTextInFile "$codeToWrite" $tmpCommandFile
	file attributes $tmpCommandFile -permissions 0755

	if {[catch {exec $tmpCommandFile} Message] && [regexp -nocase "error" $Message]} {
		puts "$codeToWrite\n"
		puts $Message
	}

	if {![file exists $tmpResultsLogFile]} {set g_AnnotSV(PhenoGeniusCli) "0"; return}

	foreach L [LinesFromFile $tmpResultsFile] {
		set Ls [split $L "\t"]
		# Header
		# gene_id gene_symbol     rank    score   hpo_implicated  hpo_description_implicated      phenotype_specificity
		if {[regexp "^#gene_id" $L]} {
			set i_ncbi_id "0"
			set i_gene [lsearch -exact $Ls "gene_symbol"]
			set i_score [lsearch -exact $Ls "score"]
			set i_phenotype [lsearch -exact $Ls "hpo_description_implicated"]
			set i_specificity [lsearch -exact $Ls "phenotype_specificity"]
			continue
		}
		# Line example:
		# 5310    PKD1    1       3.2    [{'HP:0000108': 0.7}, {'HP:0005562': 0.5}]    [{'Renal corticomedullary cysts': 0.7}, {'Multiple renal cysts': 0.5}]    A - the reported phenotype is highly specific and relatively unique to the gene (top 40, 50 perc of diagnosis in PhenoGenius cohort).
		set NCBI_ID [lindex $Ls $i_ncbi_id]
		set gene [lindex $Ls $i_gene]
        set score [lindex $Ls $i_score]
        regsub -all "': \[0-9\]\\\.\[0-9\]" [lindex $Ls $i_phenotype] "" phenotype
		regsub -all "\[\\\[\\\]{}'\]" $phenotype "" phenotype
        set specificity [string index [lindex $Ls $i_specificity] 0] ;# A, B, C or D
		if {[lsearch -exact $L_Genes $gene] ne -1} {
			set g_PhenoGeniusCli($gene) "$score\t$phenotype\t$specificity"
		}
		if {[lsearch -exact $L_NCBI_ID $NCBI_ID] ne -1} {
            set g_PhenoGeniusCli($NCBI_ID) "$score\t$phenotype\t$specificity"
        }
	}

	file delete -force $tmpCommandFile
    file delete -force $tmpResultsFile
    file delete -force $tmpResultsLogFile

    return ""
}


# g_PhenoGeniusCli($geneName) = PhenoGeniusCli_score\tPhenoGeniusCli_phenotype\tPhenoGeniusCli_specificity
# g_PhenoGeniusCli($NCBIid) = PhenoGeniusCli_score\tPhenoGeniusCli_phenotype\tPhenoGeniusCli_specificity
# Return either only the specificity or all the annotation:
# what = "specificity" or "all"
proc PhenoGeniusCliAnnotation {GeneName what} {
    
    global g_PhenoGeniusCli
    
    if {$what eq "all"} {
        # Return all the annotation (PhenoGeniusCli_score\tPhenoGeniusCli_phenotype\tPhenoGeniusCli_specificity) => for gene annotations
        if {[info exists g_PhenoGeniusCli($GeneName)]} {
            return $g_PhenoGeniusCli($GeneName)
		} else {
			set NCBI_ID [searchforGeneID $GeneName]
			if {[info exists g_PhenoGeniusCli($NCBI_ID)]} {
				return $g_PhenoGeniusCli($NCBI_ID)
			} else {
				return "-1.0\t\t"
			}
		}
    } elseif {$what eq "specificity"} {
        # Return only the specificity (PhenoGeniusCli_specificity) => for regulatory elements annotations
        if {[info exists g_PhenoGeniusCli($GeneName)]} {
            return [lindex [split $g_PhenoGeniusCli($GeneName) "\t"] end]
        } else {
            set NCBI_ID [searchforGeneID $GeneName]
            if {[info exists g_PhenoGeniusCli($NCBI_ID)]} {
                return [lindex [split $g_PhenoGeniusCli($NCBI_ID) "\t"] end]
            } else {
                return ""
            }
        }
    } else {
        puts "proc PhenoGeniusCliAnnotation: Bad option value for \"what\" ($what)"
    }
}

