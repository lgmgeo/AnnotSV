############################################################################################################
# AnnotSV 3.5.5                                                                                            #
#                                                                                                          #
# AnnotSV: An integrated tool for Structural Variations annotation and ranking                             #
#                                                                                                          #
# Copyright (C) 2017-present Veronique Geoffroy (veronique.geoffroy@inserm.fr)                             #
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

proc checkNCBIandHGNC {} {
    
    global g_AnnotSV
    
    # Check the existence of the "$NCBIandHGNCgeneDir/results.txt" file (for Human annotations)
    # HGNC ID   Approved symbol                   Previous symbols        Alias symbols   NCBI Gene ID
    # HGNC:5               A1BG                                                                      1
    # HGNC:37133       A1BG-AS1        NCRNA00181, A1BGAS, A1BG-AS             FLJ23569         503538
    # HGNC:24086           A1CF  ACF, ASP, ACF64, ACF65, APOBEC1CF                               29974
    # HGNC:6              A1S9T
    
    set NCBIandHGNCgeneDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Gene-based/NCBIandHGNCgeneID"
    if {![regexp "Human" $g_AnnotSV(organism)]} {
        ## Checked the organism: should be "Human"
        set g_AnnotSV(hpo) "" ;# used for exomiser
    } elseif {![file exists "$NCBIandHGNCgeneDir/results.txt"]} {
        ## Checked if the "results.txt" file exists
        puts "\nWARNING: No Exomiser annotations available."
        puts "...$NCBIandHGNCgeneDir/results.txt doesn't exist"
        set g_AnnotSV(hpo) "" ;# used for exomiser
    } elseif {![file exists "$NCBIandHGNCgeneDir/geneSymbol_NCBIandHGNCgeneID.tsv"]} {
        # Checked if the "geneSymbol_NCBIandHGNCgeneID.tsv" file exists
        set L_TextToWrite {"genes\tNCBI_gene_ID\tHGNC_gene_ID"}
        foreach L [LinesFromFile "$NCBIandHGNCgeneDir/results.txt"] {
            set Ls [split $L "\t"]
            
            set NCBIgeneID [lindex $Ls 4]
            if {$NCBIgeneID eq "" || [regexp -nocase "NCBI Gene ID" $NCBIgeneID]} {continue}
            
            set HGNCgeneID [lindex $Ls 0]
            
            set ApprovedSymbol [lindex $Ls 1]
            if {$ApprovedSymbol ne "" && ![info exist tmp($ApprovedSymbol)]} {set tmp($ApprovedSymbol) 1; lappend L_TextToWrite "$ApprovedSymbol\t$NCBIgeneID\t$HGNCgeneID"}
            
            regsub -all " " [lindex $Ls 2] "" PreviousSymbol
            foreach ps [split $PreviousSymbol ","] {
                if {$ps ne "" && ![info exist tmp($ps)]} {set tmp($ps) 1; lappend L_TextToWrite "$ps\t$NCBIgeneID\t$HGNCgeneID"}
            }
            
            regsub -all " " [lindex $Ls 3] "" AliasSymbol
            foreach as [split $PreviousSymbol ","] {
                if {$as ne "" && ![info exist tmp($as)]} {set tmp($as) 1; lappend L_TextToWrite "$as\t$NCBIgeneID\t$HGNCgeneID"}
            }
        }
        WriteTextInFile [join $L_TextToWrite "\n"] "$NCBIandHGNCgeneDir/geneSymbol_NCBIandHGNCgeneID.tsv"
    }
    
    return
}


proc searchforNCBIGeneID {geneName} {
    
    global g_AnnotSV
    global geneID
    
    if {![array exists geneID]} {
        
        set NCBIandHGNCgeneDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Gene-based/NCBIandHGNCgeneID"
        # Header:
        # HGNC ID   Approved symbol                   Previous symbols        Alias symbols   NCBI Gene ID
        # HGNC:5               A1BG                                                                      1
        # HGNC:37133       A1BG-AS1        NCRNA00181, A1BGAS, A1BG-AS             FLJ23569         503538
        # HGNC:24086           A1CF  ACF, ASP, ACF64, ACF65, APOBEC1CF                               29974
        # HGNC:6              A1S9T
        foreach L [LinesFromFile "$NCBIandHGNCgeneDir/results.txt"] {
            set Ls [split $L "\t"]
            if {[regexp "Approved symbol" $L]} {
                set i_approuved [lsearch -exact $Ls "Approved symbol"]
                set i_ncbi [lsearch -exact $Ls "NCBI Gene ID"]
                set i_previous [lsearch -exact $Ls "Previous symbol"]
                set i_alias [lsearch -exact $Ls "Alias symbol"]
                continue
            }
            set Approved_symbol [lindex $Ls $i_approuved]
            set NCBI_gene_ID [lindex $Ls $i_ncbi]
            if {[regexp "\\\[" $NCBI_gene_ID]} { ;# some lines have bad "NCBI_gene_ID" : "CYP4F30P   4F-se9[6:7:8]   C2orf14   100132708"
                continue
            }
            set Previous_symbol [lindex $Ls $i_previous]
            set Alias_symbol [lindex $Ls $i_alias]
            
            set geneID($Approved_symbol) $NCBI_gene_ID
            set geneID($Previous_symbol) $NCBI_gene_ID
            set geneID($Alias_symbol) $NCBI_gene_ID
        }
    }
    
    if {![info exists geneID($geneName)]} {
        set geneID($geneName) ""
    }
    
    return $geneID($geneName)
}

