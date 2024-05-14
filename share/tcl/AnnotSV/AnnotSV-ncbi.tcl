############################################################################################################
# AnnotSV 3.4.2                                                                                            #
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

proc checkNCBI {} {
    
    global g_AnnotSV
    
    # Check the existence of the "$NCBIgeneDir/results.txt" file (for Human annotations)
    # Approved symbol Alias symbol    Previous symbol NCBI gene ID
    # A1BG-AS1        FLJ23569        NCRNA00181      503538
    set NCBIgeneDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Gene-based/NCBIgeneID"
    if {![regexp "Human" $g_AnnotSV(organism)]} {
        ## Checked the organism: should be "Human"
        set g_AnnotSV(hpo) "" ;# used for exomiser
    } elseif {![file exists "$NCBIgeneDir/results.txt"]} {
        ## Checked if the "results.txt" file exists
        puts "\nWARNING: No Exomiser annotations available."
        puts "...$NCBIgeneDir/results.txt doesn't exist"
        set g_AnnotSV(hpo) "" ;# used for exomiser
    } elseif {![file exists "$NCBIgeneDir/geneSymbol_NCBIgeneID.tsv"]} {
        # Checked if the "geneSymbol_NCBIgeneID.tsv" file exists
        set L_TextToWrite {"genes\tNCBI_gene_ID"}
        foreach L [LinesFromFile "$NCBIgeneDir/results.txt"] {
            set Ls [split $L "\t"]
            set NCBIgeneID [lindex $Ls 3]
            if {$NCBIgeneID eq "" || $NCBIgeneID eq "NCBI gene ID"} {continue}
            set ApprovedSymbol [lindex $Ls 0]
            if {$ApprovedSymbol ne "" && ![info exist tmp($ApprovedSymbol)]} {set tmp($ApprovedSymbol) 1; lappend L_TextToWrite "$ApprovedSymbol\t$NCBIgeneID"}
            set AliasSymbol [lindex $Ls 1]
            if {$AliasSymbol ne "" && ![info exist tmp($AliasSymbol)]} {set tmp($AliasSymbol) 1; lappend L_TextToWrite "$AliasSymbol\t$NCBIgeneID"}
            set PreviousSymbol [lindex $Ls 2]
            if {$PreviousSymbol ne "" && ![info exist tmp($PreviousSymbol)]} {set tmp($PreviousSymbol) 1; lappend L_TextToWrite "$PreviousSymbol\t$NCBIgeneID"}
        }
        WriteTextInFile [join $L_TextToWrite "\n"] "$NCBIgeneDir/geneSymbol_NCBIgeneID.tsv"
    }
    
    return
}


proc searchforGeneID {geneName} {
    
    global g_AnnotSV
    global geneID
    
    if {![array exists geneID]} {
        
        set NCBIgeneDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Gene-based/NCBIgeneID"
        # Header:
        #Approved symbol Alias symbol    Previous symbol    	NCBI gene ID
        #A1BG                    			    	1
        #A1BG-AS1        FLJ23569        NCRNA00181         	503538
        #A1BG-AS1        FLJ23569        A1BGAS  		503538
        foreach L [LinesFromFile "$NCBIgeneDir/results.txt"] {
            set Ls [split $L "\t"]
            if {[regexp "Approved symbol" $L]} {
                set i_approuved [lsearch -exact $Ls "Approved symbol"]
                set i_ncbi [lsearch -exact $Ls "NCBI gene ID"]
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


