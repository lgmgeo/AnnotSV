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


## Define the position of some columns in an AnnotSV output line (columns requested for the ranking) 
proc SVprepareRanking {L_header} {

    global g_AnnotSV
    global g_i
    global L_Candidates

    # List of the candidate genes (given by the user)
    set L_Candidates {}
    if {$g_AnnotSV(candidateGenesFile) ne ""} {
	set L_Candidates [split [ContentFromFile $g_AnnotSV(candidateGenesFile)] " |\n|\t"]
    }

    # Check if we have all the needed information for ranking
    # Before this checking, ranking can not be done
    set g_AnnotSV(ranking) 0

    set Ls [split $L_header "\t"]

    # The user is using a VCF SV input file: the svtBEDcol is necessarily the 6th in the corresponding created annotated file
    if {[regexp ".vcf$|.vcf.gz$" $g_AnnotSV(SVinputFile)]} {
	set g_AnnotSV(svtBEDcol) 5
    } 

    # Check if we have all the needed information for ranking
    # unset g_i <=> no SV ranking  
    set g_i(gene)     [lsearch -regexp $Ls "Gene name"];    if {$g_i(gene) == -1} {unset g_i; return}  
    set g_i(full)     [lsearch -regexp $Ls "AnnotSV type"]; if {$g_i(full) == -1} {unset g_i; return}        

    set g_i(GAINtot)  [lsearch -regexp $Ls "DGV_GAIN_n_samples_tested"]; if {$g_i(GAINtot) == -1} {unset g_i; return}  
    set g_i(GAINfreq) [lsearch -regexp $Ls "DGV_GAIN_Frequency"];        if {$g_i(GAINfreq) == -1} {unset g_i; return}  
    set g_i(LOSStot)  [lsearch -regexp $Ls "DGV_LOSS_n_samples_tested"]; if {$g_i(LOSStot) == -1} {unset g_i; return}  
    set g_i(LOSSfreq) [lsearch -regexp $Ls "DGV_LOSS_Frequency"];        if {$g_i(LOSSfreq) == -1} {unset g_i; return}  

    set g_i(GDPOPMAXAF) [lsearch -regexp $Ls "GD_POPMAX_AF"];     if {$g_i(GDPOPMAXAF) == -1} {unset g_i; return}  

    set g_i(dbVar_event)   [lsearch -regexp $Ls "dbVar_event"];   if {$g_i(dbVar_event) == -1} {unset g_i; return}  
    set g_i(dbVar_status)  [lsearch -regexp $Ls "dbVar_status"];  if {$g_i(dbVar_status) == -1} {unset g_i; return}  
    
    set g_i(morbidGenes)           [lsearch -regexp $Ls "morbidGenes"];           if {$g_i(morbidGenes) == -1} {unset g_i; return}  
    set g_i(morbidGenesCandidates) [lsearch -regexp $Ls "morbidGenesCandidates"]; if {$g_i(morbidGenesCandidates) == -1} {unset g_i; return}  

    set g_i(GHgene_elite)          [lsearch -regexp $Ls "GHgene_elite"];         
    set g_i(GHgene_not_elite)      [lsearch -regexp $Ls "GHgene_not_elite"];     

    set g_i(pLI)       [lsearch -regexp $Ls "pLI"];       if {$g_i(pLI) == -1} {unset g_i; return}  
    set g_i(HI_CGscore)   [lsearch -regexp $Ls "HI_CGscore"];   if {$g_i(HI_CGscore) == -1} {unset g_i; return}  
    set g_i(TriS_CGscore) [lsearch -regexp $Ls "TriS_CGscore"]; if {$g_i(TriS_CGscore) == -1} {unset g_i; return}  

    # If we have all the needed information, ranking will be done
    set g_AnnotSV(ranking) 1
}


# If enhancer annotations is available, check the gene-enhancer relations (look at the gene, if it is in a precise list)
# Return:
#########
# - MorbidGenes            : SV overlaps an enhancer associated to a morbid gene                             (<=> ranking = 4)
# - DEL                    : SV overlaps an enhancer of a gene with a pLI > 0.9 or with HI_CGscore of 3 or 2 (<=> ranking = 4)
# - DUP                    : SV overlaps the enhancer of a gene with a TriS_CGscore of 3 or 2                (<=> ranking = 4)
# - MorbidGenesCandidates  : SV overlaps an enhancer associated to a morbid gene candidate                   (<=> ranking = 3)
# - Candidates             : SV overlaps an enhancer associated to a gene candidate (given by the user)      (<=> ranking = 3)
# - NotUsed                : No GeneHancer annotation available
proc EnhancerInformation {Ls SVtype SVtoAnn} {
    global g_AnnotSV
    global g_i

    global rankingPreparation
    global L_MorbidGenes
    global L_DEL
    global L_DUP
    global L_MorbidGenesCandidates
    global L_Candidates

    global EliteGene
    global NotEliteGene

    # Do we have enhancer information?
    ##################################
    if {$g_i(GHgene_elite) == -1 || $g_i(GHgene_not_elite) == -1} {
	return "NotUsed"
    }


    # Creation of different genes list to prepare the ranking:
    ##########################################################
    if {![info exists rankingPreparation]} {
	set rankingPreparation 1

	# List of the morbid genes:
	set omimDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Genes-based/OMIM"
	set MorbidGenesFileFormattedGzip [glob -nocomplain "$omimDir/*_morbidGenes.tsv.gz"]
	set L_MorbidGenes {}
	if {[file exists $MorbidGenesFileFormattedGzip]} {
	    foreach L [LinesFromGZFile $MorbidGenesFileFormattedGzip] {
		lappend L_MorbidGenes [lindex $L 0]
	    }
	}
	# L_DEL: List of genes with a pLI > 0.9 or with HI_CGscore of 3 or 2
	# L_DUP: List of genes with a TriS_CGscore of 3 or 2
	set ExACdir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Genes-based/ExAC"
	set ExACfile [glob -nocomplain "$ExACdir/*_GeneIntolerance.pLI-Zscore.annotations.tsv.gz"]
	set L_DEL {}
	if {[file exists $ExACfile]} {
	    foreach L [LinesFromGZFile $ExACfile] {
		set Ls [split $L "\t"]
		set pLI [lindex $Ls 3] 
		if {$pLI > 0.9} {
		    lappend L_DEL [lindex $Ls 0]
		}
	    }
	}
	set ClinGenDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Genes-based/ClinGen"
	set CGfile [glob -nocomplain "$ClinGenDir/*_ClinGenAnnotations.tsv.gz"]
	set L_DUP {}
	if {[file exists $CGfile]} {
	    foreach L [LinesFromGZFile $CGfile] {
		set Ls [split $L "\t"]
		set HI_CGscore [lindex $Ls 1] 
		if {$HI_CGscore eq "2" || $HI_CGscore eq "3"} {
		    lappend L_DEL [lindex $Ls 0]
		}
		set TriS_CGscore [lindex $Ls 2] 
		if {$TriS_CGscore eq "2" || $TriS_CGscore eq "3"} {
		    lappend L_DUP [lindex $Ls 0]
		}
 
	    }
	}
	# List of the morbid genes candidates:
	set MorbidGenesCandidatesFileFormattedGzip [glob -nocomplain "$omimDir/*_morbidGenescandidates.tsv.gz"]
	set L_MorbidGenesCandidates {}
	if {[file exists $MorbidGenesCandidatesFileFormattedGzip]} {
	    foreach L [LinesFromGZFile $MorbidGenesCandidatesFileFormattedGzip] {
		lappend L_MorbidGenesCandidates [lindex $L 0]
	    }
	}
    }


    # Listing of the elite_genes and not_elite_genes overlapped by the SV
    #####################################################################
    set L_enhancersAssociatedGenes {}
    set EG "[lindex $Ls $g_i(GHgene_elite)]"
    if {[regexp "\\.\\.\\.$" $EG] && [info exists EliteGene($SVtoAnn)]} {       ;# AnnotSV restrict the number of overlapping reported features to 20. Keep back all the genes values
	set EG $EliteGene($SVtoAnn)                                             ;# ok, there is no ";"
    } else {
	set EG [split $EG ";"]
    }
    set NEG "[lindex $Ls $g_i(GHgene_not_elite)]"
    if {[regexp "\\.\\.\\.$" $NEG] && [info exists NotEliteGene($SVtoAnn)]} {   ;# AnnotSV restrict the number of overlapping reported features to 20. Keep back all the genes values
	set NEG $NotEliteGene($SVtoAnn)                                         ;# ok, there is no ";"
    } else {
	set NEG [split $NEG ";"]
    } 
    set L_enhancersAssociatedGenes [lsort -unique [list {*}$EG {*}$NEG]]


    # Check if the SV overlaps an enhancer associated to a morbid gene 
    ##################################################################
    set thegenes "MorbidGenes"
    foreach g $L_enhancersAssociatedGenes {
	if {[lsearch -exact $L_MorbidGenes $g] ne -1} {
	    lappend thegenes $g
	}
    }	
    if {$thegenes ne "MorbidGenes"} {return "$thegenes"}


    # Check if the SV overlaps an enhancer associated to a gene with a pLI > 0.9 or with HI_CGscore of 3 or 2
    #########################################################################################################
    # Check for a del:
    set thegenes "DEL"
    if {[regexp -nocase "Del|Loss|<CN0>" $SVtype]} {
	foreach g $L_enhancersAssociatedGenes {
	    if {[lsearch -exact $L_DEL $g] ne -1} {
		lappend thegenes $g
	    }
	}	    
    }
    if {$thegenes ne "DEL"} {return "$thegenes"}

    # Check if the SV overlaps an enhancer associated to a gene with a TriS_CGscore of 3 or 2
    #########################################################################################
    # Check for a dup:
    set thegenes "DUP"
    if {[regexp -nocase "Dup|Gain|Multiplication|<CN\[2-9\]" $SVtype]} {
	foreach g $L_enhancersAssociatedGenes {
	    if {[lsearch -exact $L_DUP $g] ne -1} {
		lappend thegenes $g
	    }
	}
    }	    
    if {$thegenes ne "DUP"} {return "$thegenes"}

    # Check if the SV overlaps an enhancer associated to a morbid gene candidate
    ############################################################################
    set thegenes "MorbidGenesCandidates"
    foreach g $L_enhancersAssociatedGenes {
	if {[lsearch -exact $L_MorbidGenesCandidates $g] ne -1} {
	    lappend thegenes $g
	}
    }	
    if {$thegenes ne "MorbidGenesCandidates"} {return "$thegenes"}

    # Check if the SV overlaps an enhancer associated to a candidate gene (given by the user)
    #########################################################################################
    set thegenes "Candidates"
    foreach g $L_enhancersAssociatedGenes {
	if {[lsearch -exact $L_Candidates $g] ne -1} {
	    lappend thegenes $g
	}
    }	
    if {$thegenes ne "Candidates"} {return "$thegenes"}
}


## Rank the pathogenicity of the different SV as follows:
#########################################################
## category 1 = benign
##           > 70% SV overlapped with a benign SV + doesn't overlap with i) a morbid gene, ii) a morbid gene candidate and iii) a candidate gene
## category 2 = likely benign
##           < 70% SV overlapped with a benign SV + doesn't overlap with i) a morbid gene, ii) a morbid gene candidate and iii) a candidate gene (or their enhancers)
## category 3 = VOUS (variant of unknown significance)
##           SV that overlaps i) a morbid gene candidate  (or its enhancer) and ii) a candidate gene  (or its enhancer) (with at least 1bp)
## category 4 = likely pathogenic 
##           SV that overlaps a morbid gene (or its enhancer) (with at least 1bp)
##           or for a del: SV that overlap a gene (or its enhancer) with a pLI > 0.9 or with HI_CGscore of 3 or 2
##           or for a dup: SV that overlap a gene (or its enhancer) TriS_CGscore of 3 or 2
## category 5 = pathogenic
##           SV that overlaps a pathogenic SV (with at least 1bp)

## ClinGen HI_CGscore and TriS_CGscore explanations:
####################################################
##   Rating	Possible Clinical Interpretation
##   ------     --------------------------------
##   3	        Sufficient evidence for dosage pathogenicity
##   2	        Some evidence for dosage pathogenicity
##   1	        Little evidence for dosage pathogenicity
##   0	        No evidence for dosage pathogenicity
##   40         Evidence suggests the gene is not dosage sensitive
##
## HI = Haploinsufficiency  TriS = Triplosensitivity
proc SVranking {L_annotations ref alt} {

    global g_AnnotSV
    global g_i
    global L_Candidates
    global L_rankingExplanations

    # Check if we have enougth information to do the ranking:
    #########################################################
    set ranking ""
    if {$g_AnnotSV(svtBEDcol) eq -1} {return $ranking}
    if {!$g_AnnotSV(ranking)} {return $ranking}


    # Ranking!!
    ###########

    set Ls [split $L_annotations "\t"]
    set SVtype [lindex $Ls $g_AnnotSV(svtBEDcol)]   
    set SVtoAnn [join [lrange $Ls 1 3] ","]
    set enhancer [EnhancerInformation $Ls $SVtype $SVtoAnn]
    set AnnotSV_ID [settingOfTheAnnotSVID "[join [lrange $Ls 1 3] "_"]_$SVtype" "$ref" "$alt"]
    set AnnotSVtype [lindex $Ls $g_i(full)]
    set genes [lindex $Ls $g_i(gene)]

    ## category 5 = pathogenic
    ##              SV that overlap a pathogenic SV (with at least 1bp)
    ###################################################################
    set dbVar_status [lindex $Ls $g_i(dbVar_status)]
    if {$dbVar_status ne ""} {
	set dbVar_event [lindex $Ls $g_i(dbVar_event)]
	if {[regexp -nocase "Del|Loss|<CN0>" $SVtype] && [regexp -nocase "Del|Loss|<CN0>" $dbVar_event]} {
	    set ranking "5"	
	    lappend L_rankingExplanations "$AnnotSV_ID\t$AnnotSVtype\t$genes\t$ranking\tLOSS: pathogenic SV overlapped"
	    return $ranking
	}
	if {[regexp -nocase "Dup|Gain|Multiplication|<CN\[2-9\]" $SVtype] && [regexp -nocase "Dup|Gain|Multiplication|<CN\[2-9\]" $dbVar_event]} {
	    set ranking "5"	
	    lappend L_rankingExplanations "$AnnotSV_ID\t$AnnotSVtype\t$genes\t$ranking\tGAIN: pathogenic SV overlapped"
	    return $ranking
	}
    }

    ## category 4 = likely pathogenic 
    ##           SV that overlaps a morbid gene (or its enhancer) (with at least 1bp)
    ##           or for a del: SV that overlaps a gene (or its enhancer) with a pLI > 0.9 or with HI_CGscore of 3 or 2
    ##           or for a dup: SV that overlaps a gene (or its enhancer) with a TriS_CGscore of 3 or 2
    ###################################################################

    # Check if a SV overlaps a morbid gene 
    set morbidGenes [lindex $Ls $g_i(morbidGenes)]
    if {[regexp "yes" $morbidGenes]} {
	set ranking "4"	
	lappend L_rankingExplanations "$AnnotSV_ID\t$AnnotSVtype\t$genes\t$ranking\tmorbid gene overlapped"
	return $ranking
    }
    # Check if a SV overlaps an enhancer associated to a morbid gene 
    if {[lindex $enhancer 0] eq "MorbidGenes"} {
	set ranking "4"	
	lappend L_rankingExplanations "$AnnotSV_ID\t$AnnotSVtype\t$genes\t$ranking\tenhancer of a morbid gene overlapped ([lrange $enhancer 1 end])"
	return $ranking
    }
    
    # Check for a del:
    regsub "," [lindex $Ls $g_i(pLI)] "." pLI; # needed with -metrics=fr 
    set HI_CGscore [lindex $Ls $g_i(HI_CGscore)]
    if {[regexp -nocase "Del|Loss|<CN0>" $SVtype]} {
	# Check SV that overlap a gene with a pLI > 0.9 or with HI_CGscore of 3 or 2
	if {$pLI > 0.9} {    ; # {"">0.9} is false; code ok
	    set ranking "4"	
	    lappend L_rankingExplanations "$AnnotSV_ID\t$AnnotSVtype\t$genes\t$ranking\tLOSS: pLI ($pLI) > 0.9"
	    return $ranking
	} elseif {$HI_CGscore eq 3 || $HI_CGscore eq 2} {   
	    set ranking "4"	
	    lappend L_rankingExplanations "$AnnotSV_ID\t$AnnotSVtype\t$genes\t$ranking\tLOSS: HI_CGscore = $HI_CGscore"
	    return $ranking
	}
	# Check SV that overlap the enhancer of a gene with a pLI > 0.9 or with HI_CGscore of 3 or 2
	if {[lindex $enhancer 0] eq "DEL"} {
	    set ranking "4"	
	    lappend L_rankingExplanations "$AnnotSV_ID\t$AnnotSVtype\t$genes\t$ranking\tLOSS: overlap the enhancer of a gene with a pLI > 0.9 or with a HI_CGscore of 3 or 2 ([lrange $enhancer 1 end])"
	    return $ranking
	}
    }
    
    # Check for a dup:
    set TriS_CGscore [lindex $Ls $g_i(TriS_CGscore)]
    if {[regexp -nocase "Dup|Gain|Multiplication|<CN\[2-9\]" $SVtype]} {
	# Check SV that overlap a gene TriS_CGscore of 3 or 2
	if {$TriS_CGscore eq 3 || $TriS_CGscore eq 2} {
	    set ranking "4"	
	    lappend L_rankingExplanations "$AnnotSV_ID\t$AnnotSVtype\t$genes\t$ranking\tGAIN: TriS_CGscore = $TriS_CGscore"
	    return $ranking
	}
	# Check SV that overlap the enhancer of a gene TriS_CGscore of 3 or 2
	if {[lindex $enhancer 0] eq "DUP"} {
	    set ranking "4"	
	    lappend L_rankingExplanations "$AnnotSV_ID\t$AnnotSVtype\t$genes\t$ranking\tGAIN: overlap the enhancer of a gene with a TriS_CGscore of 3 or 2 ([lrange $enhancer 1 end])"
	    return $ranking
	}
    }
    
    
    ## category 3 = VOUS (variant of unknown significance)
    ##           SV that overlap an enhancer or a CDS from i) a morbid gene candidate and ii) a candidate gene (with at least 1bp)
    ###################################################################

    # Check if a SV overlaps a morbid gene candidate
    set morbidGenesCandidates [lindex $Ls $g_i(morbidGenesCandidates)]
    if {[regexp "yes" $morbidGenesCandidates]} {
	set ranking "3"	
	lappend L_rankingExplanations "$AnnotSV_ID\t$AnnotSVtype\t$genes\t$ranking\tmorbid gene candidate overlapped"
	return $ranking
    }
    # Check if a SV overlaps an enhancer associated to a morbid gene candidate
    if {[lindex $enhancer 0] eq "MorbidGenesCandidates"} {
	set ranking "3"	
	lappend L_rankingExplanations "$AnnotSV_ID\t$AnnotSVtype\t$genes\t$ranking\tenhancer of a morbid gene candidate overlapped ([lrange $enhancer 1 end])"
	return $ranking
    }

    if {$g_AnnotSV(candidateGenesFile) ne ""} {
	# Check if a SV overlaps a CDS from a candidate gene (with at least 1bp)
	foreach gene [split $genes "/"] {
	    if {$gene ne "" && [lsearch -exact $L_Candidates $gene] ne -1} {
		set ranking "3"	
		lappend L_rankingExplanations "$AnnotSV_ID\t$AnnotSVtype\t$genes\t$ranking\tcandidate gene overlapped ($gene)"
		return $ranking
	    }
	}
	# Check if a SV overlaps the enhancer from a candidate gene (with at least 1bp)
	if {[lindex $enhancer 0] eq "Candidates"} {
	    set ranking "3"	
	    lappend L_rankingExplanations "$AnnotSV_ID\t$AnnotSVtype\t$genes\t$ranking\tenhancer of a candidate gene overlapped ([lrange $enhancer 1 end])"
	    return $ranking
	}
    }


    ## category 1 = benign 
    ##           > 70% SV overlapped with a benign SV
    ##           > 70% SV overlapped with a frequent SV from gnomAD (GD_POPMAX_AF)
    ##           + doesn't contain CDS from i) a morbid gene, ii) a morbid gene candidate and iii) a candidate gene
    ###################################################################
    set GAINtot [lindex $Ls $g_i(GAINtot)] 
    regsub ","  [lindex $Ls $g_i(GAINfreq)] "." GAINfreq; # needed with -metrics=fr 
    set LOSStot [lindex $Ls $g_i(LOSStot)]
    regsub ","  [lindex $Ls $g_i(LOSSfreq)] "." LOSSfreq; # needed with -metrics=fr 
    regsub ","  [lindex $Ls $g_i(GDPOPMAXAF)] "." GDPOPMAXAF; # needed with -metrics=fr 

    if {[regexp -nocase "Del|Loss|<CN0>" $SVtype]} {
	if {$LOSStot > $g_AnnotSV(minTotalNumber) && $LOSSfreq > 0.01} {
	    set ranking "1"	
	    lappend L_rankingExplanations "$AnnotSV_ID\t$AnnotSVtype\t$genes\t$ranking\t> $g_AnnotSV(overlap)% SV overlapped with a frequent SV + doesn't contain CDS from i) a morbid gene, ii) a morbid gene candidate and iii) a candidate gene"
	    return $ranking
	}
    } elseif {[regexp -nocase "Dup|Gain|Multiplication|<CN\[2-9\]" $SVtype]} {
	if {$GAINtot > $g_AnnotSV(minTotalNumber) && $GAINfreq > 0.01} {
	    set ranking "1"	
	    lappend L_rankingExplanations "$AnnotSV_ID\t$AnnotSVtype\t$genes\t$ranking\t> $g_AnnotSV(overlap)% SV overlapped with a frequent SV + doesn't contain CDS from i) a morbid gene, ii) a morbid gene candidate and iii) a candidate gene"
	    return $ranking
	}
    }
    if {$GDPOPMAXAF > 0.01} {
	set ranking "1"	
	lappend L_rankingExplanations "$AnnotSV_ID\t$AnnotSVtype\t$genes\t$ranking\t> $g_AnnotSV(overlap)% SV overlapped with a frequent SV + doesn't contain CDS from i) a morbid gene, ii) a morbid gene candidate and iii) a candidate gene"
	return $ranking
    }
    
    ## category 2 = likely benign
    ##           < 70% SV overlapped with a benign SV + doesn't contain CDS from i) a morbid gene, ii) a morbid gene candidate and iii) a candidate gene
    ###################################################################    
    set ranking "2"	
    lappend L_rankingExplanations "$AnnotSV_ID\t$AnnotSVtype\t$genes\t$ranking\t< $g_AnnotSV(overlap)% SV overlapped with a benign SV + doesn't contain CDS from i) a morbid gene, ii) a morbid gene candidate and iii) a candidate gene"
    return $ranking
}
