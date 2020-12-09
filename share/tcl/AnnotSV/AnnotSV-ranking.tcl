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


## Define the position of some columns in an AnnotSV output line (columns requested for the ranking) 
proc SVprepareRanking {L_header} {

    global g_AnnotSV
    global g_i


    # Check if we have all the needed information for ranking
    # Before this checking, ranking can not be done
    set g_AnnotSV(ranking) 0
    set Ls [split $L_header "\t"]

    # The user is using a VCF SV input file: the svtTSVcol is necessarily the 6th in the corresponding output annotated file.tsv
    if {[regexp ".vcf$|.vcf.gz$" $g_AnnotSV(SVinputFile)]} {
	set g_AnnotSV(svtTSVcol) 5
	set g_AnnotSV(svtBEDcol) 3
    }

    # "SV type" is compulsory for the ranking
    if {$g_AnnotSV(svtTSVcol) eq -1} {return}

    # Check if we have all the needed information for ranking
    # unset g_i <=> no SV ranking  
    set g_i(gene)     [lsearch -regexp $Ls "Gene name"];      if {$g_i(gene) == -1} {unset g_i; return}  
    set g_i(NbGenes)  [lsearch -regexp $Ls "Number of gene"]; if {$g_i(NbGenes) == -1} {unset g_i; return}  
    set g_i(full)     [lsearch -regexp $Ls "AnnotSV type"];   if {$g_i(full) == -1} {unset g_i; return}        

    set g_i(REgene)   [lsearch -regexp $Ls "RE_gene"]; if {$g_i(REgene) == -1} {unset g_i; return}
    
    set g_i(Ploss)      [lsearch -regexp $Ls "P_loss_coord"];  if {$g_i(Ploss) == -1} {unset g_i; return}  
    set g_i(Pgain)      [lsearch -regexp $Ls "P_gain_coord"];  if {$g_i(Pgain) == -1} {unset g_i; return}  
    set g_i(Bloss)      [lsearch -regexp $Ls "B_loss_coord"];  if {$g_i(Bloss) == -1} {unset g_i; return}  
    set g_i(Bgain)      [lsearch -regexp $Ls "B_gain_coord"];  if {$g_i(Bgain) == -1} {unset g_i; return}  
    set g_i(Psnvindel)  [lsearch -regexp $Ls "P_snvindel_nb"]; if {$g_i(Psnvindel) == -1} {unset g_i; return} 

    set g_i(HI) [lsearch -regexp $Ls "HI"];     if {$g_i(HI) == -1} {unset g_i; return}  
    set g_i(TS) [lsearch -regexp $Ls "TS"];     if {$g_i(TS) == -1} {unset g_i; return}  

    set g_i(pLI)       [lsearch -regexp $Ls "pLI_gnomAD"];    if {$g_i(pLI) == -1} {unset g_i; return}  
    set g_i(loeuf)     [lsearch -regexp $Ls "LOEUF_bin"];     if {$g_i(loeuf) == -1} {unset g_i; return}  
    set g_i(HIpercent) [lsearch -regexp $Ls "HI_DDDpercent"]; if {$g_i(HIpercent) == -1} {unset g_i; return}  
   
    set g_i(exomiser)  [lsearch -regexp $Ls "EXOMISER_GENE_PHENO_SCORE"];  

    set g_i(morbidGenes) [lsearch -regexp $Ls "morbidGenes"]; if {$g_i(morbidGenes) == -1} {unset g_i; return}   

    set g_i(location)    [lsearch -regexp $Ls "location"];               if {$g_i(location) == -1} {unset g_i; return}  
    set g_i(location2)   [lsearch -regexp $Ls "location2"];              if {$g_i(location2) == -1} {unset g_i; return} 
    set g_i(CDSlength)   [lsearch -regexp $Ls "overlapped CDS length"];  if {$g_i(CDSlength) == -1} {unset g_i; return}
    set g_i(CDSpercent)  [lsearch -regexp $Ls "overlapped CDS percent"]; if {$g_i(CDSpercent) == -1} {unset g_i; return} 
    set g_i(frameshift)  [lsearch -regexp $Ls "frameshift"];             if {$g_i(frameshift) == -1} {unset g_i; return} 
    set g_i(NbExons)     [lsearch -regexp $Ls "Number of exons"];        if {$g_i(NbExons) == -1} {unset g_i; return} 

    # If we have all the needed information, ranking will be done
    set g_AnnotSV(ranking) 1
}



## Rank the pathogenicity of the "DEL" SV following the ACMG classification.
## Creation of 2 global variables:
## - g_rankingScore($AnnotSV_ID)
## - g_rankingExplanations($AnnotSV_ID) 
proc SVrankingLoss {L_annotations} {

    global g_AnnotSV
    global g_i
    global g_rankingExplanations
    global g_rankingScore

   
    set Ls [split $L_annotations "\t"]
    set AnnotSV_ID "[lindex $Ls 0]" 

    # In case of SV redundancy in the input file
    if {[info exists g_rankingScore($AnnotSV_ID)]} {return}
    
    set AnnotSVtype [lindex $Ls $g_i(full)]    
    if {$AnnotSVtype eq "full"} {
	
	set NbGenes   [lindex $Ls $g_i(NbGenes)]
	set REgene    [lindex $Ls $g_i(REgene)]
	set Ploss     [lindex $Ls $g_i(Ploss)]
	set Bloss     [lindex $Ls $g_i(Bloss)]
	
	## Section 1: Initial assessment of genomic content
	####################################################################################################################
	if {$REgene eq "" && $NbGenes eq "0"} {
	    # 1B. Does NOT contain protein-coding or any known functionally important elements (- 0.60)
	    set g_rankingScore($AnnotSV_ID) "-0.60"
	    set g_rankingExplanations($AnnotSV_ID) "1B "
	} else {
	    # 1A. Contains protein-coding or other known functionally important elements (+ 0)
	    set g_rankingScore($AnnotSV_ID) "0"
	    set g_rankingExplanations($AnnotSV_ID) "1A "
	}
	
	## Section 2: Overlap with established/predicted haploinsufficiency (HI) or established benign genes/genomic regions
	##            (Skip to section 3 if your copy-number loss DOES NOT overlap these types of genes/regions)
	#####################################################################################################################
	if {$Ploss ne ""} {
	    # 2A. Complete overlap of an established HI gene/genomic region (+ 1)
	    set g_rankingExplanations($AnnotSV_ID,2A) ""
	} elseif {$Bloss ne ""} {
	    # 2F. Completely contained within an established benign CNV region (- 1)
	    set g_rankingExplanations($AnnotSV_ID,2F) ""
	} elseif {0} {		
	    # 2G. Overlaps an established benign CNV, but includes additional genomic material (+ 0)
	    # Vero to improve
	    set g_rankingExplanations($AnnotSV_ID,2G) ""
	}
	
	## Section 3: Evaluation of gene number
	####################################################################################################################
	if {$NbGenes <= 24} {
	    # 3A. 0-24 genes (+ 0)
	    if {$NbGenes < 2} {
		set g_rankingExplanations($AnnotSV_ID,3A) "+ 3A ($NbGenes gene) "
	    } else {
		set g_rankingExplanations($AnnotSV_ID,3A) "+ 3A ($NbGenes genes) "
	    }
	} elseif {$NbGenes <= 35} {
	    # 3B. 25-34 genes (+ 0.45)
	    set g_rankingExplanations($AnnotSV_ID,3B) "+ 3B ($NbGenes genes) "
	} else {
	    # 3C. 35+ genes (+ 0.90)
	    set g_rankingExplanations($AnnotSV_ID,3C) "+ 3C ($NbGenes genes) "
	} 
	
	## Section 4: Detailed evaluation of genomic content using cases from published literature, public databases, and/or internal lab data
	##            (Skip to section 5 if either your CNV overlapped with an established HI gene/region in section 2,
	##             OR there have been no reports associating either the CNV or any genes within the CNV with human phenotypes
	##             caused by loss of function [LOF] or copy-number loss)
	####################################################################################################################
	
	
	## Section 5: Evaluation of inheritance pattern/family history for patient being studied			
	####################################################################################################################
	if {$REgene ne ""} {
	    set L_infos [split $REgene ";"]
	    set bestEx "0"
	    foreach infos $L_infos {
		if {[regexp "(.+?) \\\(.*?EX=(.+)\\\)$" $infos match g ex]} {
		    regsub "," $ex "." ex
		    if {$ex > $bestEx} {
			set bestEx $ex
			set bestG $g
		    }
		}
	    }
	    if {$bestEx >= 0.7} {
		# 5H. Inheritance information is unavailable or uninformative. The patient phenotype is highly specific and consistent with what has been described in similar cases (+ 0.30)
		lappend g_rankingExplanations($AnnotSV_ID,5H) "RE:$bestG"
	    } elseif {$bestEx >= 0.5} {
		# 5G. Inheritance information is unavailable or uninformative. The patient phenotype is nonspecific, but is consistent with what has been described in similar cases (+ 0.10)
		lappend g_rankingExplanations($AnnotSV_ID,5G) "RE:$bestG"
	    }
	}
	    
    } else {
	# Split lines
	
	## Section 2: Overlap with established/predicted haploinsufficiency (HI) or established benign genes/genomic regions
	##            (Skip to section 3 if your copy-number loss DOES NOT overlap these types of genes/regions)
	#####################################################################################################################
	set HI [lindex $Ls $g_i(HI)]
	set morbidGenes [lindex $Ls $g_i(morbidGenes)]
	set gene       [lindex $Ls $g_i(gene)]
	if {$HI eq "3" || $morbidGenes eq "yes"} {
	    set location   [lindex $Ls $g_i(location)]
	    set location2  [lindex $Ls $g_i(location2)]
	    set CDSlength  [lindex $Ls $g_i(CDSlength)]
	    set CDSpercent [lindex $Ls $g_i(CDSpercent)]
	    regsub "," $CDSpercent "." CDSpercent
	    set Psnvindel  [lindex $Ls $g_i(Psnvindel)]
	    set NbExons    [lindex $Ls $g_i(NbExons)]

	    if {[regexp "^5'UTR" $location2] && ![regexp "3'UTR$" $location2]} {
		# 2C. Partial overlap with the 5’ end of an established HI gene (3’ end of the gene not involved)…
		if {$CDSlength ne "0"} {
		    # 2C-1. …and coding sequence is involved (+ 0.90)
		    lappend g_rankingExplanations($AnnotSV_ID,2C-1) "$gene"
		} else {
		    # 2C-2. …and only the 5’ UTR is involved (+ 0)
		    lappend g_rankingExplanations($AnnotSV_ID,2C-2) "$gene"		    
		}
	    } elseif {[regexp "3'UTR$" $location2] && ![regexp "^5'UTR" $location2]} {
		# 2D. Partial overlap with the 3’ end of an established HI gene (5’ end of the gene not involved)…
		if {![regexp "exon" $location]} {
		    # 2D-1. …and only the 3’ untranslated region is involved (+ 0)
		    lappend g_rankingExplanations($AnnotSV_ID,2D-1) "$gene"		    
		} elseif {$location eq "^exon${NbExons}-3'UTR" && $Psnvindel ne ""} {
		    # 2D-2. …and only the last exon is involved. Other established pathogenic snv/indel have been reported in the SV (+ 0.90)
		    # Vero to improve: Other established pathogenic snv/indel have been reported in this exon
		    lappend g_rankingExplanations($AnnotSV_ID,2D-2) "$gene"
		} elseif {$location eq "^exon${NbExons}-3'UTR"} {
		    # 2D-3. …and only the last exon is involved. No other established pathogenic snv/indel have been reported in the SV (+ 0.30)
		    # Vero to improve: No other established pathogenic snv/indel have been reported in this exon
		    lappend g_rankingExplanations($AnnotSV_ID,2D-3) "$gene"
		} elseif {[regexp "exon" $location]} {	    
		    # 2D-4. …and it includes other exons in addition to the last exon. Nonsense-mediated decay is expected to occur (+ 0.90)
		    lappend g_rankingExplanations($AnnotSV_ID,2D-4) "$gene"		    
		}
	    } elseif {![regexp "txStart|txEnd" $location]} {
		# 2E. Both breakpoints are within the same gene (intragenic CNV; gene-level sequence variant).
		set frameshift [lindex $Ls $g_i(frameshift)]
		if {$frameshift eq "yes"} {
		    # 2E-1. frameshift (+ 0.90)
		    lappend g_rankingExplanations($AnnotSV_ID,2E-1) "$gene"		    
		} elseif {[regexp "exon" $location] && $Psnvindel ne "" && $CDSpercent >= 10} {
		    # 2E-2. exon(s) overlapped AND pathogenic snv/indel overlapped AND SV removes > 10% of protein (+ 0.45)
		    lappend g_rankingExplanations($AnnotSV_ID,2E-2) "$gene"		    
		} elseif {[regexp "exon" $location] && $Psnvindel ne "" && $CDSpercent < 10} {
		    # 2E-3. exon(s) overlapped AND pathogenic snv/indel overlapped AND SV removes < 10% of protein (+ 0.30)
		    lappend g_rankingExplanations($AnnotSV_ID,2E-3) "$gene"		    
		} 
	    }
	} else {
	    set pLI        [lindex $Ls $g_i(pLI)]
	    set loeuf      [lindex $Ls $g_i(loeuf)]
	    set HIpercent  [lindex $Ls $g_i(HIpercent)]

	    set i 0
	    if {$pLI >= 0.9} {incr i}
	    if {$HIpercent <= 10} {incr i}
	    if {$loeuf < 2} {incr i} ;# loeuf = 0 or loeuf = 1
	    if {$i >= 2} {
		# 2H. Two or more HI predictors suggest that AT LEAST ONE gene in the interval is HI (+ 0.15)
		lappend g_rankingExplanations($AnnotSV_ID,2H) "$gene"
	    } 
	}

	## Section 5: Evaluation of inheritance pattern/family history for patient being studied			
	####################################################################################################################
	set exomiser   [lindex $Ls $g_i(exomiser)]
	regsub "," $exomiser "." exomiser
	if {$exomiser >= 0.7} {
	    # 5H. Inheritance information is unavailable or uninformative. The patient phenotype is highly specific and consistent with what has been described in similar cases (+ 0.30)
	    lappend g_rankingExplanations($AnnotSV_ID,5H) "$gene"
	} elseif {$exomiser >= 0.5} {
	    # 5G. Inheritance information is unavailable or uninformative. The patient phenotype is nonspecific, but is consistent with what has been described in similar cases (+ 0.10)
	    lappend g_rankingExplanations($AnnotSV_ID,5G) "$gene"
	}
    }
   
    return
}

proc achieveSVrankingLoss {AnnotSV_ID} {

    global g_rankingScore
    global g_rankingExplanations
    global g_achieve
    
    # In case of SV redundancy in the input file
    if {[info exists g_achieve($AnnotSV_ID)]} {return}
    set g_achieve($AnnotSV_ID) ""
    
    if {![info exists g_rankingScore($AnnotSV_ID)]} {set g_rankingScore($AnnotSV_ID) ""; return}

    ## Section 2: Overlap with established/predicted haploinsufficiency (HI) or established benign genes/genomic regions
    ##            (Skip to section 3 if your copy-number loss DOES NOT overlap these types of genes/regions)
    #####################################################################################################################
    # Add the higher score of the section 2
    if {[info exists g_rankingExplanations($AnnotSV_ID,2A)]} {
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+1}]
	append g_rankingExplanations($AnnotSV_ID) "+ 2A (cf P_loss_source) "
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2C-1)]} {
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.90}]
	append g_rankingExplanations($AnnotSV_ID) "+ 2C-1 ([join $g_rankingExplanations($AnnotSV_ID,2C-1) "/"]) "
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2D-2)]} {
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.90}]
	append g_rankingExplanations($AnnotSV_ID) "+ 2D-2 ([join $g_rankingExplanations($AnnotSV_ID,2D-2) "/"]) "
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2D-4)]} {
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.90}]
	append g_rankingExplanations($AnnotSV_ID) "+ 2D-4 ([join $g_rankingExplanations($AnnotSV_ID,2D-4) "/"]) "
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2E-1)]} {
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.90}]
	append g_rankingExplanations($AnnotSV_ID) "+ 2E-1 ([join $g_rankingExplanations($AnnotSV_ID,2E-1) "/"]) "
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2E-2)]} {
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.45}]
	append g_rankingExplanations($AnnotSV_ID) "+ 2E-2 ([join $g_rankingExplanations($AnnotSV_ID,2E-2) "/"]) "
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2D-3)]} {
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.30}]
	append g_rankingExplanations($AnnotSV_ID) "+ 2D-3 ([join $g_rankingExplanations($AnnotSV_ID,2D-3) "/"]) "
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2E-3)]} {
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.30}]
	append g_rankingExplanations($AnnotSV_ID) "+ 2E-3 ([join $g_rankingExplanations($AnnotSV_ID,2E-3) "/"]) "
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2H)]} {
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.15}]
	append g_rankingExplanations($AnnotSV_ID) "+ 2H ([join $g_rankingExplanations($AnnotSV_ID,2H) "/"]) "
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2C-2)]} {
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0}]
	append g_rankingExplanations($AnnotSV_ID) "+ 2C-2 ([join $g_rankingExplanations($AnnotSV_ID,2C-2) "/"]) "
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2D-1)]} {
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0}]
	append g_rankingExplanations($AnnotSV_ID) "+ 2D-1 ([join $g_rankingExplanations($AnnotSV_ID,2D-1) "/"]) "
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2G)]} {
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0}]
	append g_rankingExplanations($AnnotSV_ID) "+ 2G "
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2F)]} {
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)-1}]
	append g_rankingExplanations($AnnotSV_ID) "+ 2F (cf B_loss_source) "
    }

    ## Section 3: Evaluation of gene number
    ####################################################################################################################
    # Add the higher score of the section 3
    if {[info exists g_rankingExplanations($AnnotSV_ID,3C)]} {
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.90}]
	append g_rankingExplanations($AnnotSV_ID) "$g_rankingExplanations($AnnotSV_ID,3C)"
     } elseif {[info exists g_rankingExplanations($AnnotSV_ID,3B)]} {
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.45}]
	append g_rankingExplanations($AnnotSV_ID) "$g_rankingExplanations($AnnotSV_ID,3B)"
     } elseif {[info exists g_rankingExplanations($AnnotSV_ID,3A)]} {
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.00}]
	append g_rankingExplanations($AnnotSV_ID) "$g_rankingExplanations($AnnotSV_ID,3A)"
     } 
    
    ## Section 4: Detailed evaluation of genomic content using cases from published literature, public databases, and/or internal lab data
    ##            (Skip to section 5 if either your CNV overlapped with an established HI gene/region in section 2,
    ##             OR there have been no reports associating either the CNV or any genes within the CNV with human phenotypes
    ##             caused by loss of function [LOF] or copy-number loss)
    ####################################################################################################################
    # Add the higher score of the section 4

    
    ## Section 5: Evaluation of inheritance pattern/family history for patient being studied			
    ####################################################################################################################
    # Add the higher score of the section 5
    if {[info exists g_rankingExplanations($AnnotSV_ID,5H)]} {
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.30}]
	append g_rankingExplanations($AnnotSV_ID) "+ 5H ([join $g_rankingExplanations($AnnotSV_ID,5H) "/"]) "
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,5G)]} {
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.10}]
	append g_rankingExplanations($AnnotSV_ID) "+ 5G ([join $g_rankingExplanations($AnnotSV_ID,5G) "/"]) "
    }

    return
}

    
## Rank the pathogenicity of the "DUP" SV following the ACMG classification.
## Creation of 2 global variables:
## - g_rankingScore($AnnotSV_ID)
## - g_rankingExplanations($AnnotSV_ID) 
proc SVrankingGain {L_annotations} {

    global g_AnnotSV
    global g_i
    global g_rankingExplanations
    global g_rankingScore

   
    set Ls [split $L_annotations "\t"]
    set AnnotSV_ID "[lindex $Ls 0]" 

    set AnnotSVtype [lindex $Ls $g_i(full)]
    
    if {$AnnotSVtype eq "full"} {
	
	set NbGenes   [lindex $Ls $g_i(NbGenes)]
	set REgene    [lindex $Ls $g_i(REgene)]
	set Pgain     [lindex $Ls $g_i(Pgain)]
	set Bgain     [lindex $Ls $g_i(Bgain)]
	
	## Section 1: Initial assessment of genomic content
	####################################################################################################################
	if {$REgene eq "" && $NbGenes eq "0"} {
	    # 1B. Does NOT contain protein-coding or any known functionally important elements (- 0.60)
	    set g_rankingScore($AnnotSV_ID) "-0.60"
	    set g_rankingExplanations($AnnotSV_ID) "1B "
	} else {
	    # 1A. Contains protein-coding or other known functionally important elements (+ 0)
	    set g_rankingScore($AnnotSV_ID) "0"
	    set g_rankingExplanations($AnnotSV_ID) "1A "
	}
	
	## Section 2: Overlap with established triplosensitive (TS), haploinsufficient (HI), or benign genes or genomic regions
	##            (Skip to section 3 if the copy-number gain DOES NOT overlap these types of genes/regions)			
	#####################################################################################################################
	if {$Pgain ne ""} {
	    # 2A. Complete overlap; the TS gene or minimal critical region is fully contained within the observed copy-number gain (+ 1)
	    set g_rankingExplanations($AnnotSV_ID,2A) "+ 2A (cf P_gain_source) "
	} elseif {$Bgain ne ""} {
	    # 2D. Smaller than established benign copy-number gain, breakpoint(s) does not interrupt protein-coding genes (- 1)
	    # Need to check with the split lines if breakpoint(s) does not interrupt protein-coding genes (cf location)
	    set g_rankingExplanations($AnnotSV_ID,2D) ""
	}
	
	## Section 3: Evaluation of gene number
	####################################################################################################################
	if {$NbGenes <= 34} {
	    # 3A. 0-34 genes (+ 0)
	    if {$NbGenes < 2} {
		set g_rankingExplanations($AnnotSV_ID,3A) "+ 3A ($NbGenes gene) "
	    } else {
		set g_rankingExplanations($AnnotSV_ID,3A) "+ 3A ($NbGenes genes) "
	    }
	} elseif {$NbGenes <= 49} {
	    # 3B. 35-49 genes (+ 0.45)
	    set g_rankingExplanations($AnnotSV_ID,3B) "+ 3B ($NbGenes genes) "
	} else {
	    # 3C. 50 or more genes (+ 0.90)
	    set g_rankingExplanations($AnnotSV_ID,3C) "+ 3C ($NbGenes genes) "
	} 
	
	##  Section 4: Detailed evaluation of genomic content using cases from published literature, public databases, and/or internal lab data
	##             (Note: If there have been no reports associating either the copy-number gain or any of the genes therein with human phenotypes caused by triplosensitivity, skip to section 5)  	
	####################################################################################################################
	
	
	## Section 5: Evaluation of inheritance pattern/family history for patient being studied			
	####################################################################################################################
	if {$REgene ne ""} {
	    set L_infos [split $REgene ";"]
	    set bestEx "0"
	    foreach infos $L_infos {
		if {[regexp "(.+?) \\\(.*?EX=(.+)\\\)$" $infos match g ex]} {
		    regsub "," $ex "." ex
		    if {$ex > $bestEx} {
			set bestEx $ex
			set bestG $g
		    }
		}
	    }
	    if {$bestEx >= 0.7} {
		# 5H. Inheritance information is unavailable or uninformative. The patient phenotype is highly specific and consistent with what has been described in similar cases (+ 0.15)
		lappend g_rankingExplanations($AnnotSV_ID,5H) "RE:$bestG"
	    } elseif {$bestEx >= 0.5} {
		# 5G. Inheritance information is unavailable or uninformative. The patient phenotype is nonspecific, but is consistent with what has been described in similar cases (+ 0.10)
		lappend g_rankingExplanations($AnnotSV_ID,5G) "RE:$bestG"
	    }
	}
	    
    } else {
	# Split lines
	
	## Section 2: Overlap with established triplosensitive (TS), haploinsufficient (HI), or benign genes or genomic regions
	##            (Skip to section 3 if the copy-number gain DOES NOT overlap these types of genes/regions)			
	#####################################################################################################################
	set HI          [lindex $Ls $g_i(HI)]
	set morbidGenes [lindex $Ls $g_i(morbidGenes)]
	set gene        [lindex $Ls $g_i(gene)]
	set location    [lindex $Ls $g_i(location)]
	set exomiser    [lindex $Ls $g_i(exomiser)]
	regsub "," $exomiser "." exomiser

	if {[info exists g_rankingExplanations($AnnotSV_ID,2D)] && ![regexp "exon|intron" $location]} {
	    # 2D. Smaller than established benign copy-number gain, breakpoint(s) does not interrupt protein-coding genes (- 1)
	    set g_rankingExplanations($AnnotSV_ID,2D) "+ 2D (cf B_gain_source) "
	}

	if {$HI eq "3" || $morbidGenes eq "yes"} {
	    set frameshift [lindex $Ls $g_i(frameshift)]
	    
	    if {![regexp "^txStart" $location] && ![regexp "txEnd$" $location]} {
		# Both breakpoints are within the same HI gene 
		if {$frameshift eq "yes"} {
		    # 2I. Both breakpoints are within the same HI gene (gene-level sequence variant, possibly resulting in loss of function [LOF]).
		    #     and frameshift (+ 0.9)
		    lappend g_rankingExplanations($AnnotSV_ID,2I) "$gene"
		}
	    } elseif {![regexp "^txStart" $location] || ![regexp "txEnd$" $location]} {
		# One breakpoint is within an established HI gene
		if {$exomiser > 0.7} {
		    # 2K. One breakpoint is within an established HI gene, patient’s phenotype is highly specific and consistent with what is expected for LOF of that gene (+ 0.45)
		    lappend g_rankingExplanations($AnnotSV_ID,2K) "$gene"
		}
	    }
	}

	## Section 5: Evaluation of inheritance pattern/family history for patient being studied			
	####################################################################################################################
	if {$exomiser >= 0.7} {
	    # 5H. Inheritance information is unavailable or uninformative. The patient phenotype is highly specific and consistent with what has been described in similar cases (+ 0.15)
	    lappend g_rankingExplanations($AnnotSV_ID,5H) "$gene"
	} elseif {$exomiser >= 0.5} {
	    # 5G. Inheritance information is unavailable or uninformative. The patient phenotype is nonspecific, but is consistent with what has been described in similar cases (+ 0.10)
	    lappend g_rankingExplanations($AnnotSV_ID,5G) "$gene"
	}
    }
   
    return
}


proc achieveSVrankingGain {AnnotSV_ID} {

    global g_rankingScore
    global g_rankingExplanations

    if {![info exists g_rankingScore($AnnotSV_ID)]} {set g_rankingScore($AnnotSV_ID) ""; return}
    
    ## Section 2: Overlap with established triplosensitive (TS), haploinsufficient (HI), or benign genes or genomic regions
    ##            (Skip to section 3 if the copy-number gain DOES NOT overlap these types of genes/regions)			
    #####################################################################################################################
    # Add the higher score of the section 2
    if {[info exists g_rankingExplanations($AnnotSV_ID,2A)]} {
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+1}]
	append g_rankingExplanations($AnnotSV_ID) "+ 2A (cf P_loss_source) "
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2D)] && $g_rankingExplanations($AnnotSV_ID,2D) eq "+ 2D (cf B_gain_source) "} {
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)-1}]
	append g_rankingExplanations($AnnotSV_ID) "+ 2D (cf B_gain_source) "
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2I)]} {
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.9}]
	append g_rankingExplanations($AnnotSV_ID) "+ 2I ([join $g_rankingExplanations($AnnotSV_ID,2I) "/"]) "
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2K)]} {
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.45}]
	append g_rankingExplanations($AnnotSV_ID) "+ 2K ([join $g_rankingExplanations($AnnotSV_ID,2K) "/"]) "
    }
	      
    ## Section 3: Evaluation of gene number
    ####################################################################################################################
    # Add the higher score of the section 3
    if {[info exists g_rankingExplanations($AnnotSV_ID,3C)]} {
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.90}]
	append g_rankingExplanations($AnnotSV_ID) "$g_rankingExplanations($AnnotSV_ID,3C)"
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,3B)]} {
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.45}]
	append g_rankingExplanations($AnnotSV_ID) "$g_rankingExplanations($AnnotSV_ID,3B)"
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,3A)]} {
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.00}]
	append g_rankingExplanations($AnnotSV_ID) "$g_rankingExplanations($AnnotSV_ID,3A)"
    } 
    
    ## Section 4: Detailed evaluation of genomic content using cases from published literature, public databases, and/or internal lab data
    ##            (Note: If there have been no reports associating either the copy-number gain or any of the genes therein with human phenotypes caused by triplosensitivity, skip to section 5)		 
    ####################################################################################################################
    # Add the higher score of the section 4
    
    
    ## Section 5: Evaluation of inheritance patterns/family history for patient being studied			
    ####################################################################################################################
    # Add the higher score of the section 5
    if {![info exists g_rankingExplanations($AnnotSV_ID,2K)]} {
	if {[info exists g_rankingExplanations($AnnotSV_ID,5H)]} {
	    set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.15}]
	    append g_rankingExplanations($AnnotSV_ID) "+ 5H ([join $g_rankingExplanations($AnnotSV_ID,5H) "/"]) "
	} elseif {[info exists g_rankingExplanations($AnnotSV_ID,5G)]} {
	    set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.10}]
	    append g_rankingExplanations($AnnotSV_ID) "+ 5G ([join $g_rankingExplanations($AnnotSV_ID,5G) "/"]) "
	}
    }
    
    return
}




## Rank the pathogenicity of the "INS" SV.
## Creation of 2 global variables:
## - g_rankingScore($AnnotSV_ID)
## - g_rankingExplanations($AnnotSV_ID) 
proc SVrankingINS {L_annotations} {

    global g_AnnotSV
    global g_i
    global g_rankingExplanations
    global g_rankingScore

    set Ls [split $L_annotations "\t"]
    set AnnotSV_ID "[lindex $Ls 0]" 

    set g_rankingScore($AnnotSV_ID) ""
    set g_rankingExplanations($AnnotSV_ID) ""
    
    return
}


## Rank the pathogenicity of the "INV" SV.
## Creation of 2 global variables:
## - g_rankingScore($AnnotSV_ID)
## - g_rankingExplanations($AnnotSV_ID) 
proc SVrankingINV {L_annotations} {

    global g_AnnotSV
    global g_i
    global g_rankingExplanations
    global g_rankingScore

    set Ls [split $L_annotations "\t"]
    set AnnotSV_ID "[lindex $Ls 0]" 

    set g_rankingScore($AnnotSV_ID) ""
    set g_rankingExplanations($AnnotSV_ID) ""
   
    return
}
