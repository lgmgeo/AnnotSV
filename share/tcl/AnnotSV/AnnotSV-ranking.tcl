############################################################################################################
# AnnotSV 3.1.3                                                                                            #
#                                                                                                          #
# AnnotSV: An integrated tool for Structural Variations annotation and ranking                             #
#                                                                                                          #
# Copyright (C) 2017-2022 Veronique Geoffroy (veronique.geoffroy@inserm.fr)                                #
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

    # The user is using a VCF SV input file: the svtTSVcol ("SV_type") is necessarily the 6th in the corresponding output annotated file.tsv
    if {[regexp ".vcf$|.vcf.gz$" $g_AnnotSV(SVinputFile)]} {
	set g_AnnotSV(svtTSVcol) 5
	set g_AnnotSV(svtBEDcol) 3
    }

    # "SV_type" is compulsory for the ranking
    if {$g_AnnotSV(svtTSVcol) eq -1} {return}

    # Check if we have all the needed information for ranking
    # unset g_i <=> no SV ranking  
    set g_i(gene)     [lsearch -regexp $Ls "Gene_name"];      if {$g_i(gene) == -1} {unset g_i; return}  
    set g_i(NbGenes)  [lsearch -regexp $Ls "Gene_count"]; if {$g_i(NbGenes) == -1} {unset g_i; return}  
    set g_i(full)     [lsearch -regexp $Ls "Annotation_mode"];   if {$g_i(full) == -1} {unset g_i; return}        

    set g_i(REgene)   [lsearch -regexp $Ls "RE_gene"]; if {$g_i(REgene) == -1} {unset g_i; return}
    
    set g_i(Ploss)      [lsearch -regexp $Ls "P_loss_coord"];  if {$g_i(Ploss) == -1} {unset g_i; return}  
    set g_i(Pgain)      [lsearch -regexp $Ls "P_gain_coord"];  if {$g_i(Pgain) == -1} {unset g_i; return}  
    set g_i(Bloss)      [lsearch -regexp $Ls "B_loss_coord"];  if {$g_i(Bloss) == -1} {unset g_i; return}  
    set g_i(Bgain)      [lsearch -regexp $Ls "B_gain_coord"];  if {$g_i(Bgain) == -1} {unset g_i; return}  
    set g_i(Psnvindel)  [lsearch -regexp $Ls "P_snvindel_nb"]; if {$g_i(Psnvindel) == -1} {unset g_i; return}
    
    set g_i(poPloss)    [lsearch -regexp $Ls "po_P_loss_coord"];  if {$g_i(poPloss) == -1} {unset g_i; return}  
    set g_i(poBlossAllG)   [lsearch -regexp $Ls "po_B_loss_allG_coord"];  if {$g_i(poBlossAllG) == -1} {unset g_i; return}
    set g_i(poBlossSomeG)  [lsearch -regexp $Ls "po_B_loss_someG_coord"]; if {$g_i(poBlossSomeG) == -1} {unset g_i; return}

    set g_i(HI) [lsearch -regexp $Ls "HI"];     if {$g_i(HI) == -1} {unset g_i; return}  
    set g_i(TS) [lsearch -regexp $Ls "TS"];     if {$g_i(TS) == -1} {unset g_i; return}  

    set g_i(pLI)       [lsearch -regexp $Ls "GnomAD_pLI"];    if {$g_i(pLI) == -1} {unset g_i; return}  
    set g_i(loeuf)     [lsearch -regexp $Ls "LOEUF_bin"];     if {$g_i(loeuf) == -1} {unset g_i; return}  
    set g_i(HIpercent) [lsearch -regexp $Ls "DDD_HI_percent"]; if {$g_i(HIpercent) == -1} {unset g_i; return}  
   
    set g_i(exomiser)  [lsearch -regexp $Ls "Exomiser_gene_pheno_score"];  

    set g_i(morbid)    [lsearch -regexp $Ls "OMIM_morbid"]; if {$g_i(morbid) == -1} {unset g_i; return}   

    set g_i(location)    [lsearch -regexp $Ls "Location"];               if {$g_i(location) == -1} {unset g_i; return}  
    set g_i(location2)   [lsearch -regexp $Ls "Location2"];              if {$g_i(location2) == -1} {unset g_i; return} 
    set g_i(CDSlength)   [lsearch -regexp $Ls "Overlapped_CDS_length"];  if {$g_i(CDSlength) == -1} {unset g_i; return}
    set g_i(CDSpercent)  [lsearch -regexp $Ls "Overlapped_CDS_percent"]; if {$g_i(CDSpercent) == -1} {unset g_i; return} 
    set g_i(frameshift)  [lsearch -regexp $Ls "Frameshift"];             if {$g_i(frameshift) == -1} {unset g_i; return} 
    set g_i(NbExons)     [lsearch -regexp $Ls "Exon_count"];        if {$g_i(NbExons) == -1} {unset g_i; return} 

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
    set AnnotationMode [lindex $Ls $g_i(full)]
    if {$AnnotationMode eq "full"} {

	# In case of SV redundancy in the input file
	if {[info exists g_i(fullDone,$AnnotSV_ID)]} {
	    set g_i(splitDone,$AnnotSV_ID) 1
	    return
	}
	set g_i(fullDone,$AnnotSV_ID) 1
	
	set NbGenes   [lindex $Ls $g_i(NbGenes)]
	set REgene    [lindex $Ls $g_i(REgene)]
	set Ploss     [lindex $Ls $g_i(Ploss)]
	set poPloss   [lindex $Ls $g_i(poPloss)]
	set Bloss     [lindex $Ls $g_i(Bloss)]
	set poBlossSomeG [lindex $Ls $g_i(poBlossSomeG)]
	set poBlossAllG  [lindex $Ls $g_i(poBlossAllG)]
	
	## Section 1: Initial assessment of genomic content
	####################################################################################################################
	if {$REgene eq "" && $NbGenes eq "0"} {
	    # 1B. Does NOT contain protein-coding or any known functionally important elements (-0.60)
	    set g_rankingScore($AnnotSV_ID) "-0.60"
	    set g_rankingExplanations($AnnotSV_ID) "1B (-0.60);"
	} else {
	    # 1A. Contains protein-coding or other known functionally important elements (+0.00)
	    set g_rankingScore($AnnotSV_ID) "0"
	    set g_rankingExplanations($AnnotSV_ID) "1A (+0.00);"
	}
	
	## Section 2: Overlap with established/predicted haploinsufficiency (HI) or established benign genes/genomic regions
	##            (Skip to section 3 if your copy-number loss DOES NOT overlap these types of genes/regions)
	#####################################################################################################################
	if {$Ploss ne ""} {
	    # 2A. Complete overlap of an established HI gene/genomic region (+1.00)
	    set g_rankingExplanations($AnnotSV_ID,2A) "2A (cf P_loss_source, +1.00);"
	} else {
	    if {$poPloss ne ""} {
		# 2B. Partial overlap of a known pathogenic Loss SV (+0.00)
		set g_rankingExplanations($AnnotSV_ID,2B) "2B (cf po_P_loss_source, +0.00);"
	    }
	    if {$Bloss ne ""} {
		# 2F. Completely contained within an established benign CNV region (-1.00)
		set g_rankingExplanations($AnnotSV_ID,2F) "2F (cf B_loss_source, -1.00);"
	    } elseif {$poBlossSomeG ne ""} {		
		# 2G. Overlaps an established benign CNV, but includes additional genomic material (+0.00)
		set g_rankingExplanations($AnnotSV_ID,2G) "2G (+0.00);"
	    }
	}
	
	## Section 3: Evaluation of gene number
	####################################################################################################################
	if {$NbGenes <= 24} {
	    # 3A. 0-24 genes (+0.00)
	    if {$NbGenes < 2} {
		set g_rankingExplanations($AnnotSV_ID,3A) "3A ($NbGenes gene, +0.00);"
	    } else {
		set g_rankingExplanations($AnnotSV_ID,3A) "3A ($NbGenes genes, +0.00);"
	    }
	} elseif {$NbGenes <= 35} {
	    # 3B. 25-34 genes (+0.45)
	    set g_rankingExplanations($AnnotSV_ID,3B) "3B ($NbGenes genes, +0.45);"
	} else {
	    # 3C. 35+ genes (+0.90)
	    set g_rankingExplanations($AnnotSV_ID,3C) "3C ($NbGenes genes, +0.90);"
	} 
	
	## Section 4: Detailed evaluation of genomic content using cases from published literature, public databases, and/or internal lab data
	##            (Skip to section 5 if either your CNV overlapped with an established HI gene/region in section 2,
	##             OR there have been no reports associating either the CNV or any genes within the CNV with human phenotypes
	##             caused by loss of function [LOF] or copy-number loss)
	####################################################################################################################

	if {$poBlossAllG ne ""} {
	    # 4O. Overlap with common population variation (completely contained within a common population CNV
            #     OR contains no additional genomic material).(-1.00)
	    set g_rankingExplanations($AnnotSV_ID,40) "40 (-1.00);"    
	}
	
	## Section 5: Evaluation of inheritance pattern/family history for patient being studied			
	####################################################################################################################
	if {$REgene ne ""} {
	    set L_infos [split $REgene "\\\);"]
	    set bestEx "0"
	    foreach infos $L_infos {
		if {[regexp "(.+?) \\\(.*?EX=(.+)$" $infos match g ex]} {
		    regsub "," $ex "." ex
		    if {$ex > $bestEx} {
			set bestEx $ex
			set bestG $g
		    }
		}
	    }
	    if {$bestEx >= 0.7} {
		# 5H. Inheritance information is unavailable or uninformative. The patient phenotype is highly specific and consistent with what has been described in similar cases (+0.30)
		lappend g_rankingExplanations($AnnotSV_ID,5H) "RE:$bestG"
	    } elseif {$bestEx >= 0.5} {
		# 5G. Inheritance information is unavailable or uninformative. The patient phenotype is nonspecific, but is consistent with what has been described in similar cases (+0.10)
		lappend g_rankingExplanations($AnnotSV_ID,5G) "RE:$bestG"
	    }
	}
	    
    } else {
	## Split lines

	# In case of SV redundancy in the input file
	if {[info exists g_i(splitDone,$AnnotSV_ID)]} {return}
	
	## Section 2: Overlap with established/predicted haploinsufficiency (HI) or established benign genes/genomic regions
	##            (Skip to section 3 if your copy-number loss DOES NOT overlap these types of genes/regions)
	#####################################################################################################################
	set HI     [lindex $Ls $g_i(HI)]
	set morbid [lindex $Ls $g_i(morbid)]
	set gene   [lindex $Ls $g_i(gene)]
	set exomiser    [lindex $Ls $g_i(exomiser)]
	regsub "," $exomiser "." exomiser

	if {$HI eq "3" || $morbid eq "yes"} {
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
		    # 2C-1. …and coding sequence is involved (+0.90)
		    lappend g_rankingExplanations($AnnotSV_ID,2C-1) "$gene"
		} else {
		    # 2C-2. …and only the 5’ UTR is involved (+0.00)
		    lappend g_rankingExplanations($AnnotSV_ID,2C-2) "$gene"		    
		}
	    } elseif {[regexp "3'UTR$" $location2] && ![regexp "^5'UTR" $location2]} {
		# 2D. Partial overlap with the 3’ end of an established HI gene (5’ end of the gene not involved)…
		if {$CDSlength eq "0"} {
		    # 2D-1. …and only the 3’ untranslated region is involved (+0.00)
		    lappend g_rankingExplanations($AnnotSV_ID,2D-1) "$gene"		    
		} elseif {$location eq "exon${NbExons}-txEnd" && $Psnvindel ne ""} {
		    # 2D-2. …and only the last exon is involved. Other established pathogenic snv/indel have been reported in the SV (+0.90)
		    # Vero to improve: Other established pathogenic snv/indel have been reported in this exon
		    lappend g_rankingExplanations($AnnotSV_ID,2D-2) "$gene"
		} elseif {$location eq "exon${NbExons}-txEnd"} {
		    # 2D-3. …and only the last exon is involved. No other established pathogenic snv/indel have been reported in the SV (+0.30)
		    # Vero to improve: No other established pathogenic snv/indel have been reported in this exon
		    lappend g_rankingExplanations($AnnotSV_ID,2D-3) "$gene"
		} elseif {$location ne "exon${NbExons}-txEnd" && $location ne "intron[expr {${NbExons}-1}]-txEnd"} {	    
		    # 2D-4. …and it includes other exons in addition to the last exon. Nonsense-mediated decay is expected to occur (+0.90)
		    lappend g_rankingExplanations($AnnotSV_ID,2D-4) "$gene"		    
		}
	    } elseif {![regexp "txStart|txEnd" $location]} {
		# 2E. Both breakpoints are within the same gene (intragenic CNV; gene-level sequence variant) AND...
		set frameshift [lindex $Ls $g_i(frameshift)]
		# At least 1 exon overlapped ?? <=> (exon'i' in location) OR (location eq intron'i'-intron'j') BUT NOT (location eq intron'i'-intron'i')
		set exonOverlapped 0
		if {[regexp "exon" $location]} {
		    set exonOverlapped 1
		} elseif {[regsub -all "intron" $location "" tutu]} {
		    if {[string index $tutu 0] ne [string index $tutu 2]} {
			set exonOverlapped 1
		    }
		}
		if {$frameshift eq "yes"} {
		    # 2E-1. ...frameshift (+0.90)
		    lappend g_rankingExplanations($AnnotSV_ID,2E-1) "$gene"
		} elseif {$exonOverlapped && $Psnvindel ne "" && $CDSpercent >= 10} {
		    # 2E-2. ...exon(s) overlapped AND pathogenic snv/indel overlapped AND SV removes > 10% of protein (+0.45)
		    lappend g_rankingExplanations($AnnotSV_ID,2E-2) "$gene"		    
		} elseif {$exonOverlapped && $Psnvindel ne "" && $CDSpercent < 10} {
		    # 2E-3. ...exon(s) overlapped AND pathogenic snv/indel overlapped AND SV removes < 10% of protein (+0.30)
		    lappend g_rankingExplanations($AnnotSV_ID,2E-3) "$gene"		    
		} elseif {$exonOverlapped && $Psnvindel eq "" && $CDSpercent > 10} {
		    # 2E-4. ...>=1 exon deleted AND no established pathogenic snv/indel have been reported in the observed CNV AND variant removes > 10% of protein (+0.20)
		    lappend g_rankingExplanations($AnnotSV_ID,2E-4) "$gene"
		}
	    }
	} else {
	    set pLI        [lindex $Ls $g_i(pLI)]
	    regsub -all "," $pLI "." pLI
	    set loeuf      [lindex $Ls $g_i(loeuf)]
	    regsub -all "," $loeuf "." loeuf
	    set HIpercent  [lindex $Ls $g_i(HIpercent)]
	    regsub -all "," $HIpercent "." HIpercent
	    set i 0
	    if {$pLI >= 0.9} {incr i}
	    if {$HIpercent <= 10} {incr i}
	    if {$i >= 2} {
		# 2H. Two or more HI predictors suggest that AT LEAST ONE gene in the interval is HI (+0.15)
		lappend g_rankingExplanations($AnnotSV_ID,2H) "$gene"
	    }
#	    if {$loeuf != ""} {incr i} 
#	    if {$i >= 2} {
#		# 2H. Two or more HI predictors suggest that AT LEAST ONE gene in the interval is HI (+0.15)
#		lappend g_rankingExplanations($AnnotSV_ID,2H) "$gene"
#		set scoreLoeuf(0) "0.075"
#		set scoreLoeuf(1) "0.0675"
#		set scoreLoeuf(2) "0.06"
#		set scoreLoeuf(3) "0.0525"
#		set scoreLoeuf(4) "0.0450"
#		set scoreLoeuf(5) "0.0375"
#		set scoreLoeuf(6) "0.03"
#		set scoreLoeuf(7) "0.0225"
#		set scoreLoeuf(8) "0.0150"
#		set scoreLoeuf(9) "0.0075"
#		set score2H [expr {0.075+$scoreLoeuf($loeuf)}]
#		if {![info exists g_rankingScore($AnnotSV_ID,maxLoeuf)]} {
#		    set g_rankingScore($AnnotSV_ID,maxLoeuf) $score2H
#		} elseif {$score2H > $g_rankingScore($AnnotSV_ID,maxLoeuf)} {
#		    set g_rankingScore($AnnotSV_ID,maxLoeuf) $score2H
#		}
#	    } 
	}

	## Section 5: Evaluation of inheritance pattern/family history for patient being studied			
	####################################################################################################################
	if {$exomiser >= 0.7} {
	    # 5H. Inheritance information is unavailable or uninformative. The patient phenotype is highly specific and consistent with what has been described in similar cases (+0.30)
	    lappend g_rankingExplanations($AnnotSV_ID,5H) "$gene"
	} elseif {$exomiser >= 0.5} {
	    # 5G. Inheritance information is unavailable or uninformative. The patient phenotype is nonspecific, but is consistent with what has been described in similar cases (+0.10)
	    lappend g_rankingExplanations($AnnotSV_ID,5G) "$gene"
	}
    }
   
    return
}

# Add the higher score of each section
# (WARNING: Evaluation of the different "if" statements should be done from the highest score to the lowest score)
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

    # Add the higher score of the section 2:
    # - Evaluate scores from the higher (positive value) to the lower (negative value); stop the evaluation at the first match (only the higher match is reported)
    #   (the "0" scores are evaluated in a second time)
    if {[info exists g_rankingExplanations($AnnotSV_ID,2A)]} {
	# 2A
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+1.00}]
	append g_rankingExplanations($AnnotSV_ID) "$g_rankingExplanations($AnnotSV_ID,2A)"
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2C-1)]} {
	# 2C-1
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.90}]
	append g_rankingExplanations($AnnotSV_ID) "2C-1 ([join [lsort -unique $g_rankingExplanations($AnnotSV_ID,2C-1)] "/"], +0.90);"
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2D-2)]} {
	# 2D-2
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.90}]
	append g_rankingExplanations($AnnotSV_ID) "2D-2 ([join [lsort -unique $g_rankingExplanations($AnnotSV_ID,2D-2)] "/"], +0.90);"
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2D-4)]} {
	# 2D-4
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.90}]
	append g_rankingExplanations($AnnotSV_ID) "2D-4 ([join [lsort -unique $g_rankingExplanations($AnnotSV_ID,2D-4)] "/"], +0.90);"
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2E-1)]} {
	# 2E-1
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.90}]
	append g_rankingExplanations($AnnotSV_ID) "2E-1 ([join [lsort -unique $g_rankingExplanations($AnnotSV_ID,2E-1)] "/"], +0.90);"
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2E-2)]} {
	# 2E-2
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.45}]
	append g_rankingExplanations($AnnotSV_ID) "2E-2 ([join [lsort -unique $g_rankingExplanations($AnnotSV_ID,2E-2)] "/"], +0.45);"
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2D-3)]} {
	# 2D-3
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.45}]
	append g_rankingExplanations($AnnotSV_ID) "2D-3 ([join [lsort -unique $g_rankingExplanations($AnnotSV_ID,2D-3)] "/"], +0.45);"
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2E-3)]} {
	# 2E-3
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.30}]
	append g_rankingExplanations($AnnotSV_ID) "2E-3 ([join [lsort -unique $g_rankingExplanations($AnnotSV_ID,2E-3)] "/"], +0.30);"
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2E-4)]} {
	# 2E-4
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.20}]
	append g_rankingExplanations($AnnotSV_ID) "2E-4 ([join [lsort -unique $g_rankingExplanations($AnnotSV_ID,2E-4)] "/"]), +0.20;"
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2H)]} {
	# 2H
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.15}]
	#set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+$g_rankingScore($AnnotSV_ID,maxLoeuf)}]
	append g_rankingExplanations($AnnotSV_ID) "2H ([join [lsort -unique $g_rankingExplanations($AnnotSV_ID,2H)] "/"], +0.15);"
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2F)]} {
	# 2F
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)-1.00}]
	append g_rankingExplanations($AnnotSV_ID) "$g_rankingExplanations($AnnotSV_ID,2F)"
    }

    # - Evaluate the "0" scores; all the match are reported.
    if {[info exists g_rankingExplanations($AnnotSV_ID,2B)] && ![regexp "2A|2C|2D|2E" $g_rankingExplanations($AnnotSV_ID)]} {
	# 2B
	append g_rankingExplanations($AnnotSV_ID) "$g_rankingExplanations($AnnotSV_ID,2B)"
    }
    if {[info exists g_rankingExplanations($AnnotSV_ID,2C-2)]} {
	# 2C-2
	append g_rankingExplanations($AnnotSV_ID) "2C-2 ([join [lsort -unique $g_rankingExplanations($AnnotSV_ID,2C-2)] "/"],+0.00);"
    }
    if {[info exists g_rankingExplanations($AnnotSV_ID,2D-1)]} {
	# 2D-1
	append g_rankingExplanations($AnnotSV_ID) "2D-1 ([join [lsort -unique $g_rankingExplanations($AnnotSV_ID,2D-1)] "/"],+0.00);"
    }
    if {[info exists g_rankingExplanations($AnnotSV_ID,2G)]} {
	# 2G
	append g_rankingExplanations($AnnotSV_ID) "$g_rankingExplanations($AnnotSV_ID,2G)"
    } 
    
    ## Section 3: Evaluation of gene number
    ####################################################################################################################

    # Add the higher score of the section 3
    if {[info exists g_rankingExplanations($AnnotSV_ID,3C)]} {
	# 3C
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.90}]
	append g_rankingExplanations($AnnotSV_ID) "$g_rankingExplanations($AnnotSV_ID,3C)"
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,3B)]} {
	# 3B
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.45}]
	append g_rankingExplanations($AnnotSV_ID) "$g_rankingExplanations($AnnotSV_ID,3B)"
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,3A)]} {
	# 3A
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.00}]
	append g_rankingExplanations($AnnotSV_ID) "$g_rankingExplanations($AnnotSV_ID,3A)"
     } 
    
    ## Section 4: Detailed evaluation of genomic content using cases from published literature, public databases, and/or internal lab data
    ##            (Skip to section 5 if either your CNV overlapped with an established HI gene/region in section 2,
    ##             OR there have been no reports associating either the CNV or any genes within the CNV with human phenotypes
    ##             caused by loss of function [LOF] or copy-number loss)
    ####################################################################################################################

    # Add the higher score of the section 4

    if {[info exists g_rankingExplanations($AnnotSV_ID,40)]} {
        # 4O
        set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)-1.00}]
        append g_rankingExplanations($AnnotSV_ID) "$g_rankingExplanations($AnnotSV_ID,4O)"
    } 
 

    ## Section 5: Evaluation of inheritance pattern/family history for patient being studied			
    ####################################################################################################################

    # Add the higher score of the section 5
    if {[info exists g_rankingExplanations($AnnotSV_ID,5H)]} {
	# 5H
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.30}]
	append g_rankingExplanations($AnnotSV_ID) "5H ([join [lsort -unique $g_rankingExplanations($AnnotSV_ID,5H)] "/"], +0.30)"
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,5G)]} {
	# 5G
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.10}]
	append g_rankingExplanations($AnnotSV_ID) "5G ([join [lsort -unique $g_rankingExplanations($AnnotSV_ID,5G)] "/"], +0.10)"
    } else {
	# 5F
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.00}]
	append g_rankingExplanations($AnnotSV_ID) "5F (+0.00)"
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

    set AnnotationMode [lindex $Ls $g_i(full)]
    
    if {$AnnotationMode eq "full"} {
	
	set NbGenes   [lindex $Ls $g_i(NbGenes)]
	set REgene    [lindex $Ls $g_i(REgene)]
	set Pgain     [lindex $Ls $g_i(Pgain)]
	set Bgain     [lindex $Ls $g_i(Bgain)]
	
	## Section 1: Initial assessment of genomic content
	####################################################################################################################
	if {$REgene eq "" && $NbGenes eq "0"} {
	    # 1B. Does NOT contain protein-coding or any known functionally important elements (-0.60)
	    set g_rankingScore($AnnotSV_ID) "-0.60"
	    set g_rankingExplanations($AnnotSV_ID) "1B (-0.60);"
	} else {
	    # 1A. Contains protein-coding or other known functionally important elements (+0.00)
	    set g_rankingScore($AnnotSV_ID) "0"
	    set g_rankingExplanations($AnnotSV_ID) "1A (+0.00);"
	}
	
	## Section 2: Overlap with established triplosensitive (TS), haploinsufficient (HI), or benign genes or genomic regions
	##            (Skip to section 3 if the copy-number gain DOES NOT overlap these types of genes/regions)			
	#####################################################################################################################
	if {$Pgain ne ""} {
	    # 2A. Complete overlap; the TS gene or minimal critical region is fully contained within the observed copy-number gain (+1.00)
	    set g_rankingExplanations($AnnotSV_ID,2A) "2A (cf P_gain_source, +1.00);"
	} elseif {$Bgain ne ""} {
	    # 2DE. Smaller than established benign copy-number gain
	    set g_rankingExplanations($AnnotSV_ID,2DE) ""
	}
	
	## Section 3: Evaluation of gene number
	####################################################################################################################
	if {$NbGenes <= 34} {
	    # 3A. 0-34 genes (+0.00)
	    if {$NbGenes < 2} {
		set g_rankingExplanations($AnnotSV_ID,3A) "3A ($NbGenes gene, +0.00);"
	    } else {
		set g_rankingExplanations($AnnotSV_ID,3A) "3A ($NbGenes genes, +0.00);"
	    }
	} elseif {$NbGenes <= 49} {
	    # 3B. 35-49 genes (+0.45)
	    set g_rankingExplanations($AnnotSV_ID,3B) "3B ($NbGenes genes, +0.45);"
	} else {
	    # 3C. 50 or more genes (+0.90)
	    set g_rankingExplanations($AnnotSV_ID,3C) "3C ($NbGenes genes, +0.90);"
	} 
	
	##  Section 4: Detailed evaluation of genomic content using cases from published literature, public databases, and/or internal lab data
	##             (Note: If there have been no reports associating either the copy-number gain or any of the genes therein with human phenotypes caused by triplosensitivity, skip to section 5)  	
	####################################################################################################################
	
	
	## Section 5: Evaluation of inheritance pattern/family history for patient being studied			
	####################################################################################################################
	if {$REgene ne ""} {
	    set L_infos [split $REgene "\\\);"]
	    set bestEx "0"
	    foreach infos $L_infos {
		if {[regexp "(.+?) \\\(.*?EX=(.+)$" $infos match g ex]} {
		    regsub "," $ex "." ex
		    if {$ex > $bestEx} {
			set bestEx $ex
			set bestG $g
		    }
		}
	    }
	    if {$bestEx >= 0.7} {
		# 5H. Inheritance information is unavailable or uninformative. The patient phenotype is highly specific and consistent with what has been described in similar cases (+0.15)
		lappend g_rankingExplanations($AnnotSV_ID,5H) "RE:$bestG"
	    } elseif {$bestEx >= 0.5} {
		# 5G. Inheritance information is unavailable or uninformative. The patient phenotype is nonspecific, but is consistent with what has been described in similar cases (+0.10)
		lappend g_rankingExplanations($AnnotSV_ID,5G) "RE:$bestG"
	    }
	}
	
    } else {
	# Split lines
	
	## Section 2: Overlap with established triplosensitive (TS), haploinsufficient (HI), or benign genes or genomic regions
	##            (Skip to section 3 if the copy-number gain DOES NOT overlap these types of genes/regions)			
	#####################################################################################################################
	set HI          [lindex $Ls $g_i(HI)]
	set morbid      [lindex $Ls $g_i(morbid)]
	set gene        [lindex $Ls $g_i(gene)]
	set location    [lindex $Ls $g_i(location)]
	set exomiser    [lindex $Ls $g_i(exomiser)]
	regsub "," $exomiser "." exomiser

        if {[info exists g_rankingExplanations($AnnotSV_ID,2DE)]} {
	    if {![regexp "exon|intron" $location]} {
		# 2D. Smaller than established benign copy-number gain, breakpoint(s) does not interrupt protein-coding genes (-1.00)
		set g_rankingExplanations($AnnotSV_ID,2D) "2D (cf B_gain_source, -1.00);"
	    } else {
		# 2E. Smaller than established benign copy-number gain, breakpoint(s) potentially interrupts protein-coding gene (+0.00)
		set g_rankingExplanations($AnnotSV_ID,2E) "2E (cf B_gain_source, +0.00);"
	    }
	} 
	
	if {$HI eq "3" || $morbid eq "yes"} {
	    set frameshift [lindex $Ls $g_i(frameshift)]
	    
	    if {![regexp "^txStart" $location] && ![regexp "txEnd$" $location]} {
		# Both breakpoints are within the same HI gene / morbid gene (gene-level sequence variant, possibly resulting in loss of function [LOF]) AND...
		if {$frameshift eq "yes"} {
		    # 2I-1. ...disrupts the reading frame (+0.45)
		    lappend g_rankingExplanations($AnnotSV_ID,2I-1) "$gene"
		} elseif {$exomiser >= 0.7} {
		    # 2I-2. ...patient phenotype is highly specific and consistent with what has been described for this HI gene / morbid gene (Exomiser_gene_pheno_score >= 0.7) (+0.45)
		    lappend g_rankingExplanations($AnnotSV_ID,2I-2) "$gene"
		} else {
		    # 2I-3. ...patient phenotype is either inconsistent with what has been described for this HI gene / morbid gene (Exomiser_gene_pheno_score < 0.7) OR unknown (+0.00)
		    lappend g_rankingExplanations($AnnotSV_ID,2I-3) "$gene"
		}
	    } elseif {![regexp "^txStart" $location] || ![regexp "txEnd$" $location]} {
		# One breakpoint is within an established HI gene
		if {$exomiser >= 0.7} {
		    # 2K. One breakpoint is within an established HI gene / morbid gene, patient’s phenotype is highly specific and consistent with what is expected
		    #     for LOF of that gene (Exomiser_gene_pheno_score >= 0.7) (+0.45)
		    lappend g_rankingExplanations($AnnotSV_ID,2K) "$gene"
		} else {
		    # 2J. One breakpoint is within an established HI gene / morbid gene, patient’s phenotype is either inconsistent with what is expected for LOF of that gene OR unknown (+0.00) 
		    lappend g_rankingExplanations($AnnotSV_ID,2J) "$gene"
		}
	    } elseif {$location eq "txStart-txEnd"} {
		# 2H. HI gene / morbid gene fully contained within observed copy-number gain AND...
		if {$exomiser >= 0.7} {
		    # 2H-1. ...patient’s phenotype is highly specific and consistent with what is expected for LOF of that gene (Exomiser_gene_pheno_score >= 0.7) (+0.45)
		    lappend g_rankingExplanations($AnnotSV_ID,2H-1) "$gene"
		} else {
		    # 2H-2. ...patient's phenotype is nonspecific with what is expected for LOF of that gene (Exomiser_gene_pheno_score < 0.7) (+0.00)
		    lappend g_rankingExplanations($AnnotSV_ID,2H-2) "$gene"
		}
	    }
	} elseif {![regexp "^txStart" $location] || ![regexp "txEnd$" $location]} {
	    # 2L. One or both breakpoints are within gene(s) of no established clinical significance (+0.00)
	    lappend g_rankingExplanations($AnnotSV_ID,2L) "$gene"
	}

	## Section 5: Evaluation of inheritance pattern/family history for patient being studied			
	####################################################################################################################
	if {$exomiser >= 0.7} {
	    # 5H. Inheritance information is unavailable or uninformative. The patient phenotype is highly specific and consistent with what has been described in similar cases (+0.15)
	    lappend g_rankingExplanations($AnnotSV_ID,5H) "$gene"
	} elseif {$exomiser >= 0.5} {
	    # 5G. Inheritance information is unavailable or uninformative. The patient phenotype is nonspecific, but is consistent with what has been described in similar cases (+0.10)
	    lappend g_rankingExplanations($AnnotSV_ID,5G) "$gene"
	}
    }
   
    return
}


# Add the higher score of each section
# (WARNING: Evaluation of the different "if" statements should be done from the highest score to the lowest score)
proc achieveSVrankingGain {AnnotSV_ID} {

    global g_rankingScore
    global g_rankingExplanations

    if {![info exists g_rankingScore($AnnotSV_ID)]} {set g_rankingScore($AnnotSV_ID) ""; return}
    
    ## Section 2: Overlap with established triplosensitive (TS), haploinsufficient (HI), or benign genes or genomic regions
    ##            (Skip to section 3 if the copy-number gain DOES NOT overlap these types of genes/regions)			
    #####################################################################################################################

    # Add the higher score of the section 2
    if {[info exists g_rankingExplanations($AnnotSV_ID,2A)]} {
	# 2A
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+1.00}]
	append g_rankingExplanations($AnnotSV_ID) "$g_rankingExplanations($AnnotSV_ID,2A)"
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2I-1)]} {
	# 2I-1
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.45}]
	append g_rankingExplanations($AnnotSV_ID) "2I-1 ([join [lsort -unique $g_rankingExplanations($AnnotSV_ID,2I-1)] "/"], +0.45);"
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2I-2)]} {
	# 2I-2
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.45}]
	append g_rankingExplanations($AnnotSV_ID) "2I-2 ([join [lsort -unique $g_rankingExplanations($AnnotSV_ID,2I-2)] "/"], +0.45);"
     } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2K)]} {
	# 2K
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.45}]
	append g_rankingExplanations($AnnotSV_ID) "2K ([join [lsort -unique $g_rankingExplanations($AnnotSV_ID,2K)] "/"], +0.45);"
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2H-1)]} {
	# 2H-1
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.45}]
	append g_rankingExplanations($AnnotSV_ID) "2H-1 ([join [lsort -unique $g_rankingExplanations($AnnotSV_ID,2H-1)] "/"], +0.45);"
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2E)]} {
	# 2E
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.00}]
	append g_rankingExplanations($AnnotSV_ID) "$g_rankingExplanations($AnnotSV_ID,2E)"
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2H-2)]} {
	# 2H-2
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.00}]
	append g_rankingExplanations($AnnotSV_ID) "2H-2 ([join [lsort -unique $g_rankingExplanations($AnnotSV_ID,2H-2)] "/"], +0.00);"
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2I-3)]} {
	# 2I-3
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.00}]
	append g_rankingExplanations($AnnotSV_ID) "2I-3 ([join [lsort -unique $g_rankingExplanations($AnnotSV_ID,2I-3)] "/"], +0.00);"
     } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2J)]} {
	# 2J
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.00}]
	append g_rankingExplanations($AnnotSV_ID) "2J ([join [lsort -unique $g_rankingExplanations($AnnotSV_ID,2J)] "/"], +0.00);"
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2L)]} {
	# 2L
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.00}]
	append g_rankingExplanations($AnnotSV_ID) "2L ([join [lsort -unique $g_rankingExplanations($AnnotSV_ID,2L)] "/"], +0.00);"
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,2D)]} {
	# 2D
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)-1.00}]
	append g_rankingExplanations($AnnotSV_ID) "$g_rankingExplanations($AnnotSV_ID,2D)"
    } 
	      
    ## Section 3: Evaluation of gene number
    ####################################################################################################################

    # Add the higher score of the section 3
    if {[info exists g_rankingExplanations($AnnotSV_ID,3C)]} {
	# 3C
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.90}]
	append g_rankingExplanations($AnnotSV_ID) "$g_rankingExplanations($AnnotSV_ID,3C)"
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,3B)]} {
	# 3B
	set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.45}]
	append g_rankingExplanations($AnnotSV_ID) "$g_rankingExplanations($AnnotSV_ID,3B)"
    } elseif {[info exists g_rankingExplanations($AnnotSV_ID,3A)]} {
	# 3A
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
	    # 5H
	    set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.15}]
	    append g_rankingExplanations($AnnotSV_ID) "5H ([join [lsort -unique $g_rankingExplanations($AnnotSV_ID,5H)] "/"], +0.15)"
	} elseif {[info exists g_rankingExplanations($AnnotSV_ID,5G)]} {
	    # 5G
	    set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.10}]
	    append g_rankingExplanations($AnnotSV_ID) "5G ([join [lsort -unique $g_rankingExplanations($AnnotSV_ID,5G)] "/"], +0.10)"
	} else {
	    # 5F
	    set g_rankingScore($AnnotSV_ID) [expr {$g_rankingScore($AnnotSV_ID)+0.00}]
	    append g_rankingExplanations($AnnotSV_ID) "5F (+0.00)"
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
