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


# No link is done without a check of the genomic location
#########################################################
# - AGS3 example:
#		grep AGS3 $ANNOTSV/share/AnnotSV/Annotations_Human/Gene-based/OMIM/20240806_morbid.tsv
#			AGS3    yes
#			=> AGS3 is morbid
#		grep AGS3 $ANNOTSV/share/AnnotSV/Annotations_Human/Genes/GRCh38/genes.RefSeq.sorted.bed
#			=> No genomic coordinates
#		grep AGS3 $ANNOTSV/share/AnnotSV/Annotations_Human/Gene-based/NCBIgeneID/geneSymbol_NCBIgeneID.tsv
#			AGS3    26086
#			=> NCBIgeneID(AGS3) = 26086
#		grep 26086 $ANNOTSV/share/AnnotSV/Annotations_Human/Gene-based/NCBIgeneID/geneSymbol_NCBIgeneID.tsv
#			GPSM1   26086
#			AGS3    26086
#			DKFZP727I051    26086
#			=> Without AGS3 coordinates, we can't link this gene to other genes (GPSM1, DKFZP727I051)
proc memorizeGeneNameAlias {} {

    global g_AnnotSV
	global g_alias
	global g_L_coord

    # Here, we memorize the RefSeq/ENSEMBL GRCh38 gene coordinates
    # => Will permit to check the different alias of a unique gene (only for overlapping genomic coordinates. cf issues 156 + 132)
    set geneCoordFile "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Genes/GRCh38/genes.$g_AnnotSV(tx).sorted.bed"
    foreach L [LinesFromFile $geneCoordFile] {
        set chrom [lindex $L 0]
        set start [lindex $L 1]
        set end [lindex $L 2]
        set gene [lindex $L 4]
		# We keep all the coordinates for each gene (different transcripts; genes present on the chrX and chrY (e.g. WASH6P); genes present at different locations (e.g. SNORD81))
		#
		# Example: WASH6P 
		# grep WASH6P $ANNOTSV/share/AnnotSV/Annotations_Human/Genes/GRCh38/genes.ENSEMBL.sorted.bed
		# X    156020825    156022415   +   WASH6P  ENST00000476066 156022415       156022415       156020825,156021505,156021693,156021998,156022322,      156021148,156021563,156021792,156022145,156022415,
		# Y    57209886     57212186    +   WASH6P  ENST00000711283 57212186        57212186        57209886,57210639,      57209980,57212186,
		#
		# Example SNORD81:
		# grep SNORD81 $ANNOTSV/share/AnnotSV/Annotations_Human/Genes/GRCh38/genes.ENSEMBL.sorted.bed
		# 1       85592279        85592356        +       SNORD81 ENST00000363064 85592356        85592356        85592279,       85592356,
		# 5       18236122        18236196        +       SNORD81 ENST00000390976 18236196        18236196        18236122,       18236196,
		# 7       137287686       137287762       +       SNORD81 ENST00000365153 137287762       137287762       137287686,      137287762,
		lappend g_L_coord($gene) "$chrom $start $end"
    }
	foreach gene [array names g_L_coord] {
		set g_L_coord($gene) [lsort -unique $g_L_coord($gene)]
	}

    # Then, for gene names not present in the $geneCoordFile (previous symbol, alias), we use the NCBIgeneID to retrieve new alias
    # Memorize the link NCBIgeneID-GeneName
    set NCBIgeneDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Gene-based/NCBIgeneID"
    foreach L [LinesFromFile "$NCBIgeneDir/geneSymbol_NCBIgeneID.tsv"] {
        set gene [lindex $L 0]
        set NCBIgeneID [lindex $L 1]
        set NCBIgeneIDfor($gene) $NCBIgeneID
        lappend L_genesFor($NCBIgeneID) $gene
    }

    # We memorize in "$g_alias" the association of the gene names overlapping the same genomic coordinates
    # Example line (genemap2.txt, GRCh38):
    # chr2    47806920        47906498        2p21    2p16.3  607871  FBXO11, FBX11, VIT1, PRMT9, IDDFBA      F-box only protein 11   FBXO11  80204   ENSG00000138081 ?2p16   Intellectual developmental disorder with dysmorphic facies and behavioral abnormalities, 618089 (3), Autosomal dominant        Fbxo11 (MGI:2147134)
    # => PRMT9 (previous symbol) is located on chromosome 4
    # => FBXO11 is not associated to PRMT9.
    foreach geneToMatch [array names g_L_coord] {

        if {![info exists NCBIgeneIDfor($geneToMatch)]} {continue}
		
		set testLappend 0
		foreach coordTrue $g_L_coord($geneToMatch) {
			set trueChrom [lindex $coordTrue 0]
			set trueStart [lindex $coordTrue 1]
			set trueEnd [lindex $coordTrue 2]
			foreach aliasGene $L_genesFor($NCBIgeneIDfor($geneToMatch)) {
				if {[info exists g_L_coord($aliasGene)]} {
					foreach coord $g_L_coord($aliasGene) {
						set chrom [lindex $coord 0]
						set start [lindex $coord 1]
						set end [lindex $coord 2]
						if {$trueChrom ne $chrom} {continue}
						if {($start < $trueStart && $end > $trueStart) || ($start > $trueStart && $start < $trueEnd)} {
						    lappend g_alias($geneToMatch) $aliasGene
							break
						}
					}
				}
		    }
            if {$testLappend} {break}
		}
	}
	foreach g [array names g_alias] {
		set g_alias($g) [lsort -unique $g_alias($g)]
	}
}



## - Check if the "genemap2.txt" file exists.
## - Check and create if necessary the 'date'_OMIM*annotations.tsv files.
proc checkOMIMfile {} {
    
    global g_AnnotSV
    global g_alias
    global g_L_coord
 

    ## Check if the OMIM file has been downloaded and create the formatted files if needed
    ######################################################################################
    set omimDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Gene-based/OMIM"

    # INFO:
    # Genomic coordinates associated to the OMIM phenotype are available in GRCh38 in "genemap2.txt" ($OMIMfileDownloaded)
    set OMIMfileDownloaded [glob -nocomplain "$omimDir/genemap2.txt"]

    # Header:
    # Chromosome    Genomic Position Start  Genomic Position End    Cyto Location   Computed Cyto Location  MIM Number      Gene/Locus And Other Related Symbols    Gene Name       Approved Gene Symbol    Entrez Gene IDEnsembl Gene ID  Comments        Phenotypes      Mouse Gene Symbol/ID
    
    set OMIMfile1FormattedGzip [glob -nocomplain "$omimDir/*_OMIM-1-annotations.tsv.gz"]
    set OMIMfile2FormattedGzip [glob -nocomplain "$omimDir/*_OMIM-2-annotations.tsv.gz"]
    
    if {$OMIMfileDownloaded eq ""} {
        # No OMIM annotation to do
        return
    }
    
    if {[llength $OMIMfile1FormattedGzip]>1} {
        puts "Several OMIM files exist:"
        puts "$OMIMfile1FormattedGzip"
        puts "Keep only one: [lindex $OMIMfile1FormattedGzip end]\n"
        foreach omim [lrange $OMIMfile1FormattedGzip 0 end-1] {
            file rename -force $omim $omim.notused
        }
        set OMIMfile1FormattedGzip [lindex $OMIMfile1FormattedGzip end]
    }
    if {[llength $OMIMfile2FormattedGzip]>1} {
        puts "Several OMIM files exist:"
        puts "$OMIMfile2FormattedGzip"
        puts "Keep only one: [lindex $OMIMfile2FormattedGzip end]\n"
        foreach omim [lrange $OMIMfile2FormattedGzip 0 end-1] {
            file rename -force $omim $omim.notused
        }
        set OMIMfile2FormattedGzip [lindex $OMIMfile2FormattedGzip end]
    }
    
    if {$OMIMfile1FormattedGzip ne "" && $OMIMfile2FormattedGzip ne ""} {return}
    
    
    ## - Create the 'date'_OMIM-1-annotations.tsv and 'date'_OMIM-2-annotations.tsv files (headers):
    ##   Header1: genes, OMIM_ID
    ##   Header2: genes, OMIM_phenotype, OMIM_inheritance
    
    set OMIMfile1Formatted "$omimDir/[clock format [clock seconds] -format "%Y%m%d"]_OMIM-1-annotations.tsv"
    set OMIMfile2Formatted "$omimDir/[clock format [clock seconds] -format "%Y%m%d"]_OMIM-2-annotations.tsv"
    puts "\t...creation of [file tail $OMIMfile1Formatted.gz] and [file tail $OMIMfile2Formatted] ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    puts "\t   (in [file dirname $OMIMfile1Formatted.gz])"
    ReplaceTextInFile "genes\tOMIM_ID" $OMIMfile1Formatted
    ReplaceTextInFile "genes\tOMIM_phenotype\tOMIM_inheritance" $OMIMfile2Formatted
    

    # Parsing of the GRCh38 $OMIMfileDownloaded
    foreach L [LinesFromFile $OMIMfileDownloaded] {
        
        set Ls [split $L "\t"]
        # Header of "genemap2.txt":
        ###########################
        # # Chromosome
        # Genomic Position Start
        # Genomic Position End
        # Cyto Location
        # Computed Cyto Location
        # MIM Number
        # Gene/Locus And Other Related Symbols
        # Gene Name
        # Approved Gene Symbol
        # Entrez Gene ID
        # Ensembl Gene ID
        # Comments
        # Phenotypes
        # Mouse Gene Symbol/ID
        if {[regexp "^# Chromosome\tGenomic" $L]} {
            set i_chrom [lsearch -exact $Ls "# Chromosome"]
            if {$i_chrom == -1} {
                puts "$OMIMfileDownloaded:\nBad header line syntax. \"# Chromosome\" column not found - Exit with error"
                exit 2
            }
            set i_start [lsearch -exact $Ls "Genomic Position Start"]
            if {$i_start == -1} {
                puts "$OMIMfileDownloaded:\nBad header line syntax. \"Genomic Position Start\" column not found - Exit with error"
                exit 2
            }
            set i_end [lsearch -exact $Ls "Genomic Position End"]
            if {$i_end == -1} {
                puts "$OMIMfileDownloaded:\nBad header line syntax. \"Genomic Position End\" column not found - Exit with error"
                exit 2
            }
            set i_gene1 [lsearch -exact $Ls "Gene/Locus And Other Related Symbols"]
            if {$i_gene1 == -1} {
                puts "$OMIMfileDownloaded:\nBad header line syntax. \"Gene/Locus And Other Related Symbols\" column not found - Exit with error"
                exit 2
            }
            set i_gene2 [lsearch -exact $Ls "Approved Gene Symbol"]
            if {$i_gene2 == -1} {
                puts "$OMIMfileDownloaded:\nBad header line syntax. \"Approved Symbol\" column not found - Exit with error"
                exit 2
            }
            set i_mim [lsearch -exact $Ls "MIM Number"]
            if {$i_mim == -1} {
                puts "$OMIMfileDownloaded:\nBad header line syntax. \"Mim Number\" column not found - Exit with error"
                exit 2
            }
            set i_pheno [lsearch -exact $Ls "Phenotypes"]
            if {$i_pheno == -1} {
                puts "$OMIMfileDownloaded\nBad header line syntax. \"Phenotypes\" column not found - Exit with error"
                exit 2
            }
            continue
        }
        if {[regexp "^#" $L]} {continue}
        
        regsub "chr" [lindex $Ls $i_chrom] "" trueChromOMIM
        set trueStartOMIM [lindex $Ls $i_start]
        set trueEndOMIM [lindex $Ls $i_end]
        
        set gene1 [lindex $Ls $i_gene1]
        regsub -all " " $gene1 "" gene1
        
        set gene2 [lindex $Ls $i_gene2]
        set allGenesTmp [split $gene1 ","]
        if {$gene2 ne ""} {lappend allGenesTmp $gene2}
        
        set allGenesTmp [lsort -unique $allGenesTmp]
        
        # In "genemap2.txt", we only keep the gene names overlapping the given genomic coordinates associated with the OMIM
        # Example line:
        # chr2    47806920        47906498        2p21    2p16.3  607871  FBXO11, FBX11, VIT1, PRMT9, IDDFBA      F-box only protein 11   FBXO11  80204   ENSG00000138081 ?2p16   Intellectual developmental disorder with dysmorphic facies and behavioral abnormalities, 618089 (3), Autosomal dominant        Fbxo11 (MGI:2147134)
        # => PRMT9 (previous symbol) is located on chromosome 4
        # => we don't keep this gene for the OMIM "607871"
        set allGenes ""
        foreach gTmp $allGenesTmp {
            if {$gTmp eq ""} {continue}
            
            if {[info exists g_L_coord($gTmp)]} {
                foreach coord $g_L_coord($gTmp) {
					set chrom [lindex $coord 0]
					set start [lindex $coord 1]
					set end [lindex $coord 2]
					if {$trueChromOMIM ne $chrom} {continue}
					if {($start < $trueStartOMIM && $end > $trueStartOMIM) || ($start > $trueStartOMIM && $start < $trueEndOMIM)} {
					    lappend allGenes $gTmp
						# Add the alias for genes validated with coordinates
			            if {[info exists g_alias($gTmp)]} {
						     lappend allGenes {*}$g_alias($gTmp)
					    }
						break
					}
				}	
            }
        }
        if {$allGenes eq ""} {continue}
        set allGenes [lsort -unique $allGenes]

        set mim [lindex $Ls $i_mim]
        set phenoTmp [lindex $Ls $i_pheno]
        regsub -all -nocase "Autosomal dominant" $phenoTmp "AD" phenoTmp
        regsub -all -nocase "Autosomal recessive" $phenoTmp "AR" phenoTmp
        regsub -all -nocase "X-linked dominant" $phenoTmp "XLD" phenoTmp
        regsub -all -nocase "X-linked recessive" $phenoTmp "XLR" phenoTmp
        regsub -all -nocase "Y-linked dominant" $phenoTmp "YLD" phenoTmp
        regsub -all -nocase "Y-linked recessive" $phenoTmp "YLR" phenoTmp
        regsub -all -nocase "X-linked" $phenoTmp "XL" phenoTmp
        regsub -all -nocase "Y-linked" $phenoTmp "YL" phenoTmp
        
        foreach p [split $phenoTmp ";"] {
            if {[regexp " *(.*\\\(\[0-9\]+\\\)),? *(.*)" $p match 1pheno 1inher]} {
                if {$1inher ne ""} {
                    regsub -all " *, +" $1inher "," 1inher
                    foreach g $allGenes {
                        regsub -all " " $g "" g
                        if {$g eq ""} {continue}
                        lappend Pheno($g) "$1pheno $1inher"
                        lappend Inheritance($g) {*}[split $1inher ";"]
                    }
                } else {
                    foreach g $allGenes {
                        regsub -all " " $g "" g
                        if {$g eq ""} {continue}
                        lappend Pheno($g) "$1pheno"
                    }
                }
                foreach g $allGenes {
                    regsub -all " " $g "" g
                    if {$g eq ""} {continue}
                    lappend L_genes $g
                    lappend Mim($g) $mim
                }
            }
        }
    }
    
    # creation of $OMIMfileFormatted.gz
    set L_genes [lsort -unique $L_genes]
    # Header :
    # Gene    OMIM_ID
    # Gene    OMIM_phenotype      OMIM_inheritance
    foreach g $L_genes {
        set Pheno($g) "[join [lsort -unique $Pheno($g)] ";"]"
        if {![info exists Inheritance($g)]} {set Inheritance($g) ""}
        set Inheritance($g) "[join [lsort -unique $Inheritance($g)] ";"]"
        if {[regexp "^;+$" $Inheritance($g)]} {set Inheritance($g) ""}
        WriteTextInFile "$g\t[join [lsort -unique $Mim($g)] ";"]" $OMIMfile1Formatted
        WriteTextInFile "$g\t$Pheno($g)\t$Inheritance($g)" $OMIMfile2Formatted
    }
    if {[catch {exec gzip $OMIMfile1Formatted} Message]} {
        puts "-- checkOMIMfile --"
        puts "gzip $OMIMfile1Formatted"
        puts "$Message\n"
    }
    if {[catch {exec gzip $OMIMfile2Formatted} Message]} {
        puts "-- checkOMIMfile --"
        puts "gzip $OMIMfile2Formatted"
        puts "$Message\n"
    }
    
    file delete -force $OMIMfileDownloaded
}

## - Check if the "morbidmap.txt" file exists.
## - Check and create if necessary the 'date'_morbid.tsv.gz file.
## - Check and create if necessary the 'date'_morbidCandidate.tsv.gz file.
proc checkMorbidfile {} {
    
    global g_AnnotSV
    global g_alias
	global g_L_coord 
    
    set omimDir "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Gene-based/OMIM"
    
    ## Check if the Morbid file has been downloaded the formatted
    ##################################################################
    set MorbidFileDownloaded [glob -nocomplain "$omimDir/morbidmap.txt"]
    # Header:
    # # Phenotype     Gene/Locus And Other Related Symbols    MIM Number      Cyto Location
    set MorbidFileFormattedGzip [glob -nocomplain "$omimDir/*_morbid.tsv.gz"]
    set MorbidCandidateFileFormattedGzip [glob -nocomplain "$omimDir/*_morbidCandidate.tsv.gz"]
    
    if {$MorbidFileDownloaded eq "" && $MorbidFileFormattedGzip eq ""} {
        # No Morbid annotation
        return
    }
    
    if {[llength $MorbidFileFormattedGzip]>1} {
        puts "Several Morbid Genes files exist:"
        puts "$MorbidFileFormattedGzip"
        puts "Keep only one: [lindex $MorbidFileFormattedGzip end]\n"
        foreach morbid [lrange $MorbidFileFormattedGzip 0 end-1] {
            file rename -force $morbid $morbid.notused
        }
        return
    }
    if {[llength $MorbidCandidateFileFormattedGzip]>1} {
        puts "Several Morbid Genes candidate files exist:"
        puts "$MorbidCandidateFileFormattedGzip"
        puts "Keep only one: [lindex $MorbidCandidateFileFormattedGzip end]\n"
        foreach morbid [lrange $MorbidCandidateFileFormattedGzip 0 end-1] {
            file rename -force $morbid $morbid.notused
        }
        return
    }
    
    if {$MorbidFileFormattedGzip eq ""} {
        ## - Create the 'date'_morbid.tsv file.
        ##   Header: genes, OMIM_morbid
        ## - Create the 'date'_morbidCandidate.tsv.gz file.
        ##   Header: genes, OMIM_morbid_candidate
       
        set MorbidFileFormatted "$omimDir/[clock format [clock seconds] -format "%Y%m%d"]_morbid.tsv"
        set MorbidCandidateFileFormatted "$omimDir/[clock format [clock seconds] -format "%Y%m%d"]_morbidCandidate.tsv"
        puts "\t...creation of [file tail $MorbidFileFormatted.gz] and [file tail $MorbidCandidateFileFormatted] ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
        puts "\t   (in [file dirname $MorbidFileFormatted.gz])"
        ReplaceTextInFile "genes\tOMIM_morbid" $MorbidFileFormatted
        ReplaceTextInFile "genes\tOMIM_morbid_candidate" $MorbidCandidateFileFormatted
        
        # Parsing of $MorbidFileDownloaded
        set L_morbid {}
        set L_morbidCandidate {}
        foreach L [LinesFromFile $MorbidFileDownloaded] {
            set Ls [split $L "\t"]
            
            # Header from the downloaded file:
            # Phenotype     Gene Symbols    MIM Number      Cyto Location
            if {[regexp "^# Phenotype\t" $L]} {
                set i_gene  [lsearch $Ls "Gene/Locus And Other Related Symbols"]
                if {$i_gene == -1} {
                    puts "Bad header line syntax. \"Gene/Locus And Other Related Symbols\" column not found - Exit with error"
                    exit 2
                }
                set i_pheno [lsearch $Ls "# Phenotype"]
                if {$i_pheno == -1} {
                    puts "Bad header line syntax. Phenotype column not found - Exit with error"
                    exit 2
                }
                continue
            }
            if {[regexp "^#" $L]} {continue}
            
            set pheno [lindex $Ls $i_pheno]
            
            ## "[ ]", indicate "nondiseases", mainly genetic variations that lead to apparently abnormal laboratory test values (e.g., dysalbuminemic euthyroidal hyperthyroxinemia).
            if {[regexp "^\\\[.+\\\]" $pheno]} {continue}
            
            ## (1) The disorder is placed on the map based on its association with a gene, but the underlying defect is not known.
            ## (2) The disorder has been placed on the map by linkage or other statistical method; no mutation has been found.
            if {[regexp "\\\(1\\\)|\\\(2\\\)" $pheno]} {continue}
            
            ## The following lines are automatically either (3) or (4):
            ## (3) indicates that the molecular basis of the disorder is known; a mutation has been found in the gene.
            ## (4) indicates that a contiguous gene deletion or duplication syndrome, multiple genes are deleted or duplicated causing the phenotype.
            set gene [lindex $Ls $i_gene]
            set allGenes [split $gene ","]
            set allGenes [lsort -unique $allGenes]
            foreach g $allGenes {
                regsub -all " " $g "" g
                if {$g eq ""} {continue}
                ## "{ }", indicate mutations that contribute to susceptibility to multifactorial disorders (e.g., diabetes, asthma) or to susceptibility to infection (e.g., malaria).
                ## "?", before the phenotype name indicates that the relationship between the phenotype and gene is provisional.

                # Add the gene names (approved symbol + other related symbol)
                if {[regexp "^{.+}|^\\\?" $pheno]} {
					if {[info exists g_alias($g)]} {
						lappend L_morbidCandidate {*}$g_alias($g)
					} else {
						lappend L_morbidCandidate $g
					}
                } else {
                    if {[info exists g_alias($g)]} {
						lappend L_morbid {*}$g_alias($g)
                    } else {
                        lappend L_morbid $g
                    }
                }
            }
        }
        
        # creation of $MorbidFileFormatted.gz
        set L_morbid [lsort -unique $L_morbid]

        ## created header: genes, OMIM_morbid
        foreach g $L_morbid {
            WriteTextInFile "$g\tyes" $MorbidFileFormatted
        }
        
		if {[catch {exec gzip $MorbidFileFormatted} Message]} {
            puts "-- checkMorbidfile --"
            puts "gzip $MorbidFileFormatted"
            puts "$Message\n"
        }
        
        # creation of $MorbidCandidateFileFormatted.gz
        set L_morbidCandidate [lsort -unique $L_morbidCandidate]
        ## creted header: genes, OMIM_morbid_candidate
        foreach g $L_morbidCandidate {
            WriteTextInFile "$g\tyes" $MorbidCandidateFileFormatted
        }
        if {[catch {exec gzip $MorbidCandidateFileFormatted} Message]} {
            puts "-- checkMorbidfile --"
            puts "gzip $MorbidCandidateFileFormatted"
            puts "$Message\n"
        }
        
        file delete -force $MorbidFileDownloaded
    }
	
	# Defined in memorizeGeneNameAlias
    unset g_alias
	unset g_L_coord
}


proc fromOMIMtoPhenotype {OMIM} {
    global g_AnnotSV
    global g_phen
    global g_gene
    
    if { ! [info exists g_phen(DONE)]} {
        set OMIMfile1 [glob -nocomplain "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Gene-based/OMIM/*_OMIM-1-annotations.tsv.gz"]
        # 20201107_OMIM-1-annotations.tsv.gz
        # Header = genes   OMIM_ID
        foreach L [LinesFromGZFile $OMIMfile1] {
            set Ls [split $L "\t"]
            set gene [lindex $Ls 0]
            set L_omim [split [lindex $Ls 1] ";"]
            foreach omim $L_omim {
                lappend g_gene($omim) $gene
            }
        }
        
        set OMIMfile2 [glob -nocomplain "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Gene-based/OMIM/*_OMIM-2-annotations.tsv.gz"]
        # 20201107_OMIM-2-annotations.tsv.gz
        # Header = genes   OMIM_phenotype      OMIM_inheritance
        foreach L [LinesFromGZFile $OMIMfile2] {
            set Ls [split $L "\t"]
            set gene [lindex $Ls 0]
            set phenotype [lindex $Ls 1]
            regsub -all "\\\{|\\\}|\\\[|\\\]" $phenotype "" phenotype
            regsub -all "/" $phenotype ";" phenotype
            set g_phen($gene) $phenotype
        }
        
        set g_phen(DONE) 1
    }
    
    set phenotype ""
    if {[info exists g_gene($OMIM)]} {
        foreach gene $g_gene($OMIM) {
            if {[info exists g_phen($gene)]} {
                lappend phenotype $g_phen($gene)
            }
        }
        set phenotype [join $phenotype ";"]
    }
    
    return $phenotype
}

proc fromOMIMgeneToPhenotype {OMIMgene} {
    global g_AnnotSV
    global g_phen
    global g_gene
    
    if { ! [info exists g_phen(DONE)]} {
        fromOMIMtoPhenotype 617631 ;# use to initialize "g_phen" and "g_gene"
    }
    
    if {[info exists g_phen($OMIMgene)]} {
        return $g_phen($OMIMgene)
    } else {
        return ""
    }
}


proc isMorbid {gene} {
    global g_AnnotSV
    global g_morbid
    
    if {![info exists g_morbid(DONE)]} {
        set g_morbid(DONE) 1
        set morbidFile [glob -nocomplain "$g_AnnotSV(annotationsDir)/Annotations_$g_AnnotSV(organism)/Gene-based/OMIM/*_morbid.tsv.gz"]
        if {$morbidFile ne ""} { ;# Doesn't exists for Mouse
            set f [open "| gzip -cd $morbidFile"]
            while {! [eof $f]} {
                set L [gets $f]
                set g_morbid([lindex $L 0]) 1
            }
            close $f
        }
    }
    
    if {[info exists g_morbid($gene)]} {
        return 1
    } else {
        return 0
    }
}

