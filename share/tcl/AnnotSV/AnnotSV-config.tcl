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


## Loading the different options in the following order:
## - Default options
## - Config file options (if file exists)
## - Options given in arguments
#
## Please note: the case in the name of options is not important. Ex: "vcffiles" = "vcfFiles"
##
proc configureAnnotSV {argv} {
    
    global g_AnnotSV
    
    puts "...downloading the configuration data ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    
    #######################
    ## Load default options
    #######################
    puts "\t...configuration data by default"
    set g_AnnotSV(annotationsDir)           ""
    set g_AnnotSV(bedtools)                 "bedtools"
    set g_AnnotSV(candidateGenesFile)       ""
    set g_AnnotSV(candidateGenesFiltering)  "no"
    set g_AnnotSV(extann)                   ""    ;# list of the “.../Annotations_$g_AnnotSV(organism)/*/*.tsv(.gz) files” <=> External genes annotation files
    set g_AnnotSV(filteredVCFfiles)         ""
    set g_AnnotSV(filteredVCFsamples)       ""
    set g_AnnotSV(genomeBuild)              "GRCh37"
    set g_AnnotSV(metrics)                  "us"
    set g_AnnotSV(minTotalNumber)           "500"
    set g_AnnotSV(outputDir)                ""
    set g_AnnotSV(outputFile)               ""
    set g_AnnotSV(overlap)                  "70"
    set g_AnnotSV(overwrite)                "yes"
    set g_AnnotSV(promoterSize)             "500"
    set g_AnnotSV(rankFiltering)            "1-5"
    set g_AnnotSV(rankOutput)               "no"
    set g_AnnotSV(reciprocal)               "no"
    set g_AnnotSV(SVinputFile)              ""
    set g_AnnotSV(SVinputInfo)              "1"
    set g_AnnotSV(SVminSize)                "50"
    set g_AnnotSV(svtBEDcol)                "-1"
    set g_AnnotSV(txFile)                   ""
    set g_AnnotSV(typeOfAnnotation)         "both"
    set g_AnnotSV(vcfFiles)                 ""
    set g_AnnotSV(vcfPASS)                  "0"
    set g_AnnotSV(vcfSamples)               ""
    set g_AnnotSV(ranking)                  "1"
    set g_AnnotSV(outputColHeader)          ""

    ###########################
    ## Load config file options
    ###########################
    set lOptionsOk "annotationsDir bedtools candidateGenesFile candidateGenesFiltering extann filteredVCFfiles filteredVCFsamples genomeBuild metrics minTotalNumber outputDir outputFile overlap overwrite promoterSize rankFiltering rankOutput reciprocal SVinputFile SVinputInfo SVminSize svtBEDcol txFile typeOfAnnotation vcfFiles vcfSamples vcfPASS"
    set configFile "$g_AnnotSV(etcDir)/configfile"
    if {[file exists "./configfile"]} {
	set configFile "./configfile"
    }
    puts "\t...configuration data from $configFile"
    set testColumnNames 0
    foreach L [LinesFromFile $configFile] {
	if {[regexp "^# AnnotSV Output columns:" $L]} {set testColumnNames 1} 
	if {[regexp "^#" $L]} {continue}
	if {$L eq ""} {continue}
	#Reading the config file 
	if { ! $testColumnNames} {
	    regsub -all "^-|:" $L "" L
	    set optionName  [lindex $L 0]
	    set optionValue [lindex $L 1]
	    set k [lsearch -exact -nocase $lOptionsOk $optionName]
	    if {$k != -1} {
		set optionName [lindex $lOptionsOk $k]
		set g_AnnotSV($optionName) $optionValue
	    } else {
		puts "############################################################################"
		puts "\"$optionName\" option not known."
		puts "For more information on the arguments, please use the -help option"
		puts "############################################################################"
		exit 2
	    }
	} else {
	    regsub "( |\t|\\*)+$" $L "" L
	    lappend g_AnnotSV(outputColHeader) $L
	}
    }		
    
    ##################################
    ## Load options given in arguments
    ##################################
    puts "\t...configuration data given in arguments"
    regsub -all "^-|:" $argv "" argv
    set i 0
    set j 1
    while {$j < [llength $argv]} {
	set optionName [lindex $argv $i]
	regsub -all "^-|:" $optionName "" optionName
	
	set optionValue [lindex $argv $j]
	set  k [lsearch -exact -nocase $lOptionsOk $optionName]
	if {$k != -1} {
	    set optionName [lindex $lOptionsOk $k]
	    set g_AnnotSV($optionName) $optionValue
	} else {
	    puts "\"$optionName\" option not known."
	    puts "For more information on the arguments, please use the -help option"
	    exit 2
	}
	
	incr i 2
	incr j 2
    }
    

    ########################################
    ## Checking of the configuration options
    ########################################
    puts "\t...checking configuration data and files\n"

    ## annotationsDir: It must be an existing directory (or "" for the default)
    if {$g_AnnotSV(annotationsDir) eq ""} {
	set g_AnnotSV(annotationsDir) "$g_AnnotSV(installDir)/share/AnnotSV"
    } elseif {![file isdirectory $g_AnnotSV(annotationsDir)]} {
	puts "AnnotSV needs in argument an existing path of the annotations directory (-annotationsDir = \"$g_AnnotSV(annotationsDir)\") - Exit with error."
	exit 2
    }

    ## SVinputFile: We should have a bed or VCF input file
    if {$g_AnnotSV(SVinputFile) eq ""} {
	puts "AnnotSV needs in argument the path of your bed or VCF input file (-SVinputFile ...) - Exit with error."
	exit 2
    }
    if {![regexp -nocase "\\.(bed|vcf(.gz)?)$" $g_AnnotSV(SVinputFile)]} {
	puts "############################################################################"
	puts "Bad option value: -SVinputFile = $g_AnnotSV(SVinputFile)"
	puts "Extension file should be \".bed\" or \".vcf\" - Exit with error."
	puts "############################################################################"
	exit 2
    }
    ## SVinputFile: It must be an existing file 
    if {![file exists $g_AnnotSV(SVinputFile)]} {
	puts "############################################################################"
	puts "Bad value for the SVinputFile option, file does not exist ($g_AnnotSV(SVinputFile)) - Exit with error."
	puts "############################################################################"
	exit 2
    }
    ## SVinputFile: no SV to annotate if it is an empty file
    ## (could happen if AnnotSV is implemented into a pipeline)
    if {[isAnEmptyFile $g_AnnotSV(SVinputFile)]} {
	puts "############################################################################"
	puts "SVinputFile ($g_AnnotSV(SVinputFile) is empty, no SV to annotate - Exit without error."
	puts "############################################################################"
	exit 0
    }

    ## metrics: us or fr
    if {$g_AnnotSV(metrics) ne "us" && $g_AnnotSV(metrics) ne "fr"} {
	puts "############################################################################"
	puts "Bad option value: -metrics = $g_AnnotSV(metrics)"
	puts "Should be \"us\" or \"fr\" - Continue with the \"us\" default."
	puts "############################################################################"
	set g_AnnotSV(metrics) "us"
    }

    ## It must be: yes or no
    foreach val {candidateGenesFiltering rankOutput} {
	if {$g_AnnotSV($val) ne "yes" && $g_AnnotSV($val) ne "no"} {
	    puts "############################################################################"
	    puts "Bad option value: -$val = $g_AnnotSV($val)"
	    puts "Should be \"no\" or \"yes\""
	    puts "############################################################################"
	    exit 2
	}
    }

    ## It must be a list between 1 and 5 
    ## e.g.: "3,4,5" or "3-5" 
    set liste {}
    while {[regexp "(\[1-5\])-(\[1-5\])" $g_AnnotSV(rankFiltering) match i j]} {
	while {$i <= $j} {
	    lappend liste $i
	    incr i
	}
	regsub "\[1-5\]-\[1-5\]" $g_AnnotSV(rankFiltering) "" g_AnnotSV(rankFiltering)
    }
    foreach i [split $g_AnnotSV(rankFiltering) ","] {
	if {[lsearch -exact {1 2 3 4 5} $i] ne -1} {lappend liste $i}	    
    }
    set liste [lsort -unique $liste]
    set g_AnnotSV(rankFiltering) $liste

    ## It must be a boolean (0 or 1) for the SVinputInfo and vcfPASS options.
    foreach val {SVinputInfo vcfPASS} {
	if {$g_AnnotSV($val) ne 0 && $g_AnnotSV($val) ne 1} {
	    puts "############################################################################"
	    puts "Bad option value: -$val = $g_AnnotSV($val)"
	    puts "Should be \"0\" or \"1\""
	    puts "############################################################################"
	    exit 2
	}
    }

    ## g_AnnotSV(outputDir) must be an existing directory. Else, we create it
    ## g_AnnotSV(outputFile) must be defined and should not already exists.
    ## g_AnnotSV(outputFile) extension must be ".tsv".
    if {$g_AnnotSV(outputDir) eq "" && [regexp "/" "$g_AnnotSV(outputFile)"]} {
	# Info: Directory and name file could be both given in the g_AnnotSV(outputFile) option
	set g_AnnotSV(outputDir) [file dirname "$g_AnnotSV(outputFile)"]
    }
    if {$g_AnnotSV(outputDir) ne ""} {
	if {![file exists $g_AnnotSV(outputDir)]} {
	    file mkdir "$g_AnnotSV(outputDir)"
	}
    } else {
	set g_AnnotSV(outputDir) "[clock format [clock seconds] -format "%Y%m%d"]_AnnotSV"
	file mkdir "$g_AnnotSV(outputDir)"
    }
    if {$g_AnnotSV(outputFile) ne ""} {
	set g_AnnotSV(outputFile) [file tail $g_AnnotSV(outputFile)]
	if {[string range $g_AnnotSV(outputFile) end-3 end] ne ".tsv"} {
	    set g_AnnotSV(outputFile) "$g_AnnotSV(outputFile).tsv"
	}
    } else {
	regsub -nocase "(\\.bed|\\.vcf(\\.gz)?)$" [file tail $g_AnnotSV(SVinputFile)] ".annotated.tsv" g_AnnotSV(outputFile)
    }
    if {$g_AnnotSV(overwrite) eq "no" && [file exists $g_AnnotSV(outputDir)/$g_AnnotSV(outputFile)]} {
	puts "############################################################################"
	puts "Bad value for the -outputFile option, file already exists ($g_AnnotSV(outputDir)/$g_AnnotSV(outputFile)) - Exit with error."
	puts "############################################################################"
	exit 2
    }

    # bedtools: It should be a good path that we can run
    if {[catch {exec $g_AnnotSV(bedtools)} Message]} {
	puts "############################################################################"
	puts "Bad value for the bedtools option ($g_AnnotSV(bedtools))"
	puts "$Message"
	puts "Exit with error."
	puts "############################################################################"	
	exit 2
    }
    
    ## txFile: It must be an existing file 
    if {$g_AnnotSV(txFile) ne "" && ![file exists $g_AnnotSV(txFile)]} {
	puts "############################################################################"
	puts "Bad value for the txFile option, file does not exist ($g_AnnotSV(txFile)) - Exit with error."
	puts "############################################################################"
	exit 2
    }

    ## candidateGenesFile: It must be an existing file 
    if {$g_AnnotSV(candidateGenesFile) ne "" && ![file exists $g_AnnotSV(candidateGenesFile)]} {
	puts "############################################################################"
	puts "Bad value for the candidateGenesFile option, file does not exist ($g_AnnotSV(candidateGenesFile)) - Exit with error."
	puts "############################################################################"
	exit 2
    }

    ## It must be an integer comprised between 0 and 100 
    if {[regexp "\[^0-9\]" $g_AnnotSV(overlap)] || $g_AnnotSV(overlap)<=0 || $g_AnnotSV(overlap)>100} {
	puts "############################################################################"
	puts "Bad value for the overlap option ($g_AnnotSV(overlap)), should be an integer ]0,100] - Exit with error."
	puts "############################################################################"
	exit 2
    }

    ## It must be an integer comprised between 100 and 1000
    if {[regexp "\[^0-9\]" $g_AnnotSV(minTotalNumber)] || $g_AnnotSV(minTotalNumber)<100 || $g_AnnotSV(minTotalNumber)>1000} {
	puts "############################################################################"
	puts "Bad value for the $val option ($g_AnnotSV($val)), should be an integer comprised between 100 and 1000 - Exit with error."
	puts "############################################################################"
	exit 2
    }

    ## It must be "no" or "yes" for the reciprocal/overwrite option.
    foreach val {reciprocal overwrite} {
	if {![regexp -nocase "^(no)|(yes)$" $g_AnnotSV($val)]} {
	    puts "############################################################################"
	    puts "Bad option value: -$val = $g_AnnotSV($val)"
	    puts "Should be \"no\" or \"yes\""
	    puts "############################################################################"
	    exit 2
	}
    }

    ## It must be an integer > 0
    foreach val {promoterSize SVminSize} {
	if {[regexp "\[^0-9\]" $g_AnnotSV($val)]} {
	    puts "############################################################################"
	    puts "Bad value for the $val option ($g_AnnotSV($val)), should be a positive integer - Exit with error."
	    puts "############################################################################"
	    exit 2
	}
    }

    ## It must be "both", "full" or "split"
    set L_typeOfAnnotation {both full split}
    if {[lsearch -exact $L_typeOfAnnotation "$g_AnnotSV(typeOfAnnotation)"] eq -1} {
	puts "############################################################################"
	puts "Bad option value: -typeOfAnnotation = $g_AnnotSV(typeOfAnnotation)"
	puts "Should be \"both\", \"full\" or \"split\""
	puts "############################################################################"
	exit 2
    }

    ## The following step could be improved: too long
    #################################################

    ## If "vcfFiles" option is defined:
    if {$g_AnnotSV(vcfFiles) ne ""} {
	## It must be existing files
	## All "$vcfSamples" should be presents in one of the VCF files
	## If the "-vcfSamples" option is not defined, we defined it with all samples from the VCF files.
	set L_samples {}
	set L_allSamplesFromVCF {}

	# Script à améliorer ici, un nom de fichier (sans expr reg) qui n'existe pas ne provoque pas de message d'erreur
	foreach vcfF [eval glob -nocomplain $g_AnnotSV(vcfFiles)] {
	    if {![file exists $vcfF]} { ;# Never used during the script. If the file doesn't exist, the foreach is empty! 
		puts "############################################################################"
		puts "Bad value for vcfFiles option, file does not exist ($vcfF) - Exit with error."
		puts "############################################################################"
		exit 2
	    }
	    if {[regexp ".gz$" $vcfF]} {
		set f [open "|gzip -cd $vcfF"]
	    } else {
		set f [open "$vcfF"]
	    }
	    while {![eof $f]} {
		set L [gets $f]
		if {[string range $L 0 5] ne "#CHROM"} {continue}
		set Ls [split $L "\t"]
		lappend L_allSamplesFromVCF {*}[lrange $Ls 9 end]
		foreach sample $g_AnnotSV(vcfSamples) {
		    if {[lsearch -exact $Ls $sample] ne -1} {lappend L_samples $sample}
		}
		# We can't put a break here, because the pipe command: "|gzip ..." still produces data.
		# break will creates a broken pipe signal: 
		#   child killed: write on pipe with no readers
		#   while executing "close $f"
	    }
	    close $f
	}
	if {$g_AnnotSV(vcfSamples) eq ""} {
	    set L_allSamplesFromVCF [lsort -unique $L_allSamplesFromVCF]
	    set g_AnnotSV(vcfSamples) $L_allSamplesFromVCF
	} else {
	    set L_samples [lsort -unique $L_samples]
	    set g_AnnotSV(vcfSamples) $L_samples; # Remove samples from g_AnnotSV(vcfSamples) that are not in a VCF file
	}
    }

    ## If "filteredVCFfiles" option is defined:
    if {$g_AnnotSV(filteredVCFfiles) ne ""} {
	## It must be existing files
	## All "$filteredVCFsamples" should be presents in one of the VCF files
	## If the "-filteredVCFsamples" option is not defined, we defined it with all samples from the VCF files.
	set L_samples {}
	set L_allSamplesFromVCF {}
	foreach vcfF [eval glob -nocomplain $g_AnnotSV(filteredVCFfiles)] {
	    if {![file exists $vcfF]} {
		puts "############################################################################"
		puts "Bad value for filteredVCFfiles option, file does not exist ($vcfF) - Exit with error."
		puts "############################################################################"
		exit 2
	    }
	    if {[regexp ".gz$" $vcfF]} {
		set f [open "|gzip -cd $vcfF"]
	    } else {
		set f [open "$vcfF"]
	    }
	    while {![eof $f]} {
		set L [gets $f]
		if {[string range $L 0 5] ne "#CHROM"} {continue}
		set Ls [split $L "\t"]
		lappend L_allSamplesFromVCF {*}[lrange $Ls 9 end]
		foreach sample $g_AnnotSV(filteredVCFsamples) {
		    if {[lsearch -exact $Ls $sample] ne -1} {lappend L_samples $sample}
		}
		# We can't put a break here, because the pipe command: "|gzip ..." still produces data.
		# break will creates a broken pipe signal: 
		#   child killed: write on pipe with no readers
		#   while executing "close $f"
	    }
	    close $f
	}
	if {$g_AnnotSV(filteredVCFsamples) eq ""} {
	    set L_allSamplesFromVCF [lsort -unique $L_allSamplesFromVCF]
	    set g_AnnotSV(filteredVCFsamples) $L_allSamplesFromVCF
	} else {
	    set L_samples [lsort -unique $L_samples]
	    set g_AnnotSV(filteredVCFsamples) $L_samples; # Remove samples from g_AnnotSV(filteredVCFsamples) that are not in a VCF file
	}
    }
    
    ## It must be "GRCh37" or "GRCh38" or "mm9" or "mm10" for the genomeBuild option.
    if {![regexp -nocase "^(GRCh37)|(GRCh38)|(mm9)|(mm10)$" $g_AnnotSV(genomeBuild)]} {
	puts "############################################################################"
	puts "Bad option value: -genomeBuild = $g_AnnotSV(genomeBuild)"
	puts "Should be \"GRCh37\", \"GRCh38\", \"mm9\" or \"mm10\""
	puts "############################################################################"
	exit 2
    }
    # Definition of the $g_AnnotSV(organism) variable
    if {[regexp -nocase "^(GRCh37)|(GRCh38)$" $g_AnnotSV(genomeBuild)]} {
	set g_AnnotSV(organism) "Human"
    } elseif {[regexp -nocase "^(mm9)|(mm10)$" $g_AnnotSV(genomeBuild)]} {
	set g_AnnotSV(organism) "Mouse"
    } 

    # Some annotation columns are essential for the ranking: can not be removed by the user
    set g_AnnotSV(genesBasedAnn) 1
    foreach col "DGV_GAIN_n_samples_tested DGV_GAIN_Frequency DGV_LOSS_n_samples_tested DGV_LOSS_Frequency dbVar_event dbVar_status morbidGenes morbidGenesCandidates GHgene_elite GHgene_not_elite pLI_ExAC HI_CGscore TriS_CGscore" {
	if {[lsearch -exact "$g_AnnotSV(outputColHeader)" $col] eq -1} {
	    lappend g_AnnotSV(outputColHeader) $col
	}
    }
 
    return	
}

proc createBEDinputHeaderFile {} {
    
    global g_AnnotSV
    global headerFileToRemove 

    ## SVinputfile is a BED
    regsub -nocase ".bed$" $g_AnnotSV(bedFile) ".header.tsv" BEDinputHeaderFile

    ## The header file doesn't exist
    if {![file exists $BEDinputHeaderFile]} {
	set f [open $g_AnnotSV(bedFile)]
	while {![eof $f]} {
	    set L [gets $f]
	    if {$L eq ""} {continue}
	    if {[regexp "^#" $L]} {
		set header $L
		continue
	    }
	    break
	}
	close $f
	if {[info exists header]} {
	    WriteTextInFile "$header" $BEDinputHeaderFile
	    set headerFileToRemove 1
	}
    }
    return
}
