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


## Loading the different options in the following order:
## - Default options
## - Config file options (if file exists)
## - Options given in arguments
#
## Please note: the case in the name of options is not important. Ex: "snvIndelFiles" = "snvindelfiles"
##
proc configureAnnotSV {argv} {
    
    global g_AnnotSV
    
    puts "...downloading the configuration data ([clock format [clock seconds] -format "%B %d %Y - %H:%M"])"
    
    #######################
    ## Load default options
    #######################
    puts "\t...configuration data by default"
    set g_AnnotSV(annotationsDir)           ""
    set g_AnnotSV(annotationMode)           "both"
    set g_AnnotSV(bcftools)                 "bcftools"
    set g_AnnotSV(bedtools)                 "bedtools"
    set g_AnnotSV(benignAF)                 "0.01"
    set g_AnnotSV(candidateGenesFile)       ""
    set g_AnnotSV(candidateGenesFiltering)  "0"
    set g_AnnotSV(candidateSnvIndelFiles)   ""
    set g_AnnotSV(candidateSnvIndelSamples) ""
    set g_AnnotSV(extann)                   ""    ;# list of the “.../Annotations_$g_AnnotSV(organism)/*/*.tsv(.gz) files” <=> External gene annotation files
    set g_AnnotSV(externalGeneFiles)        ""
    set g_AnnotSV(genomeBuild)              "GRCh38"
    set g_AnnotSV(hpo)                      ""    ;# "HP:0030684,HP:0085622"
    set g_AnnotSV(includeCI)                "1"
    set g_AnnotSV(metrics)                  "us"
    set g_AnnotSV(minTotalNumber)           "500"
    set g_AnnotSV(outputColHeader)          ""    ;# not given in parameter
    set g_AnnotSV(outputDir)                ""
    set g_AnnotSV(outputFile)               ""
    set g_AnnotSV(overlap)                  "100"
    set g_AnnotSV(overwrite)                "1"
    set g_AnnotSV(promoterSize)             "500"
    set g_AnnotSV(rankFiltering)            "1-5,NA"
    set g_AnnotSV(ranking)                  "1"   ;# not given in parameter
    set g_AnnotSV(reciprocal)               "0"
    set g_AnnotSV(REreport)                 "0"
    set g_AnnotSV(REselect1)                "1"
    set g_AnnotSV(REselect2)                "1"
    set g_AnnotSV(samplesidBEDcol)          "-1"
    set g_AnnotSV(samplesidTSVcol)          "-1"  ;# not given in parameter
    set g_AnnotSV(snvIndelFiles)            ""
    set g_AnnotSV(snvIndelPASS)             "0"
    set g_AnnotSV(snvIndelSamples)          ""
    set g_AnnotSV(SVinputFile)              ""
    set g_AnnotSV(SVinputInfo)              "1"
    set g_AnnotSV(SVminSize)                "50"
    set g_AnnotSV(svtBEDcol)                "-1"
    set g_AnnotSV(svtTSVcol)                "-1"  ;# not given in parameter
    set g_AnnotSV(tx)                       "RefSeq"
    set g_AnnotSV(txFile)                   ""
    set g_AnnotSV(variantconvertDir)        ""
    set g_AnnotSV(vcf)                      "0"
    
    
    ###########################
    ## Load config file options
    ###########################
    set lOptionsOk "annotationsDir annotationMode bcftools bedtools benignAF candidateGenesFile candidateGenesFiltering candidateSnvIndelFiles candidateSnvIndelSamples extann externalGeneFiles genomeBuild hpo includeCI metrics minTotalNumber outputDir outputFile overlap overwrite promoterSize rankFiltering reciprocal REreport REselect1 REselect2 samplesidBEDcol snvIndelFiles snvIndelPASS snvIndelSamples SVinputFile SVinputInfo SVminSize svtBEDcol tx txFile variantconvertDir vcf"
    
    # Setting of $g_AnnotSV(SVinputFile) from the command line
    set i 0
    set j 1
    while {$j < [llength $argv]} {
        set optionName [lindex $argv $i]
        regsub "^-|:\[ \t\]" $optionName "" optionName
        if {$optionName eq "SVinputFile"} {
            set g_AnnotSV(SVinputFile) [lindex $argv $j]
        }
        incr i 2
        incr j 2
    }
    
    set configFile "$g_AnnotSV(etcDir)/configfile"
    if {[file exists "[file dirname $g_AnnotSV(SVinputFile)]/configfile"]} {
        set configFile "[file dirname $g_AnnotSV(SVinputFile)]/configfile"
    }
    puts "\t...configuration data from $configFile"
    set testColumnNames 0
    foreach L [LinesFromFile $configFile] {
        if {[regexp "^# AnnotSV Output columns:" $L]} {set testColumnNames 1}
        if {[regexp "^#" $L]} {continue}
        if {$L eq ""} {continue}
        if { ! $testColumnNames} {
            # Reading the configfile options (not the column names)
            regsub -all "^-|:\[ \t\]" $L "" L
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
            # Reading the configfile column names (not the options)
            # - Users can select only a subset of the annotation columns provided by AnnotSV
            # - Essential annotations are added below (at the end of the proc)
            regsub "( |\t|\\*)+$" $L "" L
            lappend g_AnnotSV(outputColHeader) $L
        }
    }
    
    ##################################
    ## Load options given in arguments
    ##################################
    puts "\t...configuration data given in arguments"
    set i 0
    set j 1
    while {$j < [llength $argv]} {
        set optionName [lindex $argv $i]
        regsub "^-|:\[ \t\]" $optionName "" optionName
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
    puts "\t...checking all these configuration data\n"
    
    ## annotationsDir, variantconvertDir: It must be existing directories (or "" for the default)
    if {$g_AnnotSV(annotationsDir) eq ""} {
        set g_AnnotSV(annotationsDir) "$g_AnnotSV(installDir)/share/AnnotSV"
    } else {
        regsub "/+$" $g_AnnotSV(annotationsDir) "" g_AnnotSV(annotationsDir)
        if {![file isdirectory $g_AnnotSV(annotationsDir)]} {
            puts "AnnotSV needs in argument an existing path of the annotations directory (-annotationsDir = \"$g_AnnotSV(annotationsDir)\") - Exit with error."
            exit 2
        }
    }
    if {$g_AnnotSV(variantconvertDir) eq ""} {
        set g_AnnotSV(variantconvertDir) "$g_AnnotSV(installDir)/share/python3/variantconvert/"
    } else {
        regsub "/+$" $g_AnnotSV(variantconvertDir) "" g_AnnotSV(variantconvertDir)
        if {![file isdirectory $g_AnnotSV(variantconvertDir)]} {
            puts "AnnotSV needs in argument an existing path of the variantconvert directory (-variantconvertDir = \"$g_AnnotSV(variantconvertDir)\") - Exit with error."
            exit 2
        }
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
    
    ## benignAF: [0.001-0.1], default = 0.01
    if {![string is double $g_AnnotSV(benignAF)] || $g_AnnotSV(benignAF) < 0.001 || $g_AnnotSV(benignAF) > 0.1} {
        puts "############################################################################"
        puts "Bad option value: -benignAF = $g_AnnotSV(benignAF)"
        puts "Should be in the \[0.001-0.1\] range values, default = 0.01"
        puts "############################################################################"
        exit 2
    }
    
    ## metrics: us or fr
    if {$g_AnnotSV(metrics) ne "us" && $g_AnnotSV(metrics) ne "fr"} {
        puts "############################################################################"
        puts "Bad option value: -metrics = $g_AnnotSV(metrics)"
        puts "Should be \"us\" or \"fr\" - Continue with the \"us\" default."
        puts "############################################################################"
        set g_AnnotSV(metrics) "us"
    }
    
    ## It must be a boolean: 1 or 0
    foreach val {candidateGenesFiltering includeCI overwrite reciprocal REreport REselect1 REselect2 SVinputInfo snvIndelPASS vcf} {
        if {$g_AnnotSV($val) ne "1" && $g_AnnotSV($val) ne "0"} {
            puts "############################################################################"
            puts "Bad option value: -$val = $g_AnnotSV($val)"
            puts "Should be \"1\" or \"0\""
            puts "############################################################################"
            exit 2
        }
    }
    
    # g_AnnotSV(hpo)
    # It must be a comma, semicolon or space separated class values, default = ""
    # (e.g.: "HP:0001156,HP:0001363,HP:0011304")
    set g_AnnotSV(hpo) [split $g_AnnotSV(hpo) ",|;| "]
    set L_correctHPO {}
    set L_display ""
    foreach hpo $g_AnnotSV(hpo) {
        if {[regexp "^HP:\[0-9\]+$" $hpo]} {
            lappend L_correctHPO $hpo
        } else {
            lappend L_display "Bad format for the HPO term: $hpo. Not used."
        }
    }
    if {$L_display ne ""} {
        puts "############################################################################"
        puts "[join $L_display "\n"]"
        puts "############################################################################"
    }
    set g_AnnotSV(hpo) [join $L_correctHPO ","]
    
    # g_AnnotSV(rankFiltering)
    ## It must be a list between 1 and 5 and NA
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
        if {[lsearch -exact {1 2 3 4 5 NA} $i] ne -1} {lappend liste $i}
    }
    set liste [lsort -unique $liste]
    set g_AnnotSV(rankFiltering) $liste
    
    ## g_AnnotSV(outputDir) must be an existing directory. Else, we create it.
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
    regsub "/+$" $g_AnnotSV(outputDir) "" g_AnnotSV(outputDir)
    if {$g_AnnotSV(outputFile) ne ""} {
        set g_AnnotSV(outputFile) [file tail $g_AnnotSV(outputFile)]
        if {[string range $g_AnnotSV(outputFile) end-3 end] ne ".tsv"} {
            set g_AnnotSV(outputFile) "$g_AnnotSV(outputFile).tsv"
        }
    } else {
        regsub -nocase "(\\.bed|\\.vcf(\\.gz)?)$" [file tail $g_AnnotSV(SVinputFile)] ".annotated.tsv" g_AnnotSV(outputFile)
    }
    if {!$g_AnnotSV(overwrite) && [file exists $g_AnnotSV(outputDir)/$g_AnnotSV(outputFile)]} {
        puts "############################################################################"
        puts "Bad value for the -outputFile option, file already exists ($g_AnnotSV(outputDir)/$g_AnnotSV(outputFile)) - Exit with error."
        puts "############################################################################"
        exit 2
    }
    set g_AnnotSV(outputDir) [file normalize $g_AnnotSV(outputDir)]
 
    # bedtools/bcftools: It should be a good path that we can run
    foreach tool {bedtools bcftools} {
        if {[catch {eval exec $g_AnnotSV($tool) --version} Message]} {
            puts "############################################################################"
            puts "Bad value for the $tool option ($g_AnnotSV($tool))"
            puts "$Message"
            puts "Exit with error."
            puts "############################################################################"
            exit 2
        } else {
			# Minimum bedtools version compatible with AnnotSV is version 2.25
			if {$tool eq "bedtools"} {
				if {[regexp "(\[0-9\]+\\.\[0-9\]+\\.\[0-9\]+)" $Message match bedtoolsVersion]} {
					set L_bedtoolsVersion [split $bedtoolsVersion "."]
					if {[lindex $L_bedtoolsVersion 0] < 2 || ([lindex $L_bedtoolsVersion 0] eq 2 && [lindex $L_bedtoolsVersion 1] < 25)} {
					    puts "############################################################################"
					    puts "The minimum bedtools version compatible with AnnotSV is version 2.25."
						puts "You are using bedtools version $bedtoolsVersion"
					    puts "Exit with error."
					    puts "############################################################################"
						exit 2
					}
				}
			}
		}
    }
    
    # tx: It must be "RefSeq" or "ENSEMBL"
    if {$g_AnnotSV(tx) ne "RefSeq" && $g_AnnotSV(tx) ne "ENSEMBL"} {
        puts "############################################################################"
        puts "Bad option value: -tx = $g_AnnotSV(tx)"
        puts "Should be \"RefSeq\" or \"ENSEMBL\""
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
    
    ## candidateGenesFile, externalGeneFiles: It must be an existing file
    foreach annfile {candidateGenesFile externalGeneFiles} {
        if {$g_AnnotSV($annfile) ne "" && ![file exists $g_AnnotSV($annfile)]} {
            puts "############################################################################"
            puts "Bad value for the $annfile option, file does not exist ([set g_AnnotSV($annfile)]) - Exit with error."
            puts "############################################################################"
            exit 2
        }
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
    set L_annotationMode {both full split}
    if {[lsearch -exact $L_annotationMode "$g_AnnotSV(annotationMode)"] eq -1} {
        puts "############################################################################"
        puts "Bad option value: -annotationMode = $g_AnnotSV(annotationMode)"
        puts "Should be \"both\", \"full\" or \"split\""
        puts "############################################################################"
        exit 2
    }
    
    ## The following step could be improved: too long
    #################################################
    set g_AnnotSV(snvIndelSamples) [split $g_AnnotSV(snvIndelSamples)  ";|,"]
    ## If "snvIndelFiles" option is defined:
    ## It must be existing files
    if {$g_AnnotSV(snvIndelFiles) ne ""} {
        # Warning: A file that doesn't exists (or badly defined with a regexp) will not return an error message but will return an empty list:
        set g_AnnotSV(snvIndelFiles) [eval glob -nocomplain $g_AnnotSV(snvIndelFiles)]
    }
    ## After that, if "snvIndelFiles" contains at least 1 existing file:
    if {$g_AnnotSV(snvIndelFiles) ne ""} {
        ## All "$snvIndelSamples" should be presents in one of the VCF files
        ## If the "-snvIndelSamples" option is not defined, we defined it with all samples from the "snvIndelFiles".
        set L_correctSamples {}
        set L_allSamplesFromVCF {}
        foreach vcfF $g_AnnotSV(snvIndelFiles) {
            if {[regexp -nocase ".gz$" $vcfF]} {
                set f [open "|gzip -cd $vcfF"]
            } else {
                set f [open "$vcfF"]
            }
            while {![eof $f]} {
                set L [gets $f]
                if {[string range $L 0 5] ne "#CHROM"} {continue}
                set Ls [split $L "\t"]
                lappend L_allSamplesFromVCF {*}[lrange $Ls 9 end]
                foreach sample $g_AnnotSV(snvIndelSamples) {
                    if {[lsearch -exact $Ls $sample] ne -1} {lappend L_correctSamples $sample}
                }
                # We can't put a break here, because the pipe command: "|gzip ..." still produces data.
                # break will creates a broken pipe signal:
                #   child killed: write on pipe with no readers
                #   while executing "close $f"
            }
            close $f
        }
        if {$g_AnnotSV(snvIndelSamples) eq ""} {
            set L_allSamplesFromVCF [lsort -unique $L_allSamplesFromVCF]
            set g_AnnotSV(snvIndelSamples) $L_allSamplesFromVCF
        } else {
            if {$L_correctSamples ne ""} {
                set L_correctSamples [lsort -unique $L_correctSamples]
                set g_AnnotSV(snvIndelSamples) $L_correctSamples; # Remove samples from g_AnnotSV(snvIndelSamples) that are not in a VCF file
            } else {
                set L_allSamplesFromVCF [lsort -unique $L_allSamplesFromVCF]
                set g_AnnotSV(snvIndelSamples) $L_allSamplesFromVCF
            }
        }
    } else {
        set g_AnnotSV(snvIndelSamples) ""
    }
    
    
    set g_AnnotSV(candidateSnvIndelSamples) [split $g_AnnotSV(candidateSnvIndelSamples)  ";|,"]
    ## If "candidateSnvIndelFiles" option is defined:
    ## It must be existing files
    if {$g_AnnotSV(candidateSnvIndelFiles) ne ""} {
        # Warning: A file that doesn't exist (or badly defined with a regexp) will not return an error message but will return an empty list:
        set g_AnnotSV(candidateSnvIndelFiles) [eval glob -nocomplain $g_AnnotSV(candidateSnvIndelFiles)]
    }
    ## After that, if "candidateSnvIndelFiles" contains at least 1 existing file:
    if {$g_AnnotSV(candidateSnvIndelFiles) ne ""} {
        ## All "$candidateSnvIndelSamples" should be presents in one of the VCF files
        ## If the "-candidateSnvIndelSamples" option is not defined, we defined it with all samples from the "candidateSnvIndelFiles".
        set L_correctCandidateSamples ""
        set L_allSamplesFromVCF {}
        foreach vcfF $g_AnnotSV(candidateSnvIndelFiles) {
            if {[regexp -nocase ".gz$" $vcfF]} {
                set f [open "|gzip -cd $vcfF"]
            } else {
                set f [open "$vcfF"]
            }
            while {![eof $f]} {
                set L [gets $f]
                if {[string range $L 0 5] eq "#CHROM"} {
                    set Ls [split $L "\t"]
                    lappend L_allSamplesFromVCF {*}[lrange $Ls 9 end]
                    foreach sample $g_AnnotSV(candidateSnvIndelSamples) {
                        if {[lsearch -exact $Ls $sample] ne -1} {lappend L_correctCandidateSamples $sample}
                    }
                    # We can't put a break here, because the pipe command: "|gzip ..." still produces data.
                    # break will creates a broken pipe signal:
                    #   child killed: write on pipe with no readers
                    #   while executing "close $f"
                }
            }
            close $f
        }
        
        if {$g_AnnotSV(candidateSnvIndelSamples) eq ""} {
            set L_allSamplesFromVCF [lsort -unique $L_allSamplesFromVCF]
            set g_AnnotSV(candidateSnvIndelSamples) $L_allSamplesFromVCF
        } else {
            if {$L_correctCandidateSamples ne ""} {
                set L_correctCandidateSamples [lsort -unique $L_correctCandidateSamples]
                set g_AnnotSV(candidateSnvIndelSamples) $L_correctCandidateSamples; # Remove samples from g_AnnotSV(candidateSnvIndelSamples) that are not in a VCF file
            } else {
                set L_allSamplesFromVCF [lsort -unique $L_allSamplesFromVCF]
                set g_AnnotSV(candidateSnvIndelSamples) $L_allSamplesFromVCF
            }
        }
        
        # Check if the "candidateSnvIndelFiles" contains the GT field
        # (else the g_AnnotSV(candidateSnvIndelSamples) is set to "")
        filteredVCFannotation "FULL" "" "" "" ""
        
    } else {
        set g_AnnotSV(candidateSnvIndelSamples) ""
    }
    
    ## It must be "GRCh37" or "GRCh38" or "mm39" or "mm9" or "mm10" for the genomeBuild option.
    if {![regexp -nocase "^(GRCh37)|(GRCh38)|(mm39)|(mm9)|(mm10)$" $g_AnnotSV(genomeBuild)]} {
        puts "############################################################################"
        puts "Bad option value: -genomeBuild = $g_AnnotSV(genomeBuild)"
        puts "Should be \"GRCh37\", \"GRCh38\", \"mm39\", \"mm9\" or \"mm10\""
        puts "############################################################################"
        exit 2
    }
    # Definition of the $g_AnnotSV(organism) variable
    if {[regexp -nocase "^(GRCh37)|(GRCh38)$" $g_AnnotSV(genomeBuild)]} {
        set g_AnnotSV(organism) "Human"
    } elseif {[regexp -nocase "^(mm39)|(mm9)|(mm10)$" $g_AnnotSV(genomeBuild)]} {
        set g_AnnotSV(organism) "Mouse"
    }
    
    
    # Some annotation columns can not be removed by the user:
    # - Annotations which identify the SV:                     AnnotSV_ID SV_chrom SV_start SV_end SV_length SV_type
    # - Annotations essential for the ranking:                 Annotation_mode Gene_name Gene_Count Overlapped_CDS_percent Frameshift Location Location2 RE_gene Overlapped_CDS_length Exon_count
    #                                                          P_gain_coord P_loss_coord P_snvindel_nb B_gain_coord B_loss_coord
    #                                                          po_P_gain_coord po_P_loss_coord po_B_gain_allG_coord po_B_gain_someG_coord po_B_loss_allG_coord po_B_loss_someG_coord
    #                                                          HI TS GnomAD_pLI LOEUF_bin DDD_HI_percent PhenoGenius_specificity Exomiser_gene_pheno_score OMIM_morbid
    # - Annotations essential for variantconvert (VCF output): Samples_ID
    # - Ranking annotations:                                   AnnotSV_ranking_score AnnotSV_ranking_criteria ACMG_class
    # - Annotations linked to the previous essential annotations: Tx Tx_start Tx_end Overlapped_tx_length Dist_nearest_SS Nearest_SS_type Intersect_start Intersect_end
    #   (to be improved with the python reimplementation)         P_gain_phen P_gain_hpo P_gain_source P_loss_phen P_loss_hpo P_loss_source
    #                                                             po_P_gain_phen po_P_gain_hpo po_P_gain_source po_P_gain_percent po_P_loss_phen po_P_loss_hpo po_P_loss_source po_P_loss_percent
    #                                                             P_snvindel_phen B_gain_source B_gain_AFmax B_loss_source B_loss_AFmax
    #                                                             po_B_gain_allG_source po_B_gain_someG_source po_B_gain_someG_coord po_B_loss_allG_source po_B_loss_someG_source
    set g_AnnotSV(genesBasedAnn) 1
    foreach col "SV_chrom SV_start SV_end SV_length SV_type Annotation_mode Gene_name Gene_Count RE_gene P_gain_coord P_loss_coord P_snvindel_nb B_gain_coord B_loss_coord po_P_gain_coord po_P_loss_coord po_B_gain_allG_coord po_B_loss_allG_coord po_B_loss_someG_coord HI TS GnomAD_pLI LOEUF_bin DDD_HI_percent PhenoGenius_specificity Exomiser_gene_pheno_score OMIM_morbid Location Location2 Overlapped_CDS_percent Frameshift Exon_count Samples_ID" {
        if {[lsearch -exact "$g_AnnotSV(outputColHeader)" $col] eq -1} {
            lappend g_AnnotSV(outputColHeader) $col
        }
    }
    
    return
}

proc createBEDinputHeaderFile {} {
    
    global g_AnnotSV
    global headerFileToRemove
    global addNAinBED
    
    set addNAinBED 0
    
    ## SVinputfile is a BED
    regsub -nocase ".bed$" $g_AnnotSV(bedFile) ".header.tsv" BEDinputHeaderFile
    set BEDinputHeaderFile "$g_AnnotSV(outputDir)/[file tail $BEDinputHeaderFile]"
    
    ## Creation of the header file
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
        if {$g_AnnotSV(samplesidBEDcol) == -1} {
            set i_samplesid [lsearch -exact "Samples_ID" [split $header "\t"]]
            if {$i_samplesid ne -1} {
                set g_AnnotSV(samplesidBEDcol) [expr {$i_samplesid+1}]
            } else {
                append header "\tSamples_ID"
                set addNAinBED 1
                set g_AnnotSV(samplesidBEDcol) [llength [split $header "\t"]]
            }
        }
        WriteTextInFile "$header" $BEDinputHeaderFile
        set headerFileToRemove 1
    } else {
        set header "SV_chrom\tSV_start\tSV_end"
        set theBEDlength [llength [split [FirstLineFromFile $g_AnnotSV(bedFile)] "\t"]]
        set i 4
        set k 1
        while {$i <= $theBEDlength} {
            if {$i eq $g_AnnotSV(svtBEDcol)} {
                append header "\tSV_type"
            } elseif {$i eq $g_AnnotSV(samplesidBEDcol)} {
                append header "\tSamples_ID"
            } else {
                append header "\tuser#$k"; incr k
            }
            incr i
        }
        
        if {$g_AnnotSV(samplesidBEDcol) == -1} {
            append header "\tSamples_ID"
            set addNAinBED 1
            set g_AnnotSV(samplesidBEDcol) [expr {$theBEDlength+1}]
        }
        WriteTextInFile "$header" $BEDinputHeaderFile
        set headerFileToRemove 1
    }
    return
}

proc addNAinSamplesIDbedCol {} {
    
    global g_AnnotSV
    global addNAinBED
    
    ## SVinputfile is a BED, with no samplesid column
    if {$addNAinBED} {
        # Add a supplementary column in the BED input file with "NA" (corresponding to the sample name)
        # => needed for variantconvert + to create the SV database (in the future)
        set L_toWrite {}
        set i 0
        set f [open $g_AnnotSV(bedFile)]
        regsub -nocase ".bed$" $g_AnnotSV(bedFile) ".NA.bed" g_AnnotSV(NAbedFile)
        file delete -force "$g_AnnotSV(NAbedFile)"
        while {![eof $f]} {
            set L [gets $f]
            if {$L eq ""} {continue}
            if {[regexp "^#" $L]} {continue}
            set L "$L\tNA"
            #if {![info exists theNAbedLength]} {set theNAbedLength [llength [split $L "\t"]]}
            lappend L_toWrite $L
            incr i
            if {$i > 100000} {
                WriteTextInFile [join "$L_toWrite" "\n"] $g_AnnotSV(NAbedFile)
                set L_toWrite {}
                set i 0
            }
        }
        close $f
        WriteTextInFile [join "$L_toWrite" "\n"] $g_AnnotSV(NAbedFile)
        
        regsub -nocase ".bed$" $g_AnnotSV(bedFile) ".header.tsv" BEDinputHeaderFile
        set BEDinputHeaderFile "$g_AnnotSV(outputDir)/[file tail $BEDinputHeaderFile]"
        regsub -nocase ".bed$" $g_AnnotSV(bedFile) ".NA.header.tsv" NAbedinputHeaderFile
        set NAbedinputHeaderFile "$g_AnnotSV(outputDir)/[file tail $NAbedinputHeaderFile]"
        file rename -force $BEDinputHeaderFile $NAbedinputHeaderFile
        
        set g_AnnotSV(bedFile) $g_AnnotSV(NAbedFile)
        
        # Number of the "Samples_ID" column in the new "NA" BED file (not the informatic count in a list!)
        # set g_AnnotSV(samplesidBEDcol) "$theNAbedLength"
        
    }
}

