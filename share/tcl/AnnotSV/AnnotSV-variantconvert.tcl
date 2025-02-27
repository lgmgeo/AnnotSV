############################################################################################################
# AnnotSV 3.4.5                                                                                            #
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


# - Check if the reference fasta file in the variantconvert bed configfiles has the good path
# - Check if the "pip install -e ." command was already run
proc checkVariantconvertConfigfile {} {
    
    global g_AnnotSV
    
    # Use of variantconvert to create a vcf output
    if {$g_AnnotSV(vcf)} {
        
        # - Check if the reference fasta file in the variantconvert bed configfiles has the good path
        #############################################################################################
        
        foreach formatFile {bed vcf} {
            
            set configfile "$g_AnnotSV(variantconvertDir)/src/variantconvert/configs/$g_AnnotSV(genomeBuild)/annotsv3_from_$formatFile.json"

            if {$g_AnnotSV(variantconvertMode) eq "combined"} {
				set localConfigfile "$g_AnnotSV(variantconvertDir)/src/variantconvert/configs/$g_AnnotSV(genomeBuild)/annotsv3_from_$formatFile.combined.local.json"
		    } elseif {$g_AnnotSV(variantconvertMode) eq "full"} {
			    set localConfigfile "$g_AnnotSV(variantconvertDir)/src/variantconvert/configs/$g_AnnotSV(genomeBuild)/annotsv3_from_$formatFile.full.local.json"
		    } elseif {$g_AnnotSV(variantconvertMode) eq "fullsplit"} {
			    set localConfigfile "$g_AnnotSV(variantconvertDir)/src/variantconvert/configs/$g_AnnotSV(genomeBuild)/annotsv3_from_$formatFile.fullsplit.local.json"
            }
            
            if {$g_AnnotSV(annotationsDir) ne ""} {
                set pathDir "$g_AnnotSV(annotationsDir)"
            } else {
                set pathDir "$g_AnnotSV(installDir)/share/AnnotSV/"
            }
            
            set distributedPathLine "\"path\": \".*\","
            set newPathLine         "\"path\": \"$pathDir/Annotations_Human/BreakpointsAnnotations/GCcontent/$g_AnnotSV(genomeBuild)/$g_AnnotSV(genomeBuild)_chromFa.fasta\","

            set distributedRefLine  "\"##reference=file:.*\""
            set newRefLine          "\"##reference=file:$pathDir/Annotations_Human/BreakpointsAnnotations/GCcontent/$g_AnnotSV(genomeBuild)/$g_AnnotSV(genomeBuild)_chromFa.fasta\""

            set distributedModeLine  "\"mode\": \"combined\","
            set newfullModeLine      "\"mode\": \"full\","
            set newfullsplitModeLine "\"mode\": \"full\\&split\","

            # 1 - AnnotSV install with the root user
            # 2 - AnnotSV run with non-root user
            # => The $localConfigfile can not be created by a non-root user. This file should exists with 777 permissions

			# If the user defined a wrong "-annotationsDir" during the first execution of AnnotSV (using the -vcf 1 parameter), the $localConfigfile need to be removed then recomputed.
            if {[file exists $localConfigfile]} {
				set testPathExists 0
				set testRefExists 0
                foreach L [LinesFromFile $localConfigfile] {
					if {[regexp "\"path\": \"(.*)\"," $L match path]} {
						if {[file exists $path]} {set testPathExists 1}
					}
                    if {[regexp "\"##reference=file:(.*)\"" $L match ref]} {
                        if {[file exists $ref]} {set testRefExists 1}
                    }
				}	
				if {$testPathExists eq 0 || $testRefExists eq 0} {
					file delete -force $localConfigfile
				}
			}
			# If the file has been deleted or does not yet exist: 
            if {![file exists $localConfigfile]} {
                set L_Lines {}
                foreach L [LinesFromFile $configfile] {
                    if {[regexp "$distributedPathLine" $L]} {
                        regsub "$distributedPathLine" $L "$newPathLine" L
                    }
                    if {[regexp "$distributedRefLine" $L]} {
                        regsub "$distributedRefLine" $L "$newRefLine" L
                    }
                    if {[regexp "$distributedModeLine" $L]} {
						if {$g_AnnotSV(variantconvertMode) eq "full"} {
	                        regsub "$distributedModeLine" $L "$newfullModeLine" L
						} elseif {$g_AnnotSV(variantconvertMode) eq "fullsplit"} {
                            regsub "$distributedModeLine" $L "$newfullsplitModeLine" L
						}
					}
                    lappend L_Lines $L 
				}
				ReplaceTextInFile [join $L_Lines "\n"] $localConfigfile
			}
        }
        
        # - Check if the "pip install -e ." command was already run
        ###########################################################
        
        # => Can be done during the installation with the Makefile (if the python environment is OK)
        set currentDir [pwd]
        catch {
            if {![file exists $g_AnnotSV(variantconvertDir)/pipinstall.flag]} {
                cd $g_AnnotSV(variantconvertDir)
                if {[catch {exec pip3 install -e .}]} {
                    if {[catch {exec pip install -e .} Message]} {
                        WriteTextInFile "$Message" $g_AnnotSV(variantconvertDir)/pipinstall.flag
                    } else {
                        WriteTextInFile "Done" $g_AnnotSV(variantconvertDir)/pipinstall.flag
                    }
                } else {
                    WriteTextInFile "Done" $g_AnnotSV(variantconvertDir)/pipinstall.flag
                }
            }
        }
        cd $currentDir
    }
}



proc runVariantconvert {outputFile} {
    
    global g_AnnotSV
    
    regsub "\.tsv$" $outputFile ".vcf" VCFoutputFile
    
    catch {exec python3 $g_AnnotSV(variantconvertDir)/src/variantconvert --version} Message
    if {[regexp "variantconvert (\[0-9\]+\\.\[0-9\]+\\.\[0-9\]+)" $Message match version]} {
        set version "v$version "
    } else {
        set version ""
    }
    
    puts "...creation of the VCF output file: $VCFoutputFile"
    puts "   AnnotSV relies on the variantconvert tool ${version}(https://github.com/SamuelNicaise/variantconvert)"
    
    regsub "\.vcf$" $VCFoutputFile ".variantconvert.log" LogFile
    

	# run variantconvert
	####################

    if {[regexp -nocase "\\.vcf(.gz)?$" $g_AnnotSV(SVinputFile)]} {

        ## SVinputfile is a VCF
		#######################

	    # configfile definition
	    if {$g_AnnotSV(variantconvertMode) eq "combined"} {
	        set localConfigfile "$g_AnnotSV(variantconvertDir)/src/variantconvert/configs/$g_AnnotSV(genomeBuild)/annotsv3_from_vcf.combined.local.json"
	    } elseif {$g_AnnotSV(variantconvertMode) eq "full"} {
	        set localConfigfile "$g_AnnotSV(variantconvertDir)/src/variantconvert/configs/$g_AnnotSV(genomeBuild)/annotsv3_from_vcf.full.local.json"
	    } elseif {$g_AnnotSV(variantconvertMode) eq "fullsplit"} {
	        set localConfigfile "$g_AnnotSV(variantconvertDir)/src/variantconvert/configs/$g_AnnotSV(genomeBuild)/annotsv3_from_vcf.fullsplit.local.json"
	    }

        set command "python3 $g_AnnotSV(variantconvertDir)/src/variantconvert convert -i $outputFile -o $VCFoutputFile -c $localConfigfile"

    } else {

        ## SVinputfile is a BED)
		########################

        # configfile definition
        if {$g_AnnotSV(variantconvertMode) eq "combined"} {
            set localConfigfile "$g_AnnotSV(variantconvertDir)/src/variantconvert/configs/$g_AnnotSV(genomeBuild)/annotsv3_from_bed.combined.local.json"
        } elseif {$g_AnnotSV(variantconvertMode) eq "full"} {
            set localConfigfile "$g_AnnotSV(variantconvertDir)/src/variantconvert/configs/$g_AnnotSV(genomeBuild)/annotsv3_from_bed.full.local.json"
        } elseif {$g_AnnotSV(variantconvertMode) eq "fullsplit"} {
            set localConfigfile "$g_AnnotSV(variantconvertDir)/src/variantconvert/configs/$g_AnnotSV(genomeBuild)/annotsv3_from_bed.fullsplit.local.json"
        }

        if {$g_AnnotSV(svtBEDcol) == -1 } {
            puts "   WARNING: With a \"BED\" SV input file, the user has to define the -svtBEDcol option."
            puts "            => could not create the VCF output file:"
            if {$g_AnnotSV(svtBEDcol) == -1}       {puts "               -svtBEDcol $g_AnnotSV(svtBEDcol)"}
            return
        } else {
            set command "python3 $g_AnnotSV(variantconvertDir)/src/variantconvert convert -i $outputFile -o $VCFoutputFile -c $localConfigfile"
        }
    }
    
    # variantconvert output
    #######################
    catch {eval exec $command} Message
	regsub -all "FutureWarning: Setting an item of incompatible dtype is deprecated and will raise in a future error" $Message "..." MessageReg  
    if {[regexp -nocase "error" $MessageReg]} {
	    puts "   (a minimal Python 3.8 installation is required, as well as the natsort, panda and pyfaidx Python modules)"
        puts "Error:"
		puts "   => cf $LogFile"
    }
    ReplaceTextInFile "$command\n\n$Message" $LogFile

    return
}



