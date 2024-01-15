############################################################################################################
# AnnotSV 3.3.9                                                                                            #
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


# - Check if the reference fasta file in the variantconvert bed configfiles has the good path
# - Check if the "pip install -e ." command was already run
proc checkVariantconvertConfigfile {} {

    global g_AnnotSV

    # Use of variantconvert to create a vcf output
    if {$g_AnnotSV(vcf)} {
	
	# - Check if the reference fasta file in the variantconvert bed configfiles has the good path
        #############################################################################################

	## SVinputfile is a BED
	## (useful only with a BED input file, because there is no need to have a reference fasta file from a VCF SVinputfile)
	if {[regexp -nocase "\\.bed$" $g_AnnotSV(SVinputFile)]} {

            set configfile "$g_AnnotSV(variantconvertDir)/configs/$g_AnnotSV(genomeBuild)/annotsv3_from_bed.json"
            set localconfigfile "$g_AnnotSV(variantconvertDir)/configs/$g_AnnotSV(genomeBuild)/annotsv3_from_bed.local.json"

	    if {$g_AnnotSV(annotationsDir) ne ""} {
		set pathDir "$g_AnnotSV(annotationsDir)"
	    } else {
		set pathDir "$g_AnnotSV(installDir)/share/AnnotSV/"
	    }

	    set distributedPathLine "\"path\": \".*\","
	    set newPathLine         "\"path\": \"$pathDir/Annotations_Human/BreakpointsAnnotations/GCcontent/$g_AnnotSV(genomeBuild)/$g_AnnotSV(genomeBuild)_chromFa.fasta\","
	    set distributedRefLine  "\"##reference=file:.*\""
	    set newRefLine          "\"##reference=file:$pathDir/Annotations_Human/BreakpointsAnnotations/GCcontent/$g_AnnotSV(genomeBuild)/$g_AnnotSV(genomeBuild)_chromFa.fasta\""

	    # 1 - AnnotSV install with the root user
	    # 2 - AnnotSV run with non-root user
	    # => The $localconfigfile can not be created by a non-root user. This file should exists with 777 permissions
	    if {![file exists $localconfigfile] || [file size $localconfigfile] eq 0} {
		set L_Lines {}
		foreach L [LinesFromFile $configfile] {
		    if {[regexp "$distributedPathLine" $L]} {
			regsub "$distributedPathLine" $L "$newPathLine" L
		    }
		    if {[regexp "$distributedRefLine" $L]} {
			regsub "$distributedRefLine" $L "$newRefLine" L
		    }
		    lappend L_Lines $L
		}
		ReplaceTextInFile [join $L_Lines "\n"] $localconfigfile
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

    catch {exec python3 $g_AnnotSV(variantconvertDir)/variantconvert --version} Message
    if {[regexp "variantconvert (\[0-9\]+\\.\[0-9\]+\\.\[0-9\]+)" $Message match version]} {
        set version "v$version "
    } else {
        set version ""
    }

    puts "...creation of the VCF output file: $VCFoutputFile"
    puts "   AnnotSV relies on the variantconvert tool ${version}(https://github.com/SamuelNicaise/variantconvert)."
    puts "   A minimal Python 3.8 installation is required, as well as the natsort, panda and pyfaidx Python modules."

    regsub "\.vcf$" $VCFoutputFile ".variantconvert.log" LogFile

    if {[regexp -nocase "\\.vcf(.gz)?$" $g_AnnotSV(SVinputFile)]} {
        ## SVinputfile is a VCF
        set command "python3 $g_AnnotSV(variantconvertDir)/variantconvert convert -i $outputFile -o $VCFoutputFile -fi annotsv -fo vcf -c $g_AnnotSV(variantconvertDir)/configs/$g_AnnotSV(genomeBuild)/annotsv3_from_vcf.json"
    } else {
        ## SVinputfile is a BED)
        if {$g_AnnotSV(svtBEDcol) == -1 } {
            puts "   WARNING: With a \"BED\" SV input file, the user has to define the -svtBEDcol option."
            puts "            => could not create the VCF output file:"
            if {$g_AnnotSV(svtBEDcol) == -1}       {puts "               -svtBEDcol $g_AnnotSV(svtBEDcol)"}
            return
        } else {
            set command "python3 $g_AnnotSV(variantconvertDir)/variantconvert convert -i $outputFile -o $VCFoutputFile -fi annotsv -fo vcf -c $g_AnnotSV(variantconvertDir)/configs/$g_AnnotSV(genomeBuild)/annotsv3_from_bed.local.json"
        }
    }

    # variantconvert output

    catch {eval exec $command} Message
    if {[regexp -nocase "error" $Message]} {
        puts "Error:"
    }
    ReplaceTextInFile "$command\n\n$Message" $LogFile
    puts "   => cf $LogFile"

    return
}
