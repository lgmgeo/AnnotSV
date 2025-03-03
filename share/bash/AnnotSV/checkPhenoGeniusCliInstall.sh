#!/bin/bash


# Bash script: 
# A successfully executed code should exit with code 0. 
# Other values indicate an error. 


# Poetry is required
####################
poetryDir=`which poetry` 

# The "$poetryDir" directory should exist
if [ ! -e "$poetryDir" ]
then
    echo "Poetry seems not to be installed"
    echo "WARNING: No PhenoGenius annotations available."
	exit 1
fi


# Poetry version
################
# poetryVersion=`poetry --version | sed "s/Poetry (version //" | sed "s/)//"`
# e.g. Poetry (version 1.5.1)
# echo "Poetry: $poetryVersion"


# PhenogeniusCli install (if needed)
####################################
mkdir -p $ANNOTSV/share/python3/phenogeniuscli
cd $ANNOTSV/share/python3/phenogeniuscli
if [ ! -d PhenoGeniusCli ]
then
	git clone https://github.com/kyauy/PhenoGeniusCli.git >> PhenoGeniusCli.install.log 2>&1
	cd ./PhenoGeniusCli >> PhenoGeniusCli.install.log 2>&1
	# Search for the last PhenoGenius version tested/validated with AnnotSV
	git checkout tags/v.1.1.3 >> ../PhenoGeniusCli.install.log 2>&1
	rm -rf .git .gitattributes .github .gitignore >> ../PhenoGeniusCli.install.log 2>&1
	poetry config keyring.enabled false
    poetry install &> ../poetry_install.log
fi

cd $ANNOTSV/share/python3/phenogeniuscli/PhenoGeniusCli

poetry run python3 phenogenius_cli.py --help  &> ../PhenoGenius.run.test1.log
if [ `grep -c -- "--hpo_list" ../PhenoGenius.run.test1.log` == "1" ]
then
	exit 0
else
	poetry config keyring.enabled false
	poetry install &> ../poetry_install.log ;# Can be run from different servers
    poetry run python3 phenogenius_cli.py --help &> ../PhenoGenius.run.test2.log
	if [ `grep -c -- "--hpo_list" ../PhenoGenius.run.test2.log` == 1 ]
	then
		exit 0
	else
		echo "phenogenius_cli.py not installed"	
		echo "WARNING: No PhenoGenius annotations available."
		exit 1
	fi
fi


