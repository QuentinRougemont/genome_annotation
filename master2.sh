#!/bin/bash

#master file to run the pipeline :
#Date: 11-2023
#Author: QR

########################################
## Global variables
########################################
echo "source config file..."
source ./config/config 

echo "------------------------------------------------------------"
echo "-----check all variables from the configuration file  ------"
echo "*** ancestral genome is ${anestral_genome} ***"
echo "*** haplotype1 is ${haplotype1} ***"
echo "*** haplotype2 is ${haplotype2}"
echo "*** RNAseq? $rnaseq ***"
echo "*** gtf ? $gtf ***" 
echo "*** TE database is : $TEdatabase ***"
echo "*** lineage for busco analyses will be:  $busco_lineage ***"
echo -e "------------------------------------------------------------\n"


# if no TE database then exit

# if no lineage for busco then exit and ask for it 

############################################################
# Process the input options.                               #
############################################################
#a closely related lineage for busco inference is also needed!
#for TE: a database of TE is needed!
#to do: revoir le script de TE pour faire juste de-novo + one database

if [ -z "${haplotype1}" ]  || [ -z "${haplotype2}" ] ; then #|| [ -z "${folderpath}" ]  || [ -z "${ancestral_sp}" ]    ; then
	echo "Error! provide the genome of at least one species"
	Help
	exit 2
elif [ -n "${haplotype1}" ] && [ -z "${haplotype2}" ]  && [ $rnaseq = "NO"]   ; then
	cd haplo1
	if ["$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
	echo "only the genome of one species or was provided" 
	echo "we will only perform TE detection and genome annotation"
	echo "genome is ${haplotype1} "
	../run_script_05_to_08.sh
	else
		echo "error no fasta file in 03_genome"
	        echo "please copy your genome here"	
		help
		exit 1
	fi

elif [ -z "${haplotype1}" ] && [ -n "${haplotype2}" ] && [ $rnaseq = "NO"]    ; then
	cd haplo2
	if ["$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
	echo "only the genome of one species or was provided" 
	echo "we will only perform TE detection and genome annotation"
	echo "genome is ${haplotype2} "
	../run_script_05_to_08.sh
	else
		echo "error no fasta file in 03_genome"
	        echo "please copy your genome here"	
		help
		exit 1
	fi

elif [ -n "${haplotype1}" ] && [ -z "${haplotype2}" ] && [ $rnaseq = "YES"]  ; then
	cd haplo1/
	if ["$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
	    echo "only the genome of one species or was provided" 
	    echo "we will only perform TE detection and genome annotation with RNAseq "
	    echo "genome is ${haplotype1} "
	    ../run_step_rnaseq.sh
	    #check that this script was sucessfull else kill:
	    if [ $? -eq 0 ]; then
    	        echo rnaseq mapping succesffull
 	        ../run_script_05_to_08.sh

	    else
    	        echo ERROR - verfiy braker outputs!   
    	        exit 1
	    fi
	else
		echo "error no fasta file in 03_genome"
	        echo "please copy your genome here"	
		help
		exit 1
	fi

	
	
elif [ -z "${haplotype1}" ] && [ -n "${haplotype2}" ]  && [ $rnaseq = "YES" ]   ; then
	cd haplo2/
	if ["$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
	    echo "only the genome of one species or was provided" 
	    echo "we will only perform TE detection and genome annotation with RNAseq "
	    echo "genome is ${haplotype2} "
	    ../run_step_rnaseq.sh
	    #check that this script was sucessfull else kill:
	    if [ $? -eq 0 ]; then
    	        echo rnaseq mapping succesffull
 	        ../run_script_05_to_08.sh

	    else
    	        echo ERROR - verfiy braker outputs!   
    	        exit 1
	    fi
	else
		echo "error no fasta file in 03_genome"
	        echo "please copy your genome here"	
		help
		exit 1
	fi


	#if both species are provided without RNAseq:
elif [ -n "${haplotype1}" ] && [ -n "${haplotype2}" ]  && [[ $rnaseq = "NO" ]]   ; then
	echo "we will perform TE detection - genome annotation - Ds computation and plots"
	echo "genome are $haplotype1 and $haplotype2 "
	cd haplo1
	if ["$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
	    ../run_script_05_to_08.sh
	    #verify that alll run correctly 

	else
		echo "error no fasta file in 03_genome"
	        echo "please copy your genome here"	
		help
		exit 1
	fi

	cd haplo2
	if ["$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
	     ../run_script_05_to_08.sh
	    #verify that alll run correctly 

	else
		echo "error no fasta file in 03_genome"
	        echo "please copy your genome here"	
		help
		exit 1
	fi

	#if all is OK then run GeneSpace - paml etc :
	#modifiy the script RunGeneSpace etc to handle case with/without ancestral species


elif [ -n "${haplotype1}" ] && [ -n "${haplotype2}" ] && [[ $rnaseq = "YES" ]]    ; then
	echo "we will perform TE detection - genome annotation with RNAseq - Ds computation and plots"
	echo "genomes are $haplotype1 and $haplotype2"
	cd haplo1/
	if ["$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
	    ../run_step_rnaseq.sh
	    #check that this script was sucessfull else kill:
	    if [ $? -eq 0 ]; then
    	        echo rnaseq mapping succesffull
 	        ../run_script_05_to_08.sh

	    else
    	        echo ERROR - verfiy braker outputs!   
    	        exit 1
	    fi
	else
		echo "error no fasta file in 03_genome"
	        echo "please copy your genome here"	
		help
		exit 1
	fi

	cd haplo2/
	if ["$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
	    ../run_step_rnaseq.sh
	    #check that this script was sucessfull else kill:
	    if [ $? -eq 0 ]; then
    	        echo rnaseq mapping succesffull
 	        ../run_script_05_to_08.sh

	    else
    	        echo ERROR - verfiy braker outputs!   
    	        exit 1
	    fi
	else
		echo "error no fasta file in 03_genome"
	        echo "please copy your genome here"	
		help
		exit 1
	fi

	#then run GeneSpace etc :
	#modifiy the script RunGeneSpace etc to handle case with/without ancestral species

elif [ -n "${haplotype1}" ] && [ -n "${haplotype2}" ] && [[ $rnaseq = "NO" ]] && [ -n "$ancestral_sp" ]  ; then
	echo "we will perform all analyses with annotations performed without rnaseq "
	echo "genomes are $haplotype1 and $haplotype2"

	cd haplo1
	if ["$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
	    ../run_script_05_to_08.sh
	    #verify that alll run correctly 

	else
		echo "error no fasta file in 03_genome"
	        echo "please copy your genome here"	
		help
		exit 1
	fi

	cd haplo2
	if ["$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
	     ../run_script_05_to_08.sh
	    #verify that alll run correctly 

	else
		echo "error no fasta file in 03_genome"
	        echo "please copy your genome here"	
		help
		exit 1
	fi

	#then run GeneSpace etc :
	#modifiy the script RunGeneSpace etc to handle case with/without ancestral species

elif [ -n "${haplotype1}" ] && [ -n "${haplotype2}" ]  && [[ $rnaseq = "YES" ]]  && [ -n "$ancestral_sp" ]   ; then
	echo "we will perform all analyses including annotation with rnaseq"
	echo "genomes are $haplotype1 and $haplotype2"
	cd haplo1/
	if ["$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
	    ../run_step_rnaseq.sh
	    #check that this script was sucessfull else kill:
	    if [ $? -eq 0 ]; then
    	        echo rnaseq mapping succesffull
 	        ../run_script_05_to_08.sh

	    else
    	        echo ERROR - verfiy braker outputs!   
    	        exit 1
	    fi
	else
		echo "error no fasta file in 03_genome"
	        echo "please copy your genome here"	
		help
		exit 1
	fi

	cd haplo2/
	if ["$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
	    ../run_step_rnaseq.sh
	    #check that this script was sucessfull else kill:
	    if [ $? -eq 0 ]; then
    	        echo rnaseq mapping succesffull
 	        ../run_script_05_to_08.sh

	    else
    	        echo ERROR - verfiy braker outputs!   
    	        exit 1
	    fi
	else
		echo "error no fasta file in 03_genome"
	        echo "please copy your genome here"	
		help
		exit 1
	fi

	
	#then run GeneSpace etc :
	#modifiy the script RunGeneSpace etc to handle case with/without ancestral species

fi

#here handle cases where already annotated genome with gtf are provided
