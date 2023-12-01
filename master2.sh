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
echo "*** fasta for genome1 is ${genome1} **** "
echo "*** fasta for genome2 is ${genome2} **** "
echo "*** haplotype1 is ${haplotype1} ***"
echo "*** haplotype2 is ${haplotype2}"
echo "*** RNAseq? $rnaseq ***"
echo "*** gtf ? $gtf ***" 
echo "*** TE database is : $TEdatabase ***"
echo "*** lineage for busco analyses will be:  $busco_lineage ***"
echo "*** annotate? $annotate" 
echo -e "------------------------------------------------------------\n"

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo " "
   echo "Usage: $0 [-h]"
   echo "-h|--help: Print this Help."
   echo " "
   echo "TO BE FILLED LATER"
}

###########################################################
#a few minore chek here:
###########################################################
# if no TE database then exit
# if no lineage for busco then exit and ask for it 
if [ -z "${TEdatabase}" ] && [ $annotate = YES ] ;
then
	echo "WARNING NO Database for TE is provided "
	echo "I will not be able to assess mask TE prior to the genome annotation step "
	echo "this is very bad" 
	exit
	Help
fi


# if no lineage for busco then exit and ask for it 
if [ -z "${busco_lineage}" ] ;
then
	echo "WARNING NO lineage provided for busco" 
	echo "I will not be able to assess the quality of the runs, which is compulsory for these anaylses"
	exit
	Help
fi


############################################################
# Generate architecture:
############################################################

mkdir -p haplo1/03_genome 

# ----- check compression of fasta  ------ ##
#check compression
if file --mime-type "$genome1" | grep -q gzip$; then
   echo "$genome1 is gzipped"
   gunzip "$genome1"
   genome1=${genome1%.gz}
   cd haplo1/03_genome 
   cp $genome1 $haplo1.fa
else
   echo "$genome1 is not gzipped"
   genome1=$genome1
   cd haplo1/03_genome 
   cp $genome1 $haplo1.fa
fi

if [ -z "$haplotype2"] && [ -z "$genome2" ]  ; then mkdir -p haplo2/03_genome ; fi

#check genome compression:
# ----- check compression of fasta  ------ ##
#check compression
if file --mime-type "$genome2" | grep -q gzip$; then
   echo "$genome2 is gzipped"
   gunzip "$genome2"
   genome2=${genome%.gz}
   cd haplo2/03_genome 
   cp $genome2 $haplo2.fa

else
   echo "$genome2 is not gzipped"
   cd haplo2/03_genome 
   cp $genome2 $haplo2.fa
   genome2=$genome
fi

if [ -z "$ancestral_genome"] ; then mkdir ancestralsp ; fi

#artchitecture should be OK to proceed

############################################################
# Process the input options.                               #
############################################################

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
