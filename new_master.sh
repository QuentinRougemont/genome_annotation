#!/bin/bash

#master file to run the pipeline :
#Date: 11-2023
#Author: QR

source config/colors 
############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo " "
   echo "Usage: $0 [-h]"
   echo "-h|--help: Print this Help."
   echo -e "Main script to launch all subsequent analyses\n"
   echo -e "${RED}----ALL RECQUIRED PARAMETERS MUST BE PROVIDED IN THE CONFIG FILE (see ./config/config) ----${NC}\n"
   echo -e "use cases:

    ${bolditalic}${BLU}**** options **** : ${NC}${normal} \n
    -o 1 : perform whole inference where: \n
            step 1 : repeatmodeler + genome annotation + quality checks ; \n
            step 2 : GeneSpace  Gene Synteny + Whole Genome Synteny; \n 
            step 3 : Ds + related analyses  ;
    -o 2 : perform step 1 & 2 where :\n
            step 1 : repeatmodeler + genome annotation + quality checks ; \n
            step 2 : GeneSpace/Synteny 
    
    -o 3 :  perform step 2 & 3 where :\n
            step 2 : GeneSpace  Gene Synteny + Whole Genome Synteny; \n 
            step 3 : Ds + related analyses ;

    -o 4:  perform step 3 only where :\n
            step 3 : Ds + related analyses ;

    -o 5: performe step2 only where : 
            step 2: GeneSapce/Gene Synteny + Whole Genome Synteny
            
    -o 6 : perorm step1 only where :
           step  1 : repeatmodeler + genome annotation + quality checks ; \n
    
    All the details of the data (Genome/RNASeq/TEdata MUST be provided in the config file) \n
    see details in .infos/input_data_info.md
     \n
        ${bolditalic} - Analyses that will be performed: ${normal} 
        genome annotation, single copy orthologs inference, genespace, Dsplot, Rideogram, dotplot, circos, changepoint analysis
    
        "
}

############################################################
# Process the input options.                               #
############################################################

source ./config/config 

while [ $# -gt 0 ] ; do
  case $1 in
    -h  | --help ) Help ; exit 2 ;;
    -o  | --option ) option="$2" ; echo -e "options activated : ***${option}*** \n" >&2;;
   esac
   shift
done 

if [ -z $option ] ; then
    echo "Error no option provided ! I don't know what to do"
    Help
    exit 2
fi

###########################################################
#a few minor chek here:
###########################################################
#verify that at least a genome is provided - otherwise no work can be done at all:
if  [ -z "${genome1}" ]  && [ -z "${genome2}" ] ; 
then 
    echo "Error! provide the genome of at least one species"
    Help
    exit 2
else 
    if [ -n "${genome1}" ] ; 
    then
        echo "genome1 is $genome1" 
    fi
    if [ -n "${genome2}" ] ;
    then
        echo "genome 2 is $genome2 "
    fi
fi


# if no TE database then exit
# if no lineage for busco then exit and ask for it 
if [ -z "$TEdatabase" ] && [[ $annotate = YES ]] ;
then
    echo "Error NO Database for TE is provided "
    echo "I will not be able to assess mask TE prior to the genome annotation step "
    echo "this is very bad" 
    exit
    Help
fi


# if no lineage for busco then exit and ask for it 
if [ -z "${busco_lineage}" ]  && [[ $annotate = YES ]] ;
then
    echo "Error NO lineage provided for busco analyses" 
    echo "I will not be able to assess the quality of the annotation runs, which is compulsory for these anaylses"
    exit
    Help
fi

########################################
## Global variables
########################################

echo "------------------------------------------------------------"
echo "-----check all variables from the configuration file  ------"
echo "*** ancestral genome is ${ancestral_genome} ***"
echo "*** ancestral gff is ${ancestral_gff} ***"
echo "*** fasta for genome1 is ${genome1} **** "
echo "*** fasta for genome2 is ${genome2} **** "
echo "*** haplotype1 is ${haplotype1} ***"
echo "*** haplotype2 is ${haplotype2}"
echo "*** RNAseq? $rnaseq ***"
echo "*** RNAseq list is: ${RNAseqlist} ****" 
echo "*** Alternatively, BAM list might be : ${bamlist1} & ${bamlist2}****"
echo "*** gtf ? $gtf ***" 
echo "*** TE database is : $TEdatabase ***"
echo "*** lineage for busco analyses will be:  $busco_lineage ***"
echo "*** annotate? $annotate" 
echo "*** is species a fungus? $fungus ****"
echo -e "------------------------------------------------------------\n"



############################################################
# Generate architecture:
############################################################

echo "Generate Architecture"
mkdir -p haplo1/03_genome 

if [[ -z "${haplotype1}" ]] ; then
    haplotype1="haplo1"
fi

# ----- check compression of fasta  ------ ##

eval "$(conda shell.bash hook)"
conda activate braker_env

#check compression
if file --mime-type "$genome1" | grep -q gzip$; then
   echo "$genome1 is gzipped"
   gunzip "$genome1"
   genome1=${genome1%.gz}
   cd haplo1/03_genome 
   if [ -z ${haplotype1+x} ]; then
       echo "no haplotype name  provided - skip renaming scaffold ids"
       cp $genome1 .
   else
       echo "renaming scaffold ids "
       ../../00_scripts/00_rename_fasta.py $genome1 $current_name1 $haplotype1
   fi

   cd ../../
else
   echo "$genome1 is not gzipped"
   genome1=$genome1
   cd haplo1/03_genome 
   if [ -z ${haplotype1+x} ]; then
       echo "no haplotype name  provided - skip renaming scaffold ids"
       cp $genome1 .
   else
       echo "renaming scaffold ids "
       ../../00_scripts/00_rename_fasta.py $genome1 $current_name1 $haplotype1
   fi


   cd ../../
fi

#--- handling genome 2 -----
if [[ -n "$genome2" ]] ; then
    
    if [[ -n "${haplotype2}" ]] && [[ -n "${genome2}" ]]; then 
        mkdir -p haplo2/03_genome ; 
    fi
    if [[ -z "${haplotype2}" ]] && [[ -n "${genome2}" ]]; then 
        haplotype2="haplo2"
        mkdir -p haplo2/03_genome ; 
    fi
    
    #check genome compression:
    # ----- check compression of fasta  ------ ##
    #check compression
    if file --mime-type "$genome2" | grep -q gzip$; then
    echo "$genome2 is gzipped"
    gunzip "$genome2"
    genome2=${genome%.gz}
    cd haplo2/03_genome 
    if [ -z ${haplotype2+x} ]; then
        echo "no haplotype name  provided - skip renaming scaffold ids"
        cp $genome2 .
    else
        echo "renaming scaffold ids "
        ../../00_scripts/00_rename_fasta.py $genome2 $current_name2 $haplotype2
    fi
    
    cd ../../
    else
    echo "$genome2 is not gzipped"
    cd haplo2/03_genome 
    if [ -z ${haplotype2+x} ]; then
        echo "no haplotype name  provided - skip renaming scaffold ids"
        cp $genome2 .
    else
        echo "renaming scaffold ids "
        ../../00_scripts/00_rename_fasta.py $genome2 $current_name2 $haplotype2
    fi
    
    genome2=$genome2
    cd ../../
    fi
fi

conda deactivate

if [[ -n "${ancestral_genome}" ]] ; then 
    mkdir ancestral_sp  2>/dev/null
    #test if ancestral gtf is also provided in the config file:
    if [ -z "${ancestral_gff}" ] ; then 
        echo "error ! you provided an ancestral genome but no corresponding annotation (gff file)"
        echo "please correct that"
    fi

fi

#download uniprot for later checks
#./00_scripts/get_uniprot.sh #trivial as long as diamond successffuly installed 

if [[ $interpro = "YES" ]]
then
    echo "attempting to download and install inter-pro - this will take some time" 
    ./00_scripts/get_interpro.sh
fi

#artchitecture should be OK to proceed

############################################################
# Process the input options.                               #
############################################################

mkdir LOGS 2>/dev/null

echo -e "\ntesting use cases and checking inputs\n"



################################################################################
# ------------------ Section for option 1 (whole workflow) --------------------#
################################################################################

if [ $option == 1 ]; then 

    echo "----------------------------------------------------------------"
    echo "      the whole pipeline will be launched         " 
    echo "               checking config files settings                   "
    echo "----------------------------------------------------------------"


#if both species are provided without RNAseq:
elif [ -n "${genome1}" ] && [ -n "${genome2}" ]  && [[ $rnaseq = "NO" ]]  && [ -z "$ancestral_genome" ]  ; then
    echo "we will perform TE detection - genome annotation - Ds computation and plots"
    echo "genome are $genome1 and $genome2 "
    echo $(pwd)
    cd haplo1
    #partie en erreur à débuguer:
    if [ ! -z "$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
    
        echo "running TE detection and gene prediction"
        ../00_scripts/launch_step05_to_08.sh -g 03_genome/$haplotype1.fa  -s $haplotype1 -r NO -m YES -f YES  2>&1 |tee ../LOGS/log_step05_08_hap1
        #verify that alll run correctly 
    cd ../
    else
        echo "error no fasta file in 03_genome"
            echo "please copy your genome here"	
        Help
        exit 1
    fi
    
    cd haplo2
    if [ ! -z "$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
        echo "running TE detection and gene prediction"
        ../00_scripts/launch_step05_to_08.sh  -g 03_genome/$haplotype2.fa  -s $haplotype2 -r NO -m YES -f YES  2>&1 |tee ../LOGS/log_step05_08_hap2
    
        #verify that alll run correctly 
    cd ../
    else
        echo "error no fasta file in 03_genome"
            echo "please copy your genome here"	
        Help
        exit 1
    fi
    
    #if all is OK then run GeneSpace - paml etc :
    #modifiy the script RunGeneSpace etc to handle case with/without ancestral species
    ./00_scripts/11_run_GeneSpace_paml_ideogram.sh -s1 $haplotype1 -s2 $haplotype2 -c $scaffold 2>&1 |tee LOGS/log_GeneSpace_and_Co


#step1: 
elif [ -n "${genome1}" ] && [ -n "${genome2}" ] && [[ $rnaseq = "YES" ]]  && [ -n "${RNAseqlist}" ]  && [ -z "$ancestral_genome" ] && [ -z "$bamlist1" ] && [ -z "$bamlist2" ] ; then

    echo "we will perform TE detection - genome annotation with RNAseq - Ds computation and plots"
    echo "genomes are ${genome1} and ${genome2}"
    cd haplo1/
    if [ ! -z "$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
        ../00_scripts/launch_rnaseq.sh ${haplotype1}   2>&1 |tee ../LOGS/log_rna_haplo1
        #check that this script was sucessfull else kill:
        if [ $? -eq 0 ]; then
                echo rnaseq mapping succesffull
            echo "running TE detection and gene prediction"
            ../00_scripts/launch_step05_to_08.sh  -g 03_genome/$haplotype1.fa  -s $haplotype1 -r YES -m YES -f YES  2>&1 |tee ../LOGS/log_step05_08_hap1
    
    cd ../
        else
                echo ERROR - verfiy braker outputs!   
                exit 1
        fi
    else
        echo "error no fasta file in 03_genome"
            echo "please copy your genome here"	
        Help
        exit 1
    fi
    
    cd haplo2/
    if [ ! -z "$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
        ../00_scripts/launch_rnaseq.sh ${haplotype2}   2>&1 |tee ../LOGS/log_rna_haplo2
    
        #check that this script was sucessfull else kill:
        if [ $? -eq 0 ]; then
                echo rnaseq mapping succesffull
            echo "running TE detection and gene prediction"
            ../00_scripts/launch_step05_to_08.sh  -g 03_genome/$haplotype2.fa  -s $haplotype2 -r YES -m YES -f YES  2>&1 |tee ../LOGS/log_step05_08_hap2
    
    cd ../
        else
                echo ERROR - verfiy braker outputs!   
                exit 1
        fi
    else
        echo "error no fasta file in 03_genome"
            echo "please copy your genome here"	
        Help
        exit 1
    fi
    
    #then run GeneSpace etc :
    ./00_scripts/11_run_GeneSpace_paml_ideogram.sh -s1 $haplotype1 -s2 $haplotype2 -c $scaffold   2>&1 |tee LOGS/log_GeneSpace_and_Co


#step1: 
elif [ -n "${genome1}" ] && [ -n "${genome2}" ] && [[ $rnaseq = "YES" ]]  && [ -n "${bamlist1}" ] && [ -n "${bamlist2}" ] && [ -z "$ancestral_genome" ] ; then

    echo "we will perform TE detection - genome annotation with RNAseq - Ds computation and plots - bam files already provided"
    echo "genomes are ${genome1} and ${genome2}"
    cd haplo1/
    if [ ! -z "$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
        echo "running TE detection and gene prediction"
        ../00_scripts/launch_step05_to_08.sh  -g 03_genome/$haplotype1.fa  -s $haplotype1 -r YES -m YES -f YES  -b YES -b $bamlist1 2>&1 |tee ../LOGS/log_step05_08_hap1
    
        cd ../
    
    else
        echo "error no fasta file in 03_genome"
            echo "please copy your genome here"	
        Help
        exit 1
    fi
    
    cd haplo2/
    if [ ! -z "$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
        echo "running TE detection and gene prediction"
        ../00_scripts/launch_step05_to_08.sh  -g 03_genome/$haplotype2.fa  -s $haplotype2 -r YES -m YES -f YES  -b $bamlist2 2>&1 |tee ../LOGS/log_step05_08_hap2
    
    cd ../
    else
        echo "error no fasta file in 03_genome"
        echo "please copy your genome here"	
        Help
        exit 1
    fi
    
    #then run GeneSpace etc :
    ./00_scripts/11_run_GeneSpace_paml_ideogram.sh -s1 $haplotype1 -s2 $haplotype2 -c $scaffold   2>&1 |tee LOGS/log_GeneSpace_and_Co


#step1: 
elif [ -n "${genome1}" ] && [ -n "${genome2}" ] && [[ $rnaseq = "NO" ]] && [ -n "$ancestral_genome" ]  ; then
    echo "we will perform all analyses with annotations performed without rnaseq "
    echo "genomes are ${genome1} and ${genome2}"

    cd haplo1
    if [ ! -z "$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
        echo "running TE detection and gene prediction"
        ../00_scripts/launch_step05_to_08.sh  -g 03_genome/$haplotype1.fa  -s $haplotype1 -r NO -m YES -f YES  2>&1 |tee ../LOGS/log_step05_08_hap1
    
        #verify that alll run correctly 
    cd ../
    else
        echo "error no fasta file in 03_genome"
            echo "please copy your genome here"	
        help
        exit 1
    fi
    
    cd haplo2
    if [ ! -z "$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
        echo "running TE detection and gene prediction"
        ../00_scripts/launch_step05_to_08.sh  -g 03_genome/$haplotype2.fa  -s $haplotype2 -r NO -m YES -f YES  2>&1 |tee ../LOGS/log_step05_08_hap2
    
        #verify that alll run correctly 
    cd ../
    else
        echo "error no fasta file in 03_genome"
            echo "please copy your genome here"	
        help
        exit 1
    fi
    
    #then run GeneSpace etc :
    #modifiy the script RunGeneSpace etc to handle case with/without ancestral species
    ./00_scripts/11_run_GeneSpace_paml_ideogram.sh -s1 $haplotype1 -s2 $haplotype2 -a $ancestral_genome -g $ancestral_gff -c $scaffold   2>&1 |tee LOGS/log_GeneSpace_and_Co
    
    
#step1: 
elif [ -n "${genome1}" ] && [ -n "${genome2}" ] && [[ $rnaseq = "YES" ]] && [ -n "${bamlist1}" ] && [ -n "${bamlist2}" ] && [ -n "$ancestral_genome" ]  ; then
    echo "we will perform all analyses with annotations performed with rnaseq - list of bam provided "
    echo "genomes are ${genome1} and ${genome2}"

    cd haplo1
    if [ ! -z "$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
        echo "running TE detection and gene prediction"
        ../00_scripts/launch_step05_to_08.sh  -g 03_genome/$haplotype1.fa  -s $haplotype1 -r YES -m YES -f YES -b $bamlist1  2>&1 |tee ../LOGS/log_step05_08_hap1
    
        #verify that alll run correctly 
    cd ../
    else
        echo "error no fasta file in 03_genome"
            echo "please copy your genome here"	
        help
        exit 1
    fi
    
    cd haplo2
    if [ ! -z "$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
        echo "running TE detection and gene prediction"
        ../00_scripts/launch_step05_to_08.sh  -g 03_genome/$haplotype2.fa  -s $haplotype2 -r NO -m YES -f YES -b $bamlist2 2>&1 |tee ../LOGS/log_step05_08_hap2
    
        #verify that alll run correctly 
    cd ../
    else
        echo "error no fasta file in 03_genome"
            echo "please copy your genome here"	
        help
        exit 1
    fi
    
    #then run GeneSpace etc :
    #modifiy the script RunGeneSpace etc to handle case with/without ancestral species
    ./00_scripts/11_run_GeneSpace_paml_ideogram.sh -s1 $haplotype1 -s2 $haplotype2 -a $ancestral_genome -g $ancestral_gff -c $scaffold   2>&1 |tee LOGS/log_GeneSpace_and_Co
    
    

#step1: 
elif [ -n "${genome1}" ] && [ -n "${genome2}" ]  && [[ $rnaseq = "YES" ]] && [ -n "${RNAseqlist}" ]  && [ -n "$ancestral_genome" ] && [ -z "$bamlist1" ] && [ -z "$bamlist2" ]   ; then
    echo "we will perform all analyses including annotation with rnaseq"
    echo "genomes are ${genome1} and ${genome2}"
    cd haplo1/
    if [ ! -z "$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
        ../00_scripts/launch_rnaseq.sh ${haplotype1}   2>&1 |tee ../LOGS/log_rna_haplo1
    
        #check that this script was sucessfull else kill:
        if [ $? -eq 0 ]; then
                echo rnaseq mapping succesffull
                echo "running TE detection and gene prediction"
            ../00_scripts/launch_step05_to_08.sh  -g 03_genome/$haplotype1.fa  -s $haplotype1 -r YES -m YES -f YES  2>&1 |tee ../LOGS/log_step05_08_hap1
    
            cd ../
    
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
    if [ ! -z "$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
        ../00_scripts/launch_rnaseq.sh ${haplotype2}   2>&1 |tee ../LOGS/log_rna_haplo1
    
        #check that this script was sucessfull else kill:
        if [ $? -eq 0 ]; then
                echo rnaseq mapping succesffull
                echo "running TE detection and gene prediction"
            ../00_scripts/launch_step05_to_08.sh  -g 03_genome/$haplotype2.fa  -s $haplotype2 -r YES -m YES -f YES  2>&1 |tee ../LOGS/log_step05_08_hap2
    
            cd ../
    
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
    ./00_scripts/11_run_GeneSpace_paml_ideogram.sh -s1 $haplotype1 -s2 $haplotype2 -a $ancestral_genome -g $ancestral_gff -c $scaffold   2>&1 |tee LOGS/log_GeneSpace_and_Co


fi


################################################################################
# ------section for option 2 : TE + gene prediction + GeneSpace/Synteny--------#
################################################################################

if [ $option == 2 ] ; 
then 

    echo "----------------------------------------------------------------"
    echo "  TE + gene prediction + GeneSpace/Synteny will be launched     " 
    echo "               checking config files settings                   "
    echo "----------------------------------------------------------------"

    #insert here the test cases below: 

    #the test case are the same as option 1 for the gene prediction part
    
    #then we only to the geneSpace - so we must insert an option in script 11
    # to avoid running script 12 from within it. 
    # we still want minimap output for synteny so minimap should really be 
    # launched earlier


fi


################################################################################
# ------section for option 3 : GeneSpace/Synteny + Ds analyses ----------------#
################################################################################

################################################################################
# ----------------- section for option 4 : Ds analyses  -----------------------#
################################################################################

################################################################################
# -------------- section for option 5 : GeneSpace/Synteny analyses ------------#
################################################################################

if [ $option == 3 ] || [ $option == 4 ] || [ $option == 5 ]; then 

    #test input data:
    if [ -n "${gtf1}" ] && [ -n "${gtf2}" ] && [ -n "${genome1}" ] && [ -n "${genome2}" ] ; then
        echo "gtf and genomes file were provided" 
        echo "we will run geneSpace, compute Ds and other kind of analyses"
        mkdir -p haplo1/08_best_run haplo1/03_genome 
        mkdir -p haplo2/08_best_run haplo2/03_genome
        cp $gtf1 haplo1/08_best_run/${haplotype1}.gtf
        cp $gtf2 haplo2/08_best_run/${haplotype2}.gtf
 
        gffread -g $genome1 -w haplo1/08_best_run/"$haplotype1".spliced_cds.fa  $gtf1 
        gffread -g $genome2 -w haplo2/08_best_run/"$haplotype2".spliced_cds.fa  $gtf2 
        
        transeq -sequence haplo1/08_best_run/$haplotype1.spliced_cds.fa -outseq haplo1/08_best_run/"$haplotype1"_prot.fa
        transeq -sequence haplo2/08_best_run/$haplotype2.spliced_cds.fa -outseq haplo2/08_best_run/"$haplotype2"_prot.fa

        ##if [ -z "$ancestral_genome" ] ; then
        #
        #     #leave the variable empty
        #     
        #else [ -n "$ancestral_genome" ] ; then
        #    #do stuff related to ancestral genome
        #    #gffread transeq....
        #    #
        #fi
    fi
fi

#        echo ancestral_genome is "$ancestral_genome"
        
if [ $option == 3 ] ; then 

    echo "----------------------------------------------------------------"
    echo "       GeneSpace/Synteny + Ds analyses will be launched         " 
    echo "               checking config files settings                   "
    echo "----------------------------------------------------------------"

    opt = "snyteny_and_Ds"
    ./00_scripts/11_run_GeneSpace_paml_ideogram.sh -s1 $haplotype1 -s2 $haplotype2  -a $ancestral_genome -g $ancestral_gff -c $scaffold  -o $opt 2>&1 |tee LOGS/log_GeneSpace_and_Co

elif [ $option == 4 ] ; then 
            
    echo "----------------------------------------------------------------"
    echo "         only Ds + associated analyses will be launched         " 
    echo "               checking config files settings                   "
    echo "----------------------------------------------------------------"

    opt = "Ds_only"
    #note: in script 11 it would be important to check if a geneSpace run exist already
    ./00_scripts/11_run_GeneSpace_paml_ideogram.sh -s1 $haplotype1 -s2 $haplotype2  -a $ancestral_genome -g $ancestral_gff -c $scaffold  -o $opt 2>&1 |tee LOGS/log_GeneSpace_and_Co

elif [ $option == 5 ] ; then 
        
    echo "----------------------------------------------------------------"
    echo "          only GeneSpace/Synteny  analyses will be launched     " 
    echo "               checking config files settings                   "
    echo "----------------------------------------------------------------"

    opt = "snyteny_ony"
    ./00_scripts/11_run_GeneSpace_paml_ideogram.sh -s1 $haplotype1 -s2 $haplotype2  -a $ancestral_genome -g $ancestral_gff -c $scaffold  -o $opt 2>&1 |tee LOGS/log_GeneSpace_and_Co
    #check success here
            
fi


################################################################################
# -------------- section for option 6: TE + gene Prediction only --------------#
################################################################################

if [ $option = 6 ]; then 

    echo "----------------------------------------------------------------"
    echo "          Only TE + gene prediction  will be performed          " 
    echo "              checking config files settings                    "
    echo "----------------------------------------------------------------"


    # option = 6 - only genome 1 - No RNAseq : 
    if [ -n "${genome1}" ] && [ -z "${genome2}" ]  && [[ $rnaseq = "NO" ]]   && [ -z "$ancestral_genome" ] ; then
        cd haplo1
        if [ ! -z "$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
            echo "only the genome of one species was provided" 
            echo "genome is ${genome1} "
            echo "running TE detection and gene prediction"
            ../00_scripts/launch_step05_to_08.sh -g 03_genome/$haplotype1.fa  -s $haplotype1 -r YES -m YES -f YES  2>&1 |tee ../LOGS/log_step05_08_hap1
        
        cd ../
        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"	
            Help
            exit 1
        fi
    
        if [ -n "${genome1}" ] && [ -z "${genome2}" ]  && [[ $rnaseq = "NO" ]]   && [ -n "$ancestral_genome" ] ; then
        cd haplo1
        if [ ! -z "$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
            echo "only the genome of one species was provided" 
            echo "genome is ${genome1} "
            echo "running TE detection and gene prediction"
            ../00_scripts/launch_step05_to_08.sh -g 03_genome/$haplotype1.fa  -s $haplotype1 -r YES -m YES -f YES  2>&1 |tee ../LOGS/log_step05_08_hap1
        
        cd ../
        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"	
            Help
            exit 1
        fi
        
    # option = 6 - only genome 1 - No RNAseq : 
    elif [ -z "${genome1}" ] && [ -n "${genome2}" ] && [[ $rnaseq = "NO" ]] && [ -z "$ancestral_genome" ]   ; then
        cd haplo2
        if [ ! -z "$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
            echo "only the genome of one species was provided" 
            echo "genome is ${genome2} "
            echo "running TE detection and gene prediction"
            ../00_scripts/launch_step05_to_08.sh -g 03_genome/$haplotype2.fa  -s $haplotype2 -r YES -m YES -f YES  2>&1 |tee ../LOGS/log_step05_08_hap2
        
        cd ../
        else
            echo "error no fasta file in 03_genome"
                echo "please copy your genome here"	
            Help
            exit 1
        fi
    
    elif [ -z "${genome1}" ] && [ -n "${genome2}" ] && [[ $rnaseq = "NO" ]] && [ -n "$ancestral_genome" ]   ; then
        cd haplo2
        if [ ! -z "$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
            echo "only the genome of one species was provided" 
            echo "genome is ${genome2} "
            echo "running TE detection and gene prediction"
            ../00_scripts/launch_step05_to_08.sh -g 03_genome/$haplotype2.fa  -s $haplotype2 -r YES -m YES -f YES  2>&1 |tee ../LOGS/log_step05_08_hap2
        
        cd ../
        else
            echo "error no fasta file in 03_genome"
                echo "please copy your genome here"	
            Help
            exit 1
        fi

    # option = 6 - only genome 2 & RNAseq : 
    elif [ -n "${genome1}" ] && [ -z "${genome2}" ] && [[ $rnaseq = "YES" ]] && [ -n "${RNAseqlist}" ]  && [ -z "$ancestral_genome" ] && [ -z "$bamlist1" ] && [ -z "$bamlist2" ] ; then
        cd haplo1/
        if [ ! -z "$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
            echo "only the genome of one species was provided with RNAseq data" 
            echo "genome is ${genome1} "
            ../00_scripts/launch_rnaseq.sh ${haplotype1}    2>&1 |tee ../LOGS/log_rna_haplo1
    
            #check that this script was sucessfull else kill:
            if [ $? -eq 0 ]; then
                echo rnaseq mapping succesffull
                echo "running TE detection and gene prediction"
                ../00_scripts/launch_step05_to_08.sh -g 03_genome/$haplotype1.fa  -s $haplotype1 -r YES -m YES -f YES  2>&1 |tee ../LOGS/log_step05_08_hap1
    
                cd ../
    
            else
                echo ERROR - verfiy braker outputs!   
                exit 1
            fi
        else
                echo "error no fasta file in 03_genome"
                echo "please copy your genome here"	
            Help
            exit 1
        fi
    
    
    # option = 6 - only genome 1 & RNAseq from bamfiles: 
    elif [ -n "${genome1}" ] && [ -z "${genome2}" ]  && [[ $rnaseq = "YES" ]]  && [ -z "$ancestral_genome" ] && [ -n "$bamlist1" ] && [ -z "$bamlist2" ] ; then
        cd haplo1/
        if [ ! -z "$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
            echo "only the genome of one species was provided along with a list of bam" 
            echo "genome is ${genome2} "
            echo "running TE detection and gene prediction"
            ../00_scripts/launch_step05_to_08.sh -g 03_genome/$haplotype1.fa  -s $haplotype2 -r YES -m YES -f YES -b $bamlist1  2>&1 |tee ../LOGS/log_step05_08_hap2
        
            cd ../
        
        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"	
            Help
            exit 1
        fi
    
    
    # option = 6 - only genome 2 & RNAseq :
    elif [ -z "${genome1}" ] && [ -n "${genome2}" ]  && [[ $rnaseq = "YES" ]]  && [ -z "$ancestral_genome" ] && [ -z "$bamlist1" ] && [ -z "$bamlist2" ]  ; then
        cd haplo2/
        if [ ! -z "$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
            echo "only the genome of one species was provided with RNAseq data" 
            echo "genome is ${genome2} "
            ../00_scripts/launch_rnaseq.sh ${haplotype2}   2>&1 |tee ../LOGS/log_rna_haplo2
        
        
                #check that this script was sucessfull else kill:
            if [ $? -eq 0 ]; then
                    echo rnaseq mapping succesffull
                echo "running TE detection and gene prediction"
            ../00_scripts/launch_step05_to_08.sh -g 03_genome/$haplotype2.fa  -s $haplotype2 -r YES -m YES -f YES  2>&1 |tee ../LOGS/log_step05_08_hap2
        
                cd ../
        
            else
                    echo ERROR - verfiy braker outputs!   
                    exit 1
            fi
        else
            echo "error no fasta file in 03_genome"
                echo "please copy your genome here"	
            Help
            exit 1
        fi
    
    # option = 6 - only genome 2 & RNAseq from bam file:
    elif [ -z "${genome1}" ] && [ -n "${genome2}" ]  && [[ $rnaseq = "YES" ]] && [ -z "$ancestral_genome" ]  && [ -z "$bamlist1" ] && [ -n "${bamlist2}" ]  ; then
        cd haplo2/
        if [ ! -z "$(ls -A 03_genome/ |grep -v Readme )"  ] ; then
            echo "only the genome of one species was provided along with a list of bam" 
            echo "genome is ${genome2} "
            echo "running TE detection and gene prediction"
            ../00_scripts/launch_step05_to_08.sh -g 03_genome/$haplotype2.fa  -s $haplotype2 -r YES -m YES -f YES  -b $bamlist2 2>&1 |tee ../LOGS/log_step05_08_hap2
        
            cd ../
        
        else
            echo "error no fasta file in 03_genome"
            echo "please copy your genome here"	
            Help
            exit 1
        fi
    
    fi

fi
