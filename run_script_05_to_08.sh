#!/bin/bash
#Date:13-02-23
#Author QR
# WARNING STILL BETA #! 

#deprecated -- this needs to be fully updated


##  ------------------------ general parameters --------------------------------  ##
while [ $# -gt 0 ] ; do
  case $1 in
    -g | --genome) genome="$2" ;echo "the genome file  is: $genome" >&2;;
    -s | --haplotype) haplotype="$2" ;echo "the haplotype name will be $haplotype" >&2;;
    -m | --mask )   Mask="$2" ; echo "unknown TE will be removed after repeatemasking" >&2;;
    -r | --rna )    RNAseq="$2" ; echo "Is RNAseq provided ? $RNAseq " >&2;; 
    -f | --fungus ) fungus="$2" ; echo "Do species belong to fungi? $fungus" >&2;;
    -h | --help) echo -e "Option required:
    -g/--genome \t the reference genome file 
    -s/--haplotype\t the haplotype name (used for database building and basename in several steps)
    Optional:
    -m/--mask\t a string stating wether unknown TE should be removed for genome annotation (YES/NO) -- default: YES
    -r/--ref \t a string stating wether RNAseq data should be used (YES/NO) -- default: NO
    " >&2;exit 1;;
    esac
    shift
done

if [ -z "$genome" ] || [ -z "$haplotype" ] ; then
	echo >&2 "Fatal error: Ref genome (-g), and haplotype name (-s) not defined\n
	see manual with -h or --help"
exit 2
fi

#optional parameters:
if [ -z "$Mask" ] ; then
    Mask=YES
fi

if [ -z "$RNAseq" ] ; then 
  RNAseq=NO
fi

# ------------------ run RepeatModeler and RepeatMasker ------------------  ##

../00_scripts/05_repeatmodeler.sh "$genome" "$haplotype" "$Mask" 2>&1 |tee log_rm
if [ $? -eq 0 ]; then
    echo repeatmodeler run successfull
else
    echo repeatmodeler failed. check the provided libraries and whether all software and dependancies are correctly installed   
    exit 
fi

# -------------------- run Braker  ---------------------------- #
../00_scripts/06_braker.sh 03_genome/genome.wholemask_no_unknown.fa $haplotype $RNAseq $fungus 2>&1 |tee braker_log  #NO for no rnaseq  

# -------------------- run Busco  ---------------------------- #
if [[ $RNAseq = "YES" ]]
then
    ../00_scripts/07_busco_after_braker.sh $busco_lineage YES
else
   ../00_scripts/07_busco_after_braker.sh $busco_lineage 
fi

# ---------------------reshape Braker output ----------------- #
if [[ $RNAseq = "YES" ]]
then
    ../00_scripts/08_braker_reshaping.sh -s  $haplotype -r YES
else
   ../00_scripts/08_braker_reshaping.sh -s $haplotype -r NO
fi
if [ $? -eq 0 ]; then
    echo braker output successfully processed
else
    echo ERROR - verfiy braker outputs!   
    exit 1
fi


