#!/bin/bash

#Author: QR
#Date: 11-2022
#script to detect repeated sequences


#------------- EXTERNAL VARIABLE FROM CONFIG FILE -------------- #
source ../config/config

eval "$(conda shell.bash hook)"
conda activate repeatmodeler_env


echo TEdatabase is $TEdatabase
echo "NCBI species is $ncbi_species"
echo "genome is $genome"
#echo "database is $haplotype1"

if [ $# -ne 3  ]; then
    echo "USAGE: $0 reference_genome database name rm_unknwon(yes/no)"
    echo -e "Expecting the following parameters:\n
          1 - the reference genome\n
          2 - a basename for database building (e.g. 'myfavoritespecies')\n
          3 - a string YES/NO about wether unknwon repeat should be removed\n\n"
    exit 1
else
    genome=$1
    database=$2
    rm_unknown=$3
    echo -e "reference genome is $genome \n"
    echo -e "database is : ${database}\n"
    echo -e "rm_unknown option is set to $rm_unknown"
    echo " "
fi


base=$(basename $genome)
#--------------DECLARE THE USUAL GENERIQ STUFF: -----------------#
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="99_log_files"
#create log folder
mkdir $LOG_FOLDER 2>/dev/null

cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"

# ----- check compression of fasta  ------ ##
#check compression
if file --mime-type "$genome" | grep -q gzip$; then
   echo "$genome is gzipped"
   gunzip "$genome"
   genome=${genome%.gz}
else
   echo "$genome is not gzipped"
   genome=$genome
fi

#this should be deprecated:
#sed 's/ [0-9A-Za-z=-]*//g' $genome > ${genome%.fa}.simpl.fa
#genome=${genome%.fa}.simpl.fa

base=$(basename $genome )

mkdir 05_TE 2>/dev/null
cd 05_TE

#--------------STEP1 : RUN REPEATMODELER  -----------------------#

##build db:
#BuildDatabase -name $database -engine ncbi ../$genome 2>&1 | tee ../$LOG_FOLDER/buildDatabase.$base.$TIMESTAMP.log

#de novo TE annotations:
#RepeatModeler -threads 18 -engine ncbi -database $database 2>&1 | tee ../$LOG_FOLDER/repeatmodeler_$base.$TIMESTAMP.log


#--------------STEP2 : RUN REPEATMASKER -------------------------#

# BASED ON DATABASE : 
FOLDER1=FOLDER1_"${base}"_mask.$TIMESTAMP
mkdir $FOLDER1
lib1=$TEdatabase 
RepeatMasker -pa 18 -e ncbi -lib $lib1 -xsmall -dir "$FOLDER1" ../$genome 2>&1 | tee ../$LOG_FOLDER/F1_repeatmasker_$base.$TIMESTAMP.log


# Based on de-novo repeat + database:
FOLDER2=FOLDER2_"${base}"_mask.$TIMESTAMP

# test if we keep Unknwon repeat or not
## without Unknwon repeat ##
if [[ $rm_unknown = "YES" ]]
then
        echo "removing Unknown TE Repeats ..."
        awk '$0~/^>/{if(NR>1){print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' $database-families.fa |\
        sed -e '/Unknown/,+1d' |\
        cat $lib1 - > $base.no_unknown.repbase.fa
	libcat="$base".no_unknown.repbase.fa
else
	#with known repeat
        echo "keep all candidate TEs... "
        cat $database-families.fa $lib1 > $base.repbase.fa
	libcat="$base".repbase.fa
fi

#run repeatmasker:
RepeatMasker -pa 18 -e ncbi -lib $libcat -xsmall -dir "$FOLDER2" "$FOLDER1"/"$base".masked 2>&1 |\
	tee ../$LOG_FOLDER/F2_repeatmasker_$base.$TIMESTAMP.log


## ----- step 2.3: based on online data ----- ## 
#online database
FOLDER3=FOLDER3_"${base}"_mask.$TIMESTAMP
mkdir "$FOLDER3"

#run repeatmasker:
RepeatMasker -pa 18 -e ncbi -species  "${ncbi_species}" -xsmall -dir "$FOLDER3"   "$FOLDER2"/"$base".masked.masked 2>&1 | \
	tee ../$LOG_FOLDER/F3_repeatmasker_$base.$TIMESTAMP.log

cd ../03_genome

if [[ $rm_unknown = "YES" ]]
then
        ln -s ../05_TE/$FOLDER3/$base.masked.masked.masked genome.wholemask_no_unknown.fa
else
        ln -s ../05_TE/$FOLDER3/$base.masked.masked.masked genome.wholemask.fa
fi

