
cd 01_GTF
mkdir START_REMOVED
mkdir TMP

for i in $(ls *gtf ) ; do 

        n=$(basename $i ) ;
        awk '$3=="transcript" {print $1"_"$4"_"$5"\t"$1"_"$4"\t"$1"_"$5"\t"$0}'  $i > TMP/tmp.$n 
        awk 'NR == FNR {count[$2]++; next} count[$2]>0 {print $2"\t"$7"\t"$8"\t"$13"\t"$8-$7}' TMP/tmp.$n TMP/tmp.$n |awk '$5>max[$1]{max[$1]=$5; row[$1]=$4} END{for (i in row) print row[i]}' > TMP/longest.to.keep.tmp.$n ;
        grep -Ff TMP/longest.to.keep.tmp.$n $i > START_REMOVED/${n%.longest_transcript.gtf}.nonduplicated.gtf 
done 

mkdir ALL_REMOVED 2>/dev/null
mkdir TMP2
for i in $(ls START_REMOVED/*gtf ) ; do 
        #remove gene with exact same chr and end!
        n=$(basename $i ) ;
        awk '$3=="transcript" {print $1"_"$4"_"$5"\t"$1"_"$4"\t"$1"_"$5"\t"$0}'  $i > TMP2/tmp.$n 
        awk 'NR == FNR {count[$3]++; next} count[$3]>0 {print $3"\t"$7"\t"$8"\t"$13"\t"$8-$7}' TMP2/tmp.$n TMP2/tmp.$n |awk '$5>max[$1]{max[$1]=$5; row[$1]=$4} END{for (i in row) print row[i]}' > TMP2/longest.to.keep.tmp.$n ;
        grep -Ff TMP2/longest.to.keep.tmp.$n $i > ALL_REMOVED/${n%.gtf}_2.gtf 
done

cd ALL_REMOVED

ls *gtf |sed 's/.nonduplicated_2.gtf//g' > list.of.gtf

for i in $(cat list.of.gtf ) ; 
do 
        gffread -g ../../02_Genome/$i.fa -x $i.spliced_cds.fa $i.nonduplicated_2.gtf  
        transeq -sequence $i.spliced_cds.fa -outseq "$i"_prot.fa ; 
done
