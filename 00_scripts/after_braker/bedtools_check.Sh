bedtools intersect -a <(awk '$3=="transcript" {print $1,$4,$5,$7,$10}' OFS="\t" *gtf | sort -k1,1 -k2,2n) -b <(awk '$3=="transcript" {print $1,$4,$5,$7,$10}' OFS="\t" *gtf | sort -k1,1 -k2,2n) -wa -wb | awk '$5!=$NF && $4==$(NF-1) && ($2==$(NF-3) || $3==$(NF-2))' |wc -l

