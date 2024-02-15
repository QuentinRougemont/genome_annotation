#!/usr/bin/env python3

##RENAME THE FASTA ID IN THE GENOME##
import sys
from Bio import SeqIO
import re
input_file=sys.argv[1]
current_basename=sys.argv[2]
new_basename=sys.argv[3]

output_file=f'{new_basename}.fa'

w1=open(f'{output_file}','w')
#reads the fasta file
for contig in SeqIO.parse(f'{input_file}', "fasta"):
    contig_name=contig.id
    print(f'old contig name: {contig_name}')
    contig_number=re.split(current_basename, contig_name)[1]
    new_contig_name=new_basename+contig_number
    print(f'new contig name: {new_contig_name}')

    print('>'+new_contig_name, file=w1)
    print(contig.seq, file=w1)
w1.close()
