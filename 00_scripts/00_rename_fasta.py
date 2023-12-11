#!/usr/bin/env python3

##RENAME THE FASTA ID IN THE GENOME##
import sys
from Bio import SeqIO
import re
input_file=sys.argv[1]
species_name=sys.argv[2]

print(species_name)
species_name=re.split('[-._]',species_name)
print(species_name)
species_name=''.join(species_name)
print(species_name)

output_file=f'{species_name}.fa'

w1=open(f'{output_file}','w')
#reads the fasta file
for contig in SeqIO.parse(f'{input_file}', "fasta"):
	#for each contig, contig_name = the first element encoded after the ">" (if there is len etc, not included)
	contig_name=contig.id
	print(f'old contig name: {contig_name}')
	#contig_name=re.split('[-._]',contig_name)
	contig_name=re.split('[-._]',contig_name)[-1]
	contig_name=''.join(contig_name)
	# in case there is informations on length etc, removes it to keep just the contig name
	new_contig_name=species_name+'_'+contig_name
	print(f'new contig name: {new_contig_name}')

	print('>'+new_contig_name, file=w1)
	print(contig.seq, file=w1)
w1.close()


