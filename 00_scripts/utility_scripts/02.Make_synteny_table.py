#!/usr/bin/env python3

#Make a big table with all the DS to make a plot in R
import sys

## 1st arg = haplo1
## 2nd arg = haplo2
## 3rd arg = path to orthofinder results
## 4rd arg = path to bed
## 5th arg = ancestral species ==TRUE ?
## 6th arg = ancestral species name
hap1=sys.argv[1]
hap2=sys.argv[2]
path=sys.argv[3]
path_bed=sys.argv[4]
is_Sp_anc=sys.argv[5]

### 1)  GET THE GENE INFORMATION FROM THE BED FILES	
###---------------------------------------------------------------------------------------
##get the info from bed for ancestral species
if is_Sp_anc=='TRUE':
	sp_anc=sys.argv[6]
	d_gene={}
	bed=open(f'{path_bed}{sp_anc}.bed','r')
	for line in bed:
		line=line.split()
		gene=line[3]
		chrom=line[0]
		start=line[1]
		end=line[2]
		d_gene[gene]=[chrom,start,end]
	bed.close()

##get the info from bed for hap1
d_gene1={}
bed1=open(f'{path_bed}{hap1}.bed','r')
for line in bed1:
	line=line.split()
	gene=line[3]
	chrom=line[0]
	start=line[1]
	end=line[2]
	d_gene1[gene]=[chrom,start,end]
bed1.close()

##get the info from bed for hap2
d_gene2={}
bed2=open(f'{path_bed}{hap2}.bed','r')
for line in bed2:
	line=line.split()
	gene=line[3]
	chrom=line[0]
	start=line[1]
	end=line[2]
	d_gene2[gene]=[chrom,start,end] 
bed2.close()

#print(d_gene)


### 2)  GET THE SINGLE COPY ORTHOLOGS LIST
###---------------------------------------------------------------------------------------
## PUT THE SINGLE ORTHOGROUPS INTO A SET
set_single_ortho=set()
f0=open(f'{path}/Orthogroups/Orthogroups_SingleCopyOrthologues.txt','r')
for line in f0:
	line=line.strip()
	set_single_ortho.add(line)
f0.close()

### 3)  MAKE 1-1 TABLES FOR SYNTENY
###---------------------------------------------------------------------------------------
### A- hap 1 vs hap 2
###-------------------
d={}
f1=open(f'{path}/Orthologues/Orthologues_{hap1}/{hap1}__v__{hap2}.tsv','r')
f1.readline()
for line in f1:
	line=line.split('\t')
	OG=line[0]
	col_hap1=line[1].strip()
	col_hap2=line[2].strip()
	if len(col_hap1.split(', ')) !=1 or len(col_hap2.split(', ')) !=1:
		continue
	contig_hap1=col_hap1.split('_')[1]
	gene_hap1=col_hap1
	contig_hap2=col_hap2.split('_')[1]
	gene_hap2=col_hap2

	if OG in set_single_ortho:
		try:
			d[OG].append([gene_hap1,contig_hap1,d_gene1[gene_hap1][1],d_gene1[gene_hap1][2],gene_hap2,contig_hap2,d_gene2[gene_hap2][1],d_gene2[gene_hap2][2]])
		except KeyError:
			try:
				d[OG]=[]
				d[OG].append([gene_hap1,contig_hap1,d_gene1[gene_hap1][1],d_gene1[gene_hap1][2],gene_hap2,contig_hap2,d_gene2[gene_hap2][1],d_gene2[gene_hap2][2]])
			except KeyError:
				print('no key', gene_hap2)
	else:
		continue

f1.close()


w1=open(f'synteny_{hap1}_{hap2}.txt','w')
header=['ortho_gp','gene1','chrom1','start1','end1','gene2','chrom2','start2','end2']
print('\t'.join(header),file=w1)
for OG in d.keys():
	for i in range(0,len(d[OG])):
		new_line=[OG]+d[OG][i]
		if len(new_line)!=9:
			print('problem',new_line)
			continue
		print('\t'.join(new_line), file=w1)

### B- ancestral species vs hap 1 and hap 2
###-----------------------------------------
if is_Sp_anc=='TRUE':
	sp_anc=sys.argv[6]
	d_all={}
	for hap in [hap1, hap2]:
		d={}
		f1=open(f'{path}/Orthologues/Orthologues_{sp_anc}/{sp_anc}__v__{hap}.tsv','r')
		f1.readline()
		for line in f1:
			line=line.split('\t')
			OG=line[0]
			col_sp_anc=line[1].strip()
			col_hap=line[2].strip()
			if len(col_sp_anc.split(', ')) !=1 or len(col_hap.split(', ')) !=1:
				continue
			contig_sp_anc=col_sp_anc.split('_')[1]
			gene_sp_anc=col_sp_anc
			contig_hap=col_hap.split('_')[1]
			gene_hap=col_hap
			if OG in set_single_ortho:
				try:
					if hap==hap1:
						d[OG].append([gene_sp_anc,contig_sp_anc,d_gene[gene_sp_anc][1],d_gene[gene_sp_anc][2],gene_hap,contig_hap,d_gene1[gene_hap][1],d_gene1[gene_hap][2]])
					elif hap==hap2:
						d[OG].append([gene_sp_anc,contig_sp_anc,d_gene[gene_sp_anc][1],d_gene[gene_sp_anc][2],gene_hap,contig_hap,d_gene2[gene_hap][1],d_gene2[gene_hap][2]])

				except KeyError:
					try:
						if hap==hap1:
							d[OG]=[]
							d[OG].append([gene_sp_anc,contig_sp_anc,d_gene[gene_sp_anc][1],d_gene[gene_sp_anc][2],gene_hap,contig_hap,d_gene1[gene_hap][1],d_gene1[gene_hap][2]])
						elif hap==hap2:
							d[OG]=[]
							d[OG].append([gene_sp_anc,contig_sp_anc,d_gene[gene_sp_anc][1],d_gene[gene_sp_anc][2],gene_hap,contig_hap,d_gene2[gene_hap][1],d_gene2[gene_hap][2]])
					except KeyError:
						print('no key', gene_hap)
				try:
					if hap==hap1:
						d_all[OG].extend([hap1,gene_hap,contig_hap,d_gene1[gene_hap][1],d_gene1[gene_hap][2]])
					elif hap==hap2:
						d_all[OG].extend([hap2,gene_hap,contig_hap,d_gene2[gene_hap][1],d_gene2[gene_hap][2]])
				except KeyError:
						d_all[OG]=[sp_anc,gene_sp_anc,contig_sp_anc,d_gene[gene_sp_anc][1],d_gene[gene_sp_anc][2]]
						if hap==hap1:
							d_all[OG].extend([hap1,gene_hap,contig_hap,d_gene1[gene_hap][1],d_gene1[gene_hap][2]])
						elif hap==hap2:
							d_all[OG].extend([hap2,gene_hap,contig_hap,d_gene2[gene_hap][1],d_gene2[gene_hap][2]])
			else:
				continue

		f1.close()

		w1=open(f'synteny_{sp_anc}_{hap}.txt','w')
		header=['ortho_gp','gene1','chrom1','start1','end1','gene2','chrom2','start2','end2']
		print('\t'.join(header),file=w1)
		for OG in d.keys():
			for i in range(0,len(d[OG])):
				new_line=[OG]+d[OG][i]
				if len(new_line)!=9:
					print('problem',new_line)
					continue
				print('\t'.join(new_line), file=w1)
		w1.close()

	### C-  Make synteny tables with the 3 haplotypes	
	###---------------------------------------------------------------------------------------
	w1=open(f'synteny_{sp_anc}_{hap1}_{hap2}.txt','w')
	header=['OG','hap0','gene0','chrom0','start0','end0','hap1','gene1','contig1','start1','end1','hap2','contig2','start2','end2']
	print('\t'.join(header), file=w1)
	for OG in d_all.keys():
		new_line=[OG]+d_all[OG]
		print('\t'.join(new_line), file=w1)
	w1.close()
