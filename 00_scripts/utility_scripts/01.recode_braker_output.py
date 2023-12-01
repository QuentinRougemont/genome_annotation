#!/usr/bin/env python3

import sys

input_file=sys.argv[1]
species_name=sys.argv[2]
w1=open(f'{species_name}.IDchecked.gtf','w')

f1=open(input_file,'r')
for line in f1:
	line=line.strip() 	#removes the \n
	line=line.split('\t') #splits the line at the '\t' into a list
	contig=line[0].split('_')[1]
	if species_name != line[0].split('_')[0]:
		print('problem of concordance between species name and the first column in gtf')
	#changes the name of the contig here:
	new_contig=species_name+'_'+contig

	#spliting the last columns according to ";"
	categories=line[-1].split(';')

	new_categories=[]

	for i in range(0,len(categories)):
		category=line[-1].split(';')[i]
		if category=='': #doesn't read the empty element when the line finishes by ";"
			continue
		category_id=category.split('"')[1] #keeps only the id
		category_name=category.split('"')[0] #keeps only the name of the category ("transcript_id" or "gene_id")
		new_category_id=species_name+'_'+contig+'_'+category_id.split('_')[-1]
		if len(categories)==1: #if only one element
			new_categories.append(category_name+' "'+new_category_id+'"')
		else: #if 2 elements, put the ";" after each element
			new_categories.append(category_name+' "'+new_category_id+'"'+';')

	#prints the line into the file
	print('\t'.join([new_contig]+line[1:8]+[(' ').join(new_categories)]), file=w1)



