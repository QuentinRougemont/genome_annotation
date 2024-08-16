#!/usr/bin/env python3

#purpose: easy code to rename braker output \
#        by inserting the species name in the last column

"""Module providing a function printing python version."""
import sys


INPUT_FILE = sys.argv[1]
SPECIES_NAME = sys.argv[2]
W1 = open(f'{SPECIES_NAME}.IDchecked.gtf', 'w')

F1 = open(INPUT_FILE, 'r')
for line in F1:
    line = line.strip()     #removes the \n
    line = line.split('\t') #splits the line at the '\t' into a list
    contig = line[0].split('_')[1]
    if SPECIES_NAME != line[0].split('_')[0]:
        print('problem of concordance between species name \
              and the first column in gtf')
    #changes the name of the contig here:
    new_contig = SPECIES_NAME+'_'+contig

    #spliting the last columns according to ";"
    categories = line[-1].split(';')

    new_categories = []

    for i in range(0, len(categories)):
        category = line[-1].split(';')[i]
        if category == '': #doesn't read the empty element when the line finishes by ";"
            continue
        category_id = category.split('"')[1] #keeps only the id
        category_name = category.split('"')[0] #keeps only the name of the category \
                #("transcript_id" or "gene_id")
        new_category_id = SPECIES_NAME+'_'+contig+'_'+category_id.split('_')[-1]
        if len(categories) == 1: #if only one element
            new_categories.append(category_name+' "'+new_category_id+'"')
        else: #if 2 elements, put the ";" after each element
            new_categories.append(category_name+' "'+new_category_id+'"'+';')

    #prints the line into the file
    print('\t'.join([new_contig]+line[1:8]+[(' ').join(new_categories)]), file=W1)
