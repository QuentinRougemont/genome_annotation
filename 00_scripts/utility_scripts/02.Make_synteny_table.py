#!/usr/bin/env python3

#Make a big table with all the DS to make a plot in R
"""Module providing a function printing python version."""
import sys

## 1st arg = haplo1
## 2nd arg = haplo2
## 3rd arg = PATH to orthofinder results
## 4rd arg = PATH to bed
## 5th arg = ancestral species  == TRUE ?
## 6th arg = ancestral species name
HAP1 = sys.argv[1]
HAP2 = sys.argv[2]
PATH = sys.argv[3]
PATH_BED = sys.argv[4]
IS_SP_ANC = sys.argv[5]

### 1)  GET THE GENE INFORMATION FROM THE BED FILES
###---------------------------------------------------------------------------------------
##get the info from bed for ancestral species
if IS_SP_ANC == 'TRUE':
    SP_ANC = sys.argv[6]
    D_GENE = {}
    bed = open(f'{PATH_BED}{SP_ANC}.bed', 'r')
    for line in bed:
        line = line.split()
        gene = line[3]
        chrom = line[0]
        start = line[1]
        end = line[2]
        D_GENE[gene] = [chrom, start, end]
    bed.close()

##get the info from bed for HAP1
D_GENE1 = {}
bed1 = open(f'{PATH_BED}{HAP1}.bed', 'r')
for line in bed1:
    line = line.split()
    gene = line[3]
    chrom = line[0]
    start = line[1]
    end = line[2]
    D_GENE1[gene] = [chrom, start, end]
bed1.close()

##get the info from bed for HAP2
D_GENE2 = {}
bed2 = open(f'{PATH_BED}{HAP2}.bed', 'r')
for line in bed2:
    line = line.split()
    gene = line[3]
    chrom = line[0]
    start = line[1]
    end = line[2]
    D_GENE2[gene] = [chrom, start, end]
bed2.close()

#print(D_GENE)


### 2)  GET THE SINGLE COPY ORTHOLOGS LIST
###---------------------------------------------------------------------------------------
## PUT THE SINGLE ORTHOGROUPS INTO A SET
SET_SINGLET_ORTHO = set()
F0 = open(f'{PATH}/Orthogroups/Orthogroups_SingleCopyOrthologues.txt', 'r')
for line in F0:
    line = line.strip()
    SET_SINGLET_ORTHO.add(line)
F0.close()

### 3)  MAKE 1-1 TABLES FOR SYNTENY
###---------------------------------------------------------------------------------------
### A- hap 1 vs hap 2
###-------------------
D = {}
F1 = open(f'{PATH}/Orthologues/Orthologues_{HAP1}/{HAP1}__v__{HAP2}.tsv', 'r')
F1.readline()
for line in F1:
    line = line.split('\t')
    OG = line[0]
    col_HAP1 = line[1].strip()
    col_HAP2 = line[2].strip()
    if len(col_HAP1.split(', ')) != 1 or len(col_HAP2.split(', ')) != 1:
        continue
    contig_HAP1 = col_HAP1.split('_')[1]
    gene_HAP1 = col_HAP1
    contig_HAP2 = col_HAP2.split('_')[1]
    gene_HAP2 = col_HAP2

    if OG in SET_SINGLET_ORTHO:
        try:
            D[OG].append([gene_HAP1, contig_HAP1, D_GENE1[gene_HAP1][1], D_GENE1[gene_HAP1][2],
                          gene_HAP2, contig_HAP2, D_GENE2[gene_HAP2][1], D_GENE2[gene_HAP2][2]])
        except KeyError:
            try:
                D[OG] = []
                D[OG].append([gene_HAP1, contig_HAP1, D_GENE1[gene_HAP1][1], D_GENE1[gene_HAP1][2],
                              gene_HAP2, contig_HAP2, D_GENE2[gene_HAP2][1], D_GENE2[gene_HAP2][2]])
            except KeyError:
                print('no key', gene_HAP2)
    else:
        continue

F1.close()


W1 = open(f'synteny_{HAP1}_{HAP2}.txt', 'w')
HEADER = ['ortho_gp', 'gene1', 'chrom1', 'start1', 'end1', 'gene2', 'chrom2', 'start2', 'end2']
print('\t'.join(HEADER), file=W1)
for OG in D.keys():
    for i in range(0, len(D[OG])):
        new_line = [OG]+D[OG][i]
        if len(new_line) != 9:
            print('problem', new_line)
            continue
        print('\t'.join(new_line), file=W1)

### B- ancestral species vs hap 1 and hap 2
###-----------------------------------------
if IS_SP_ANC == 'TRUE':
    SP_ANC = sys .argv[6]
    D_ALL = {}
    for hap in [HAP1, HAP2]:
        D = {}
        f1 = open(f'{PATH}/Orthologues/Orthologues_{SP_ANC}/{SP_ANC}__v__{hap}.tsv', 'r')
        f1.readline()
        for line in f1:
            line = line.split('\t')
            OG = line[0]
            col_SP_ANC = line[1].strip()
            col_hap = line[2].strip()
            if len(col_SP_ANC.split(', ')) != 1 or len(col_hap.split(', ')) != 1:
                continue
            contig_SP_ANC = col_SP_ANC.split('_')[1]
            gene_SP_ANC = col_SP_ANC
            contig_hap = col_hap.split('_')[1]
            gene_hap = col_hap
            if OG in SET_SINGLET_ORTHO:
                try:
                    if hap == HAP1:
                        D[OG].append([gene_SP_ANC, contig_SP_ANC, D_GENE[gene_SP_ANC][1],
                                      D_GENE[gene_SP_ANC][2], gene_hap,
                                      contig_hap, D_GENE1[gene_hap][1], D_GENE1[gene_hap][2]])
                    elif hap == HAP2:
                        D[OG].append([gene_SP_ANC, contig_SP_ANC, D_GENE[gene_SP_ANC][1],
                                      D_GENE[gene_SP_ANC][2], gene_hap,
                                      contig_hap, D_GENE2[gene_hap][1], D_GENE2[gene_hap][2]])

                except KeyError:
                    try:
                        if hap == HAP1:
                            D[OG] = []
                            D[OG].append([gene_SP_ANC, contig_SP_ANC, D_GENE[gene_SP_ANC][1],
                                          D_GENE[gene_SP_ANC][2], gene_hap,
                                          contig_hap, D_GENE1[gene_hap][1], D_GENE1[gene_hap][2]])
                        elif hap == HAP2:
                            D[OG] = []
                            D[OG].append([gene_SP_ANC, contig_SP_ANC, D_GENE[gene_SP_ANC][1],
                                          D_GENE[gene_SP_ANC][2], gene_hap,
                                          contig_hap, D_GENE2[gene_hap][1], D_GENE2[gene_hap][2]])
                    except KeyError:
                        print('no key', gene_hap)
                try:
                    if hap == HAP1:
                        D_ALL[OG].extend([HAP1, gene_hap,
                                          contig_hap, D_GENE1[gene_hap][1], D_GENE1[gene_hap][2]])
                    elif hap == HAP2:
                        D_ALL[OG].extend([HAP2, gene_hap,
                                          contig_hap, D_GENE2[gene_hap][1], D_GENE2[gene_hap][2]])
                except KeyError:
                        D_ALL[OG] = [SP_ANC, gene_SP_ANC,
                                     contig_SP_ANC, D_GENE[gene_SP_ANC][1], D_GENE[gene_SP_ANC][2]]
                        if hap == HAP1:
                            D_ALL[OG].extend([HAP1, gene_hap,
                                              contig_hap, D_GENE1[gene_hap][1],
                                              D_GENE1[gene_hap][2]])
                        elif hap == HAP2:
                            D_ALL[OG].extend([HAP2, gene_hap,
                                              contig_hap, D_GENE2[gene_hap][1],
                                              D_GENE2[gene_hap][2]])
            else:
                continue

        f1.close()

        W1 = open(f'synteny_{SP_ANC}_{hap}.txt', 'w')
        HEADER = ['ortho_gp', 'gene1', 'chrom1', 'start1', 'end1',
                  'gene2', 'chrom2', 'start2', 'end2']

        print('\t'.join(HEADER), file=W1)
        for OG in D.keys():
            for i in range(0, len(D[OG])):
                new_line = [OG]+D[OG][i]
                if len(new_line) != 9:
                    print('problem', new_line)
                    continue
                print('\t'.join(new_line), file=W1)
        W1.close()

    ### C-  Make synteny tables with the 3 haplotypes
    ###---------------------------------------------------------------------------------------
    W1 = open(f'synteny_{SP_ANC}_{HAP1}_{HAP2}.txt', 'w')
    HEADER = ['OG', 'hap0', 'gene0', 'chrom0', 'start0', 'end0', 'HAP1',
              'gene1', 'contig1', 'start1', 'end1', 'HAP2', 'contig2', 'start2', 'end2']

    print('\t'.join(HEADER), file=W1)
    for OG in D_ALL.keys():
        new_line = [OG]+D_ALL[OG]
        print('\t'.join(new_line), file=W1)
    W1.close()
