#!/usr/bin/python3

##THIS SCRIPT COUNTS THE NUMBER OF MAIR IN GRCH37 AND GRCh38 FILES

CORRECTED=open("/lustre/projects/exome_variants/april/1000G_analysis/CORRECTED_IN_GRCH38.tsv", "r")
UNCORRECTED=open("/lustre/projects/exome_variants/april/1000G_analysis/UNCORRECTED_IN_GRCH38.tsv", "r")

CORRECTED_MAIR=[]
UNCORRECTED_MAIR=[]

for line in CORRECTED:
    line=line.rstrip('\n').split('\t')
    if line[0]=='Y':
        AF=((((line[5].split(';'))[2]).split('='))[1]).split(',')                                             
    else:
        AF=((((line[5].split(';'))[1]).split('='))[1]).split(',')                                               
    for i in AF:                                                                                               
        i=float(i)
        if i>0.5:
            CORRECTED_MAIR.append(i)   
print (len(CORRECTED_MAIR))

for line in UNCORRECTED:
    line=line.rstrip('\n').split('\t')
    if line[0]=='Y':
        AF=((((line[5].split(';'))[2]).split('='))[1]).split(',')
    else:
        AF=((((line[5].split(';'))[1]).split('='))[1]).split(',')
    for i in AF:
        i=float(i)
        if i>0.5:
            UNCORRECTED_MAIR.append(i)
print (len(UNCORRECTED_MAIR))

CORRECTED.close()
UNCORRECTED.close()
