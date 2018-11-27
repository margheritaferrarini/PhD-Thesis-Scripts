#!/usr/bin/env python3 

to_correct=open("/lustre/projects/exome_variants/april/1000G_analysis/comparing_rs/alternative_allele_to_correct_NEW_GRCh37.tsv", 'r')

corrected=open("/lustre/projects/exome_variants/april/1000G_analysis/comparing_rs/corrected_variants_NEW.tsv", 'r')
uncorrected=open("/lustre/projects/exome_variants/april/1000G_analysis/comparing_rs/uncorrected_variants_NEW.tsv", 'r')
only_GRCH37=open("/lustre/projects/exome_variants/april/1000G_analysis/comparing_rs/variants_only_GRCh37_NEW.tsv", 'r')

alternative_corrected=open("/lustre/projects/exome_variants/april/1000G_analysis/comparing_rs/alternative_allele_corrected_variants_NEW.tsv", 'w')
alternative_uncorrected=open("/lustre/projects/exome_variants/april/1000G_analysis/comparing_rs/alternative_allele_uncorrected_variants_NEW.tsv", 'w')
alternative_only_GRCH37=open("/lustre/projects/exome_variants/april/1000G_analysis/comparing_rs/alternative_allele_only_GRCh37_NEW.tsv", 'w')

corrected_list=[]
uncorrected_list=[]
only_GRCH37_list=[]

for line in corrected:
    line=line.rstrip('\n').split('\t')
    if line[6]=='7':
        corrected_list.append('{}|{}'.format(line[0], line[1]))

for line in uncorrected:
    line=line.rstrip('\n').split('\t')
    if line[6]=='7':
        uncorrected_list.append('{}|{}'.format(line[0], line[1]))

for line in only_GRCH37:
    line=line.rstrip('\n').split('\t')
    if line[6]=='7':
        only_GRCH37_list.append('{}|{}'.format(line[0], line[1]))

for line in to_correct:
    line=line.rstrip('\n').split('\t')
    var='{}|{}'.format(line[0], line[1])
    if var in corrected_list and var not in uncorrected_list and var not in only_GRCH37_list:
        alternative_corrected.write('{}\n'.format('\t'.join(line)))
    if var in uncorrected_list and var not in corrected_list and var not in only_GRCH37_list:
        alternative_uncorrected.write('{}\n'.format('\t'.join(line)))
    if var in only_GRCH37_list and var not in corrected_list and var not in uncorrected_list:
        alternative_only_GRCH37.write('{}\n'.format('\t'.join(line)))
    if var in corrected_list and var in uncorrected_list:
        print ('\t'.join(line))
    if var in corrected_list and var in only_GRCH37_list:
        print ('\t'.join(line))
    if var in uncorrected_list and var in only_GRCH37_list:
        print ('\t'.join(line))

to_correct.close()
corrected.close()
uncorrected.close()
only_GRCH37.close()
alternative_corrected.close()
alternative_uncorrected.close()
alternative_only_GRCH37.close()
