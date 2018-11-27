#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy

corrected=open("Supplemental_VCF_S1.vcf", 'r')

var_list=[]

for line in corrected:
    if not line.startswith('#'):
        line=line.rstrip('\n')
        line=line.split('\t')
        var=('{}\t{}\t{}\t{}'.format(line[0],line[1],line[3],line[4]))
        var_list.append(var)

var_set=set(var_list)

corrected.close()

all_var=open("ALL.wgs.integrated_phase1_v3.20101123.snps_indels_sv.sites.vcf", 'r')

AF_list=[]

for line in all_var:
    if not line.startswith('#'):
        line=line.rstrip('\n')
        line=line.split('\t')
        var=('{}\t{}\t{}\t{}'.format(line[0],line[1],line[3],line[4]))
        if var in var_set:
            info=line[7].split(';')
            for item in info:
                if item.startswith('AF='):
                    item=item.split('=')
                    AF=float(item[1])
                    AF_list.append(AF)
print(len(AF_list))
all_var.close()

plt.hist(AF_list, 100, log=True, color='steelblue', edgecolor='none')
plt.xlabel('Frequency')
plt.xticks(numpy.arange(0,1.1,step=0.1))
plt.ylabel('Updated Bases')
plt.title('Frequency of Bases Updated in GRCh38')
plt.minorticks_off()
plt.savefig('corrected_bases_frequences_log')
plt.savefig('corrected_bases_frequences_log.pdf', format='pdf')

