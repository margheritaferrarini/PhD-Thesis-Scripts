#!/usr/bin/python3

from pyliftover import LiftOver
lo = LiftOver('hg19', 'Hg38')

file_to_convert=open("/lustre/projects/exome_variants/april/1000G_analysis/comparing_rs/total_GRCh37_only_GRCh37.tsv", 'r')
file_converted=open("/lustre/projects/exome_variants/april/1000G_analysis/comparing_rs/total_GRCh37_only_GRCh37_convertedCoordinates.tsv", 'w')
file_errors=open("/lustre/projects/exome_variants/april/1000G_analysis/comparing_rs/total_GRCh37_only_GRCh37_notConvertible.tsv", 'w')
for line in file_to_convert:
    line=line.rstrip('\n').split('\t')
    chrom='chr{}'.format(line[0])
    pos=int(line[1])
    conversion_list=lo.convert_coordinate(chrom, pos)
    if not conversion_list:  ##The list may be empty (locus is deleted in the new assembly)
	file_errors.write('{}\n'.format('\t'.join(line))) 
    else:
        new_pos=str(conversion_list).split(', ')[1]
        file_converted.write('{}\t{}\t{}\n'.format(line[0], new_pos, '\t'.join(line[2:])))
file_to_convert.close()
file_converted.close()
file_errors.close()
