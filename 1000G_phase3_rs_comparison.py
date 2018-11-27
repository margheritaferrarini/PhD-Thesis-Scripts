#!/usr/bin/env python3 

##THIS SCRIPTS COMPARES rs OF 1000G PHASE3 VARIANTS IN GRCh37 e GRCh38

f=open("/lustre/projects/exome_variants/april/1000G_analysis/1000G_phase3_comparison.tsv", 'r')
only37 = open("/lustre/projects/exome_variants/april/1000G_analysis/variants_only_GRCh37.tsv", 'w')
only38 = open("/lustre/projects/exome_variants/april/1000G_analysis/variants_only_GRCh38.tsv", 'w')
corrected_var = open("/lustre/projects/exome_variants/april/1000G_analysis/corrected_variants.tsv", 'w')
uncorrected_var = open("/lustre/projects/exome_variants/april/1000G_analysis/uncorrected_variants.tsv", 'w')
alternative_allele = open("/lustre/projects/exome_variants/april/1000G_analysis/alternative_allele_to_correct.tsv", 'w')

previousLine = f.readline()                                   ##first line, equal to previous
previousLine = previousLine.rstrip('\n').split('\t')          ##split first line


while True:
    previousCHROM = previousLine[0]                           ##save previous chromosome
    previousPOS = previousLine[1]                             ##save previous position  
    previousRS = previousLine[2]                              ##save previous rs                                                                                                                 
    previousREF = previousLine[3]                             ##save previous reference   
    previousGRCh = previousLine[6]                            ##save column with 7 or 8 where 7=GRCh37 and 8=GRCh38

    currentLine = f.readline()                                ##read following line
    if not currentLine:    ##at the end of file, get out of loop
        break
    currentLine = currentLine.rstrip('\n').split('\t')        ##split current line
    currentCHROM = currentLine[0]                             ##save current chromosome
    currentPOS = currentLine[1]                               ##save current position 
    currentRS = currentLine[2]                                ##save current rs                                                                                                                     
    currentREF = currentLine[3]                               ##save current reference   
    currentGRCh = currentLine[6]                              ##save column with 7 or 8 where 7=GRCh37 and 8=GRCh38


    if (currentRS == previousRS) and (currentRS != '.'):      ##compare rs
        
        ##rs equal in different reference 
        if currentGRCh != previousGRCh:                       ##compare reference (NB: unexpectedly the same rs can be found in different chromosomes)
            if currentREF != previousREF:                         ##reference alleles are different, corrected variant
                corrected_var.write('{}\n{}\n'.format('\t'.join(previousLine), '\t'.join(currentLine)))
            if currentREF == previousREF:                         ##reference alleles are equal, not corrected variant 
                uncorrected_var.write('{}\n{}\n'.format('\t'.join(previousLine), '\t'.join(currentLine))) 
            previousLine = f.readline()                           ##read following line
            if not previousLine:
                break
            previousLine = previousLine.rstrip('\n').split('\t')

        ##rs equal in the same reference due to errors or to aleternative alleles
        else:
            ##alternative alleles of the same variant
            if (currentCHROM == previousCHROM) and (currentPOS == previousPOS):    ##alternative alleles in the same reference with the same rs, print in a new file
                alternative_allele.write('{}\n'.format('\t'.join(previousLine)))   ##print first line
                previousLine = currentLine   ##current line is now the previous line

            ##errors with different chromosome or position, but the same rs
            else:
                if currentGRCh=="7":
                    only37.write('{}\n{}\n'.format('\t'.join(previousLine), '\t'.join(currentLine)))                                                                                 
                if currentGRCh=="8":
                    only38.write('{}\n{}\n'.format('\t'.join(previousLine), '\t'.join(currentLine))) 
                previousLine = f.readline()
                if not previousLine:
                    break
                previousLine = previousLine.rstrip('\n').split('\t')
#                currentLine = currentLine.rstrip('\n').split('\t')
#                currentRS = currentLine[2]
#                currentREF = currentLine[3]
#                if (currentRS == previousRS) and (currentRS != '.'):
#                    if currentREF != previousREF:
#                        corrected_var.write('{}\n{}\n'.format('\t'.join(previousLine), '\t'.join(currentLine)))
#                    if currentREF == previousREF:
#                        uncorrected_var.write('{}\n{}\n'.format('\t'.join(previousLine), '\t'.join(currentLine)))
#                    previousLine = f.readline()
#                    if not previousLine:
#                        break
#                    previousLine = previousLine.rstrip('\n').split('\t')
    else:
        if previousGRCh=="7":
            only37.write('{}\n'.format('\t'.join(previousLine)))     ##rs only in GRCh37                                                                                                 
        if previousGRCh=="8":
            only38.write('{}\n'.format('\t'.join(previousLine)))     ##rs only in GRCh38
        previousLine = currentLine

            
f.close()
only37.close()
only38.close()
corrected_var.close()
uncorrected_var.close()
alternative_allele.close()



