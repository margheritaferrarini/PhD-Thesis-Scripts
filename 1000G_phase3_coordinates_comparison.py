#!/usr/bin/env python3 

##THIS SCRIPTS COMPARES COORDINATES OF 1000G PHASE3 VARIANTS IN GRCh37 e GRCh38 

f=open("/lustre/projects/exome_variants/april/1000G_analysis/comparing_rs/1000G_phase3_NEWcomparison.tsv", 'r')
only37 = open("/lustre/projects/exome_variants/april/1000G_analysis/comparing_rs/variants_only_GRCh37_NEW.tsv", 'w')
only38 = open("/lustre/projects/exome_variants/april/1000G_analysis/comparing_rs/variants_only_GRCh38_NEW.tsv", 'w')
corrected_var = open("/lustre/projects/exome_variants/april/1000G_analysis/comparing_rs/corrected_variants_NEW.tsv", 'w')
uncorrected_var = open("/lustre/projects/exome_variants/april/1000G_analysis/comparing_rs/uncorrected_variants_NEW.tsv", 'w')
alternative_allele = open("/lustre/projects/exome_variants/april/1000G_analysis/comparing_rs/alternative_allele_to_correct_NEW.tsv", 'w')

previousLine = f.readline()                                   ##first line, equal to previous
previousLine = previousLine.rstrip('\n').split('\t')          ##split first line 


while True:
    previousCHROM = previousLine[0]                           ##save previous chromosome
    previousPOS = previousLine[1]                             ##save previous position  
    previousREF = previousLine[3]                             ##save previous reference   
    previousGRCh = previousLine[6]                            ##save column with 7 or 8 where 7=GRCh37 and 8=GRCh38

    currentLine = f.readline()                                ##read following line
    if not currentLine:    ##at the end of file, get out of loop
        break
    currentLine = currentLine.rstrip('\n').split('\t')        ##split current line
    currentCHROM = currentLine[0]                             ##save current chromosome
    currentPOS = currentLine[1]                               ##save current position 
    currentREF = currentLine[3]                               ##save current reference  
    currentGRCh = currentLine[6]                              ##save column with 7 or 8 where 7=GRCh37 and 8=GRCh38

    ##if positions are equal
    if (currentCHROM == previousCHROM) and (currentPOS == previousPOS):  

        ##alternative alleles if the same reference --> save in alternative_allele
        if currentGRCh == previousGRCh: 
            alternative_allele.write('{}\n'.format('\t'.join(previousLine)))   ##print first line 
            previousLine = currentLine   ##current line is now the previous line

        ##variants in the two references --> compare reference alleles 
        else:
            if currentREF != previousREF:                         ##reference alleles are different, corrected variant                                                                    
                corrected_var.write('{}\n{}\n'.format('\t'.join(previousLine), '\t'.join(currentLine)))
            if currentREF == previousREF:                         ##reference alleles are equal, not corrected variant                                                                 
                uncorrected_var.write('{}\n{}\n'.format('\t'.join(previousLine), '\t'.join(currentLine)))
            previousLine = f.readline()                           ##read following line                                                                                                     
            if not previousLine:
                break
            previousLine = previousLine.rstrip('\n').split('\t')

    ##if posotions are different --> save in 37 or 38
    else:
        if previousGRCh=="7":
            only37.write('{}\n'.format('\t'.join(previousLine)))     ##present only in GRCh37                                                                                                   
        if previousGRCh=="8":
            only38.write('{}\n'.format('\t'.join(previousLine)))     ##present only in GRCh38                                                                                                   
        previousLine = currentLine

f.close()
only37.close()
only38.close()
corrected_var.close()
uncorrected_var.close()
alternative_allele.close()



