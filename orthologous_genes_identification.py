#!/usr/bin/env python3
# orthologous_genes_identification.py    ##script name

import argparse, os, Bio
from Bio.Blast.Applications import NcbiblastnCommandline

## INPUT FILES
parser=argparse.ArgumentParser(description="This script identifies orthologous genes in different 
                                           species starting from the FASTA sequence of a human gene")
parser.add_argument("-q", "--query_sequence", help="GRCh38 gene sequence in FASTA format", required=True)
parser.add_argument("-c", "--coordinates", help="GRCh38 gene coordinates in BED format", required=True)
parser.add_argument("-g", "--genomes", help="list of genome files, each in FASTA format", nargs='+', required=True) 
args=parser.parse_args()
gene_positions=args.coordinates


## For each genome in command line
for genome in args.genomes:

## 1. Save organism name by reading first line of genome file 
    genome_file=open(genome, 'r')                 ## open genome file
    first_line=genome_file.readline()             ## read first line
    first_line=first_line.rstrip('\n')            
    first_line=first_line.split(' ')              ## split first line 
    genus=first_line[1]                           ## save genus
    species=first_line[2]                         ## save species
    organism=("{}_{}").format(genus, species)     ## name
    genome_file.close()
    print (organism)

## 2. Blast on genome file
    blastn_cline=NcbiblastnCommandline(query=args.query_sequence, db=genome, outfmt=" '7 std sseq' ", out=organism + "_blast_results.tsv", soft_masking=True) 
    ## '7 std sseq' to obtain a tabular output file with comments line, with standard information and with sequence of aligned part of subject sequence
    ## max_target_seqs=10 save only first 10 best results
    ## soft_masking true as suggested by literature
    os.system(str(blastn_cline))                  ## save output 

## 3. Compare line to save in FASTA format only the result with highest bit-score or lowest e-value
    blast_results=open(organism + "_blast_results.tsv", 'r')      ## open output file
    fasta_sequence=open(organism + "_ortholog_sequence.fa", 'w')  ## open a new file
    for line in blast_results:
        line=line.rstrip('\n')
        if not line.startswith('#'):                              ## if not header
            first_line=line.split('\t')                          
            first_e_value=float(first_line[10])                   ## save lowest e-value
            first_bit_score=float(first_line[11])                 ## save highest bit-score
            subject_id=organism + "," + first_line[1]             ## save ID
            best_hit_sequence=first_line[12]                      ## save sequence
            best_hit_sequence_without_gaps=best_hit_sequence.replace("-", "")                       ## remove gap
            fasta_sequence.write(">{}\n{}\n".format(subject_id, best_hit_sequence_without_gaps))    ## print sequence
            for line in blast_results:                            ## read following lines
                if not line.startswith('#'):                      ## last line starts with #
                    line=line.split('\t')
                    e_value=float(line[10])                       ## save e-value
                    bit_score=float(line[11])                     ## save bit-score
                    if bit_score > first_bit_score:               ## if bit-score if higher than the first one
                        subject_id=organism + "," + line[1]                ## save ID
                        best_hit_sequence=line[12]                         ## save sequence
                        best_hit_sequence_without_gaps=best_hit_sequence.replace("-", "")                       ## remove gap
                        fasta_sequence.write(">{}\n{}\n".format(subject_id, best_hit_sequence_without_gaps))    ## print sequence      
                    if bit_score == first_bit_score:              ## if bit-score equal to the first one
                        if e_value == first_e_value:              ## if e-value equal to the first one
                            subject_id=organism + "," + line[1]            ## save ID
                            best_hit_sequence=line[12]                     ## save sequence
                            best_hit_sequence_without_gaps=best_hit_sequence.replace("-", "")                   ## remove gap
                            fasta_sequence.write(">{}\n{}\n".format(subject_id, best_hit_sequence_without_gaps))## print sequence
    blast_results.close()
    fasta_sequence.close()

## 4. Blast best hit sequence to human genome
    hsa_blastn_cline=NcbiblastnCommandline(query=organism + "_ortholog_sequence.fa", db="/lustre/database/genomes/homo_sapiens/GRCh38/INDEXED/BLAST/Homo_sapiens.GRCh38.91.dna.primary_assembly.fa", outfmt=" '7 std sseq' ", out=organism + "_hsa_blast_results.tsv", soft_masking=True)
    os.system(str(hsa_blastn_cline))                                       ## save output

## 5. Compare coordinates of gene and s. start e s. end of first hit 
    ## Read output
    hsa_blast_results=open(organism + "_hsa_blast_results.tsv", 'r')       ## open output file                      
    for line in hsa_blast_results:
        if not line.startswith('#'):                                       ## pass line starting with #
            line=line.rstrip('\n')
            line=line.split('\t')
            subject_chromosome=line[1]                                     ## save chromosome
            subject_start=int(line[8])                                     ## save starting position
            subject_end=int(line[9])                                       ## salva ending position
            break                                                          ## get out of cicle
    hsa_blast_results.close()

    ## Save gene coordinates
    coordinates=open(gene_positions, 'r')                                  ## open file with coordinates
    for line in coordinates:
        line=line.rstrip('\n')
        line=line.split('\t')
        gene_chromosome=line[0][-1]                                        ## save chromosome
        gene_start=int(line[1])                                            ## save starting position
        gene_end=int(line[2])                                              ## save ending position
    coordinates.close()

    ## Compare chromosome and coordinates
    if subject_chromosome==gene_chromosome:                                ## equal chromosome
        if subject_start>subject_end:                                      ## if gene on reverse strand
            if (subject_start<=gene_end) and (subject_end>=gene_start):    ## if coordinates are included 
                print ("the blast result is the orthologous gene")
            else:
                print ("the blast result is NOT the orthologous gene")
        if subject_start<subject_end:                                      ## if gene on forward strand                                                                       
            if (subject_start>=gene_start) and (subject_end<=gene_end):    ## if coordinates are inclued                                                                      
                print ("the blast result is the orthologous gene")
            else:
                print ("the blast result is NOT the orthologous gene")
    else:
        print ("the blast result is NOT the orthologous gene")
