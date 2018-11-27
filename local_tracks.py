#!/usr/bin/env python3

import sys
import re
import math

def cigar(cigar):    # it returns the actual length of the alignment
	match_len=0
	curr_val=0
	for p in range(len(cigar)):
		c = cigar[p]
		if c == '*':
			return 0
		v=ord(c)-48
		if re.search(c, "MDN=X"):
			match_len += curr_val
			curr_val=0
		elif v>=0 and v<=9:
			curr_val = curr_val * 10 + v
		else:
			curr_val=0
	return match_len

############### default parameters ###############
infile=""
indexfile=""
max_insert_len = 8000
min_insert_len = 200
skip=0

################ get parameters ##################
if len(sys.argv) <3:
	print("Sintax: " + sys.argv[0] + " -f samfile -c chromosome -s start -e end -m max_insert_len -l min_insert_len")

	sys.exit()

for item in range(len(sys.argv)):
	if skip>0:
		skip=skip-1
		continue
	if sys.argv[item] == "-f":
		infile=sys.argv[item+1]
		skip=1
	if sys.argv[item] == "-c":
		chrom=sys.argv[item+1]
		skip=1
	if sys.argv[item] == "-s":
		target_start=int(sys.argv[item+1])
		skip=1
	if sys.argv[item] == "-e":
		target_end=int(sys.argv[item+1])
		skip=1
	if sys.argv[item] == "-m":
		max_insert_len=int(sys.argv[item+1])
		skip=1
	if sys.argv[item] == "-l":
		min_insert_len=int(sys.argv[item+1])
		skip=1
	if sys.argv[item] == "-i":
		indexfile=sys.argv[item+1]
		linenumber=int(sys.argv[item+2])
		skip=2

################# variable initialization ####################
physical_coverage_change = [0] * (1 + target_end - target_start)
sequence_coverage_change = [0] * (1 + target_end - target_start)
sum_insert_length_change = [0] * (1 + target_end - target_start)
single_mapper_change = [0] * (1 + target_end - target_start)
mates_in_change = [0] * (1 + target_end - target_start) # mates directed in (like in a pcr)
mates_out_change = [0] * (1 + target_end - target_start) # mates directed out
mates_same_dir_change = [0] * (1 + target_end - target_start)
reads_for_change = [0] * (1 + target_end - target_start)
reads_rev_change = [0] * (1 + target_end - target_start)
reads_with_HS_change = [0] * (1 + target_end - target_start)
mate_to_long = [0] * (1 + target_end - target_start)
flag_256 = [0] * (1 + target_end - target_start)

################# open input file or stdin ####################
if len(infile)>0:
    inf = open(infile)
else:
    inf = sys.stdin

################# read sam file ####################

if len(indexfile) >0:
	print("indexfile not yet implemented")
	exit(0)

for line in inf:
	line = line.rstrip()
	if line.startswith('@'):
		continue
	field = line.split("\t")
	if field[2] != chrom :      # wrong chromosome
		continue
	insert_len = int(field[8]) + cigar(field[5])
	start_pos = int(field[3])
	end_pos = start_pos + insert_len
	flag=int(field[1])

	############# calculate physical coverage and insert length #################
	if ((end_pos>=target_start) and (start_pos<=target_end) and ((flag&3)==3) and (insert_len>min_insert_len) and (insert_len<max_insert_len)):
		p1=start_pos
		p2=end_pos
		if p1 < target_start:
			p1 = target_start
		if p2 > target_end:
			p2 = target_end
		physical_coverage_change[p1 - target_start] += 1    # increment start position by one
		physical_coverage_change[p2 - target_start] -= 1    # decrement end position by one
		sum_insert_length_change[p1 - target_start] += insert_len # increment start position by ins len
		sum_insert_length_change[p2 - target_start] -= insert_len # decrement end position by ins len
		if((flag&48) == 16): # mates directed out
			mates_out_change[p1 - target_start] += 1
			mates_out_change[p2 - target_start] -= 1
		elif((flag&48) == 32): # mates directed in
			mates_in_change[p1 - target_start] += 1
			mates_in_change[p2 - target_start] -= 1
		else:
			mates_same_dir_change[p1 - target_start] += 1
			mates_same_dir_change[p2 - target_start] -= 1

	############# calculate sequence coverage etc on individual reads #################
	p1=start_pos
	p2=start_pos + cigar(field[5])
	if ((p2 >= target_start) and (p1 <= target_end) ):   # insert maps within target region
		if p1 < target_start:
			p1 = target_start
		if p2 > target_end:
			p2 = target_end
		if(insert_len > max_insert_len):  # weirdo1 mate_to_long
			mate_to_long[p1 - target_start] += 1
			mate_to_long[p2 - target_start] -= 1
		if(flag > 255):  # weirdo2   flag_256
			flag_256[p1 - target_start] += 1
			flag_256[p2 - target_start] -= 1
		if ((flag & 4) == 0):  # this read maps properly, update sequence coverage
			sequence_coverage_change[p1 - target_start] += 1    # increment start position by one
			sequence_coverage_change[p2 - target_start] -= 1    # decrement end position by one
		if((flag & 12) == 8):
			single_mapper_change[p1 - target_start] += 1    # increment start position by one
			single_mapper_change[p2 - target_start] -= 1    # decrement end position by one
		if((flag & 16) == 0):
			reads_for_change[p1 - target_start] += 1
			reads_for_change[p2 - target_start] -= 1
		else:
			reads_rev_change[p1 - target_start] += 1
			reads_rev_change[p2 - target_start] -= 1
		if (re.search("H", field[5]) or re.search("S", field[5])):
			reads_with_HS_change[p1 - target_start] += 1
			reads_with_HS_change[p2 - target_start] -= 1
inf.closed


################################### output wig files #####################################

# physical coverage
f1 = open('physical.wig', 'w')
f1.write("fixedStep chrom=" + chrom + " start=" + str(target_start) + " step=1 span=1\n")
current_coverage = 0
for position in range(target_end-target_start):
	current_coverage += physical_coverage_change[position]
	f1.write(str(current_coverage) + '\n')
f1.flush()
f1.closed

# mates_in coverage
f1 = open('mates_in.wig', 'w')
f1.write("fixedStep chrom=" + chrom + " start=" + str(target_start) + " step=1 span=1\n")
current_coverage = 0
for position in range(target_end-target_start):
	current_coverage += mates_in_change[position]
	f1.write(str(current_coverage) + '\n')
f1.flush()
f1.closed

# mates_out coverage
f1 = open('mates_out.wig', 'w')
f1.write("fixedStep chrom=" + chrom + " start=" + str(target_start) + " step=1 span=1\n")
current_coverage = 0
for position in range(target_end-target_start):
	current_coverage += mates_out_change[position]
	f1.write(str(current_coverage) + '\n')
f1.flush()
f1.closed

# mates_same_direction
f1 = open('mates_same_direction.wig', 'w')
f1.write("fixedStep chrom=" + chrom + " start=" + str(target_start) + " step=1 span=1\n")
current_coverage = 0
for position in range(target_end-target_start):
	current_coverage += mates_same_dir_change[position]
	f1.write(str(current_coverage) + '\n')
f1.flush()
f1.closed

# average insert length coverage
f1 = open('av_insert_length.wig', 'w')
f1.write("fixedStep chrom=" + chrom + " start=" + str(target_start) + " step=1 span=1\n")
current_coverage = 0
current_sum=0
for position in range(target_end-target_start):
	current_coverage += physical_coverage_change[position]
	current_sum += sum_insert_length_change[position]
	if current_coverage>0:
		average=current_sum/current_coverage
		f1.write(str(average) + '\n')
	else:
		f1.write('0\n')
f1.flush()
f1.closed

# seq coverage
f1 = open('seq_cov.wig', 'w')
f1.write("fixedStep chrom=" + chrom + " start=" + str(target_start) + " step=1 span=1\n")
current_coverage = 0
for position in range(target_end-target_start):
	current_coverage += sequence_coverage_change[position]
	f1.write(str(current_coverage) + '\n')
f1.flush()
f1.closed

# single mappers
f1 = open('singlemappers.wig', 'w')
f1.write("fixedStep chrom=" + chrom + " start=" + str(target_start) + " step=1 span=1\n")
current_coverage = 0
for position in range(target_end-target_start):
	current_coverage += single_mapper_change[position]
	f1.write(str(current_coverage) + '\n')
f1.flush()
f1.closed

# single mappers
f1 = open('10xlog2_ratio_for_rev.wig', 'w')
f1.write("fixedStep chrom=" + chrom + " start=" + str(target_start) + " step=1 span=1\n")
current_for = 0
current_rev = 0
for position in range(target_end-target_start):
	current_for += reads_for_change[position]
	current_rev += reads_rev_change[position]
	logratio=int(10*math.log2((current_for+0.1)/(current_rev+0.1)))
	f1.write(str(logratio) + '\n')
f1.flush()
f1.closed

# HS reads
f1 = open('reads_with_HS.wig', 'w')
f1.write("fixedStep chrom=" + chrom + " start=" + str(target_start) + " step=1 span=1\n")
current_coverage = 0
for position in range(target_end-target_start):
	current_coverage += reads_with_HS_change[position]
	f1.write(str(current_coverage) + '\n')
f1.flush()
f1.closed

# mates too long
f1 = open('mate_to_long.wig', 'w')
f1.write("fixedStep chrom=" + chrom + " start=" + str(target_start) + " step=1 span=1\n")
current_coverage = 0
for position in range(target_end-target_start):
	current_coverage += mate_to_long[position]
	f1.write(str(current_coverage) + '\n')
f1.flush()
f1.closed

# flag_256
f1 = open('flag_256.wig', 'w')
f1.write("fixedStep chrom=" + chrom + " start=" + str(target_start) + " step=1 span=1\n")
current_coverage = 0
for position in range(target_end-target_start):
	current_coverage += flag_256[position]
	f1.write(str(current_coverage) + '\n')
f1.flush()
f1.closed





