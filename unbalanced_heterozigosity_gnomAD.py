#!/usr/bin/python3

import sys
import re
from IPython import embed
import argparse
import pandas
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy import stats
import multiprocessing
import os
import time
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import gzip
import time

def get_header(vcf):
	tmp_list = []
	try:
		v=open(vcf, 'r')
		for line in v:
			if line.startswith('#'):
				tmp_list.append(line.rstrip())  
			else: 
				break
	except:
		v=gzip.open(vcf, 'rb')
		for line in v:
			line = line.decode('utf-8')
			if line.startswith('#'):
				tmp_list.append(line.rstrip())  
			else: 
				break
			
	v.close()
	return tmp_list	

def pseudo_normatization(info):
	chrom, pos, ref, alt = info
	# modification from right to left ATTT > ATT --> AT > A
	while len(ref)>1 and len(alt)>1 and ref[-1] == alt[-1]:
		ref = ref[:-1]
		alt = alt[:-1]
	# modification from left to right TGT > TGA --> T > A --> here position changes
	while len(ref)>1 and len(alt)>1 and ref[0] == alt[0]:
		ref = ref[1:]
		alt = alt[1:]
		pos += 1
	return chrom, pos, ref, alt 

def binomial_test(ntimes_freq, test):	#test_type should be greater [G] if we want one-tailed test or two-sided [T] for two-tailed test 
	#print(ntimes_freq, test)
	ntimes, freq = ntimes_freq
	test_type = test
	# one-tailed test vs two-sided 
	pvalue = stats.binom_test(ntimes, n=222, p=float(freq), alternative=test_type)
	#return pvalue if pvalue != 1 else numpy.nan 
	# freq == numpy.nan --> pvalue equal to nan if test_type = 'greater' 
	# freq == numpy.nan --> pvalue equal to 1.0 if test_type = 'two-tailed' 
	return pvalue #if freq != numpy.nan else numpy.nan

def binomial_test_mod(ntimes_sampcount_freq, test):	#test_type should be greater [G] if we want one-tailed test or two-sided [T] for two-tailed test 
	#print(ntimes_freq, test)
	ntimes, sample_count, freq = ntimes_sampcount_freq
	test_type = test
	# one-tailed test vs two-sided 
	pvalue = stats.binom_test(ntimes, n=sample_count, p=float(freq), alternative=test_type)
	#return pvalue if pvalue != 1 else numpy.nan 
	# freq == numpy.nan --> pvalue equal to nan if test_type = 'greater' 
	# freq == numpy.nan --> pvalue equal to 1.0 if test_type = 'two-tailed' 
	return pvalue #if freq != numpy.nan else numpy.nan

def split4parallelization(df, cpu):
	total_length = len(df)
	fraction = int(total_length/cpu)
	df_list = []
	for i in range(cpu):
		if i == cpu - 1:
			df_list.append(pandas.DataFrame(df.iloc[i*fraction:total_length]))
		else:
			df_list.append(pandas.DataFrame(df.iloc[i*fraction:(i+1)*fraction]))
	return df_list

def isNaN(num):
	return num != num

def apply_binomial_test_for_het(df):
	return df['het_count expected_het_frequency'.split()].apply(binomial_test, args=(['greater']), axis='columns')

def apply_binomial_test_for_het_mod(df):
	return df['GC_Het sample_count expected_het_frequency'.split()].apply(binomial_test_mod, args=(['greater']), axis='columns')

def apply_binomial_test_for_homo_alt(df):
	return df['homo_alt_count expected_homo_alt_frequency'.split()].apply(binomial_test, args=(['two-sided']), axis='columns')

def apply_binomial_test_for_homo_ref(df):
	return df['homo_ref_count expected_homo_ref_frequency'.split()].apply(binomial_test, args=(['two-sided']), axis='columns')

def apply_binomial_test(df, genotype):
	columns = [genotype + "_count", "expected_" + genotype + "_frequency"]
	return df[columns].apply(binomial_test, args=(['two-sided']), axis='columns')

def all_snp(allele_list):
	all_snp = True
	for allele in allele_list:
		if len(allele) > 1:
			return False
	return all_snp

def allele_count(samples):
	primary_allele_count = 0
	other_allele_count = 0
	ref_allele_count = 0
	tot_allele_count = len(samples)*2
	for sample in samples:
		if pandas.isnull(sample):
			ref_allele_count += 2
		elif sample == "./.":
			ref_allele_count += 2
		elif sample == '.' or sample == '0|0':				#0/0
			ref_allele_count += 2
		elif sample.startswith('0/1') or sample.startswith('0|1') or sample.startswith('1|0'):	#0/1
			primary_allele_count += 1
			ref_allele_count += 1
		elif sample.startswith('1/1') or sample.startswith('1|1'):	#1/1
			primary_allele_count += 2
		elif sample.startswith('0/') or sample.startswith('0|'):	#0/2, 0/3 ...
			ref_allele_count += 1
			other_allele_count += 1
		elif sample.startswith('1/') or sample.startswith('1|'):	#1/2, 1/3 ...
			primary_allele_count += 1
			other_allele_count += 1
		elif ('/1' in sample) or ('|1' in sample):			#2/1, 3/1 ...
			primary_allele_count += 1
			other_allele_count += 1
		elif '/0' in sample:			#hope never happens
			ref_allele_count += 1
			other_allele_count += 1
		else:							#2/2, 2/3, 3/3 ...
			other_allele_count += 2
	return ref_allele_count, ref_allele_count/tot_allele_count, primary_allele_count, primary_allele_count/tot_allele_count, other_allele_count, other_allele_count/tot_allele_count

def extract_allele_and_genotype_info(info):
	AC = AC_pattern.search(info).group(0)
	AC_AFR = AC_AFR_pattern.search(info).group(0)
	AC_AMR = AC_AMR_pattern.search(info).group(0)
	AC_ASJ = AC_ASJ_pattern.search(info).group(0)
	AC_EAS = AC_EAS_pattern.search(info).group(0)
	AC_FIN = AC_FIN_pattern.search(info).group(0)
	AC_NFE = AC_NFE_pattern.search(info).group(0)
	AC_OTH = AC_OTH_pattern.search(info).group(0)
	AN = AN_pattern.search(info).group(0)
	AN_AFR = AN_AFR_pattern.search(info).group(0)
	AN_AMR = AN_AMR_pattern.search(info).group(0)
	AN_ASJ = AN_ASJ_pattern.search(info).group(0)
	AN_EAS = AN_EAS_pattern.search(info).group(0)
	AN_FIN = AN_FIN_pattern.search(info).group(0)
	AN_NFE = AN_NFE_pattern.search(info).group(0)
	AN_OTH = AN_OTH_pattern.search(info).group(0)

	GC = GC_pattern.search(info).group(0)
	GC_Hom_Ref, GC_Het, GC_Hom_Alt = GC.split(',')

	GC_AFR = GC_AFR_pattern.search(info).group(0)
	GC_AFR_Hom_Ref, GC_AFR_Het, GC_AFR_Hom_Alt = GC_AFR.split(',')

	GC_AMR = GC_AMR_pattern.search(info).group(0)
	GC_AMR_Hom_Ref, GC_AMR_Het, GC_AMR_Hom_Alt = GC_AMR.split(',')

	GC_ASJ = GC_ASJ_pattern.search(info).group(0)
	GC_ASJ_Hom_Ref, GC_ASJ_Het, GC_ASJ_Hom_Alt = GC_ASJ.split(',')

	GC_EAS = GC_EAS_pattern.search(info).group(0)
	GC_EAS_Hom_Ref, GC_EAS_Het, GC_EAS_Hom_Alt = GC_EAS.split(',')

	GC_FIN = GC_FIN_pattern.search(info).group(0)
	GC_FIN_Hom_Ref, GC_FIN_Het, GC_FIN_Hom_Alt = GC_FIN.split(',')

	GC_NFE = GC_NFE_pattern.search(info).group(0)
	GC_NFE_Hom_Ref, GC_NFE_Het, GC_NFE_Hom_Alt = GC_NFE.split(',')

	GC_OTH = GC_OTH_pattern.search(info).group(0)
	GC_OTH_Hom_Ref, GC_OTH_Het, GC_OTH_Hom_Alt = GC_OTH.split(',')
	

	return {'AC':AC, 'AN':AN, 'GC_Hom_Ref':GC_Hom_Ref, 'GC_Het':GC_Het, 'GC_Hom_Alt':GC_Hom_Alt}
#	return {'AC':AC, 'AC_AFR':AC_AFR, 'AC_AMR':AC_AMR, 'AC_ASJ':AC_ASJ, 'AC_EAS':AC_EAS, 'AC_FIN':AC_FIN, 'AC_NFE':AC_NFE, 'AC_OTH':AC_OTH, 'AN':AN, 'AN_AFR':AN_AFR, 'AN_AMR':AN_AMR, 'AN_ASJ':AN_ASJ, 'AN_EAS':AN_EAS, 'AN_FIN':AN_FIN, 'AN_NFE':AN_NFE, 'AN_OTH':AN_OTH, 'GC_Hom_Ref':GC_Hom_Ref, 'GC_Het':GC_Het, 'GC_Hom_Alt':GC_Hom_Alt, 'GC_AFR_Hom_Ref':GC_AFR_Hom_Ref, 'GC_AFR_Het':GC_AFR_Het, 'GC_AFR_Hom_Alt':GC_AFR_Hom_Alt, 'GC_AMR_Hom_Ref':GC_AMR_Hom_Ref, 'GC_AMR_Het':GC_AMR_Het, 'GC_AMR_Hom_Alt':GC_AMR_Hom_Alt, 'GC_ASJ_Hom_Ref':GC_ASJ_Hom_Ref, 'GC_ASJ_Het':GC_ASJ_Het, 'GC_ASJ_Hom_Alt':GC_ASJ_Hom_Alt, 'GC_EAS_Hom_Ref':GC_EAS_Hom_Ref, 'GC_EAS_Het':GC_EAS_Het, 'GC_EAS_Hom_Alt':GC_EAS_Hom_Alt,  'GC_FIN_Hom_Ref':GC_FIN_Hom_Ref, 'GC_FIN_Het':GC_FIN_Het, 'GC_FIN_Hom_Alt':GC_FIN_Hom_Alt, 'GC_NFE_Hom_Ref':GC_NFE_Hom_Ref, 'GC_NFE_Het':GC_NFE_Het, 'GC_NFE_Hom_Alt':GC_NFE_Hom_Alt, 'GC_OTH_Hom_Ref':GC_OTH_Hom_Ref, 'GC_OTH_Het':GC_OTH_Het, 'GC_OTH_Hom_Alt':GC_OTH_Hom_Alt}      

def annotate_block(df):
	res = df['INFO'].apply(extract_allele_and_genotype_info)	# return a pandas series containing in each row a dictionary with all values
	return pandas.DataFrame.from_records(res.reset_index()['INFO'], index=res.index)	# return a pandas DF with 40 columns each corresponding to a particular AC, AN or GC

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="A smart tool for counting alternative allele within [multisample] vcf file")
	parser.add_argument("-i", "--input-file", metavar="VCF", required=True, help="VCF file.")
	parser.add_argument("-c", "--cpu", metavar="CPU", type=int, help="Number of CPUs to be used [Default 1, All -1].")
	#parser.add_argument("-o", "--output-file", metavar="OUT", required=True, help="File to which write the TSV output.")
	args = parser.parse_args()
	if args.cpu == -1:
		cpu_num = os.cpu_count()
	else:
		cpu_num = args.cpu
	# patterns for alternative allele count
	AC_pattern = re.compile('(?<=AC=)[0-9]+')	# alternative allele count
	AC_AFR_pattern = re.compile('(?<=AC_AFR=)[0-9]+')
	AC_AMR_pattern = re.compile('(?<=AC_AMR=)[0-9]+')
	AC_ASJ_pattern = re.compile('(?<=AC_ASJ=)[0-9]+')
	AC_EAS_pattern = re.compile('(?<=AC_EAS=)[0-9]+')
	AC_FIN_pattern = re.compile('(?<=AC_FIN=)[0-9]+')
	AC_NFE_pattern = re.compile('(?<=AC_NFE=)[0-9]+')
	AC_OTH_pattern = re.compile('(?<=AC_OTH=)[0-9]+')
	# patterns for allele count
	AN_pattern = re.compile('(?<=AN=)[0-9]+')	# allele number
	AN_AFR_pattern = re.compile('(?<=AN_AFR=)[0-9]+')
	AN_AMR_pattern = re.compile('(?<=AN_AMR=)[0-9]+')
	AN_ASJ_pattern = re.compile('(?<=AN_ASJ=)[0-9]+')
	AN_EAS_pattern = re.compile('(?<=AN_EAS=)[0-9]+')
	AN_FIN_pattern = re.compile('(?<=AN_FIN=)[0-9]+')
	AN_NFE_pattern = re.compile('(?<=AN_NFE=)[0-9]+')
	AN_OTH_pattern = re.compile('(?<=AN_OTH=)[0-9]+')
	# patterns for genotype count
	GC_pattern = re.compile('(?<=GC=)([0-9]+),([0-9]+),([0-9]+)')
	GC_AFR_pattern = re.compile('(?<=GC_AFR=)([0-9]+),([0-9]+),([0-9]+)')
	GC_AMR_pattern = re.compile('(?<=GC_AMR=)([0-9]+),([0-9]+),([0-9]+)')
	GC_ASJ_pattern = re.compile('(?<=GC_ASJ=)([0-9]+),([0-9]+),([0-9]+)')
	GC_EAS_pattern = re.compile('(?<=GC_EAS=)([0-9]+),([0-9]+),([0-9]+)')
	GC_FIN_pattern = re.compile('(?<=GC_FIN=)([0-9]+),([0-9]+),([0-9]+)')
	GC_NFE_pattern = re.compile('(?<=GC_NFE=)([0-9]+),([0-9]+),([0-9]+)')
	GC_OTH_pattern = re.compile('(?<=GC_OTH=)([0-9]+),([0-9]+),([0-9]+)')
	
	final_df = pandas.DataFrame()
	header = get_header(args.input_file)
	print('INFO: Header length is {}'.format(len(header)))
	file_is_finished = False
	chunck_number = 1
	chunk_size = 200000
	header_len = len(header)
	columns = header[-1].split('\t')
	starting_time = time.time()
	total_rows_processed = 0
	while not file_is_finished:
		rows_to_skip = header_len + (chunck_number - 1)*chunk_size
		current_df = pandas.read_table(args.input_file, header=None, names=columns, skiprows=rows_to_skip, nrows=chunk_size)
		biallelic_df = current_df[~current_df.ALT.str.contains(',')].reset_index()
		biallelic_df['#CHROM POS REF ALT'.split()] = biallelic_df['#CHROM POS REF ALT'.split()].apply(pseudo_normatization, axis='columns').tolist()	#maybe these variants are already normalized
		filtered_biallelic_df = biallelic_df.loc[biallelic_df['#CHROM POS REF ALT'.split()].drop_duplicates().index].reset_index()[columns]
		current_blocks = split4parallelization(filtered_biallelic_df, cpu_num)
		annotation_start = time.time()
		with multiprocessing.Pool(cpu_num) as pool:
			results = pool.map(annotate_block, current_blocks)	# here df is passed to annotate_block function
			final_df = final_df.append(pandas.concat([filtered_biallelic_df['#CHROM POS REF ALT'.split()], pandas.concat(results)], axis='columns'), ignore_index=True)
		annotation_time = '{0:0.3f}'.format(time.time() - annotation_start)
		rows_x_sec = str(int(len(current_df) / (time.time() - annotation_start)))
		total_rows_processed += len(current_df)
		total_time =  '{0:0.3f}'.format(time.time() - starting_time)
		total_rows_x_sec = str(int(total_rows_processed/(time.time() - starting_time)))
		print('INFO: Block of {} rows processed in {}s [{} rows/sec]. Processed {} rows in a total time of {}s [{} rows/sec].'.format(len(current_df),  annotation_time, rows_x_sec, total_rows_processed, total_time, total_rows_x_sec), file=sys.stderr)	
		if len(current_df) < chunk_size:
			file_is_finished = True
		chunck_number += 1
	
		
	for field in 'AC AN GC_Het GC_Hom_Alt GC_Hom_Ref'.split():
		final_df[field] = final_df[field].apply(int)
	final_df['sample_count'] = final_df['GC_Het'] + final_df['GC_Hom_Alt'] + final_df['GC_Hom_Ref']
	final_df = final_df[final_df.sample_count > 0]

	final_df['alt_allele_frequency'] = final_df['AC']/final_df['AN']
	final_df['ref_allele_frequency'] = 1 - final_df['alt_allele_frequency'] 

	# expected genotype frequencies starting from allele freq
	final_df['expected_het_frequency'] = 2 * final_df['ref_allele_frequency'] * final_df['alt_allele_frequency']	# 2pq
	final_df['expected_homo_alt_frequency'] = final_df['alt_allele_frequency']**2	# q2	
	final_df['expected_homo_ref_frequency'] = final_df['ref_allele_frequency']**2	# p2

	# real genotype frequencies derived from our data
	final_df['het_frequency'] = final_df['GC_Het']/final_df['sample_count']
	final_df['homo_ref_frequency'] = final_df['GC_Hom_Ref']/final_df['sample_count']
	final_df['homo_alt_frequency'] = final_df['GC_Hom_Alt']/final_df['sample_count']

	#embed()
	# pvalue computation for het genotype using one-sided binomial test --> greater than
	df_blocks = split4parallelization(final_df, cpu_num)
	print('\rINFO: Binomial test on Heterozygous genotype...', end='', file=sys.stderr)
	starting_binomial_time = time.time()
	try:
		with multiprocessing.Pool(cpu_num) as pool:
			results_het = pool.map(apply_binomial_test_for_het_mod, df_blocks)
			final_df['het_one_tail_pvalue'] = pandas.concat(results_het)
	except:
		embed()
	final_binomial_time = '{0:0.3f}'.format(time.time() - starting_binomial_time)
	print('\rINFO: Binomial tests on Heterozygous genotype... Done ({}s).'.format(final_binomial_time), file=sys.stderr)

	# Benjamini/Hochberg FDR pvalue correction
	starting_correction_time = time.time()
	print("\rINFO: P-value correction for multiple tests using Benjamini/Hochberg FDR...", end='', file=sys.stderr)
	#for column in 'het homo_alt homo_ref het_one_tail'.split():
	column = 'het_one_tail'
	final_df[column+'_corrected_pvalue'] = multipletests(final_df[column+'_pvalue'], alpha=0.01, method='fdr_bh')[1] # 0 --> array with True and False, 1 --> array with corrected pvalue
	final_correction_time = '{0:0.3f}'.format(time.time() - starting_correction_time)
	print("\rINFO: P-value correction for multiple tests using Benjamini/Hochberg FDR...  Done ({}s).".format(final_correction_time), file=sys.stderr)

	embed()
	sys.exit()		




	args = parser.parse_args()
	allele_count_dict = {}
	allele_row_number = {}
	all_rows = []
	line_counter = 0
	allele_number = 0
	sample_number = 0
	variant_number = 0
	homo_ref_dict = {}
	homo_alt_dict = {}
	het_dict = {}
	multi_alts = {}
	multi_rows = {}
	multi_ref = {}
	header = []
	starting_reading_time = time.time()
	with open(args.input_file, 'r') as VCF:
		# VCF reading and computing allele count
		for line in VCF:
			line_counter+=1
			print('\r{} vcf lines processed'.format(line_counter), end='', file=sys.stderr)
			if not line.startswith('#'):
				#all_rows.append([])
				line = line.rstrip('\n')
				columns = line.split('\t')
				sample_number = len(columns[9:])
				allele_number = len(columns[9:]) * 2
				ALTS = columns[4].split(',')
				#is_multi_rows_variant = False
				#multi_alts = True if len(ALTS)>1 else False
				for i, alt in enumerate(ALTS):
					# pseudo normalization of variants
					POS, REF, ALT = pseudo_normatization(int(columns[1]), columns[3], alt)
					#key = '|'.join(columns[:2] + [columns[3]] + [alt])
					key = '|'.join([columns[0], str(POS), REF, ALT])
					key_ref = '|'.join([columns[0], str(POS), REF,''])	# added | at the end to have an equal sorting between multi_alts_keys and multi_ref_keys
					#all_rows[variant_number].append(key)
					if not key_ref in multi_ref:
						multi_ref[key_ref] = 1
					else:
						multi_ref[key_ref] += 1
					if not key in allele_count_dict:
						allele_count_dict[key] = 0
						homo_ref_dict[key] = 0
						homo_alt_dict[key] = 0
						het_dict[key] = 0
						multi_alts[key] = len(ALTS)>1 
						multi_rows[key] = False
					elif key in allele_count_dict:
						if not multi_alts[key]:	# multi_alts[key] == False --> if once is multi_alt the variant is always multi alt
							multi_alts[key] = len(ALTS)>1
						multi_rows[key] = True
						is_multi_rows_variant = True		
					for sample in columns[9:]:
						if sample.startswith('./.'):
							homo_ref_dict[key] += 1 
							#continue
						elif sample.startswith('{}/{}'.format(i+1,i+1)): 	# es 1/1, 2/2
							allele_count_dict[key] += 2
							homo_alt_dict[key] += 1 
						elif sample.startswith('{}/'.format(i+1)): 		# es 1/2, 1/3, 1/0
							allele_count_dict[key] += 1				
							het_dict[key] += 1 
						elif sample[2] == str(i+1):						# es 0/1, 3/1
							allele_count_dict[key] += 1
							het_dict[key] += 1 
					#	old_key = key
					variant_number += 1 
			else:
				header.append(line.rstrip('\n'))

		final_reading_time = '{0:0.3f}'.format(time.time() - starting_reading_time)
		print('\r{} vcf lines processed in {}s.'.format(line_counter, final_reading_time), file=sys.stderr)
		print('\rVery long cycle due to very ugly programming...', file=sys.stderr, end='')
		i=0
		multi_alts_sorted_keys = sorted(multi_alts)
		multi_ref_sorted_keys = sorted(multi_ref)

#		for key_ref in multi_ref_sorted_keys:
#			if multi_ref[key_ref] > 1:
#				for key in multi_alts:
#					if key.startswith(key_ref+'|'):
#						multi_alts[key] = True
#
		for key_ref in multi_ref_sorted_keys:
			if multi_ref[key_ref] > 1:
				is_explored = False
				for key in multi_alts_sorted_keys[i:]:
					if key.startswith(key_ref):
						multi_alts[key] = True
						is_explored = True
					elif is_explored:
						break
					i+=1
		print('\rVery long cycle due to very ugly programming... Done', file=sys.stderr)

		header_fields = header[-1].split('\t')
		second_ALT_index = len(header_fields) - header_fields[::-1].index('ALT') - 1
		header_fields[second_ALT_index] = 'ALT_sample'

		vcf = pandas.read_table(args.input_file, comment='#', header=None, names=header_fields)
		filtered_vcf = vcf[vcf.ALT.str.contains(',')].copy()
		multiallelic_snp = filtered_vcf[filtered_vcf.ALT.str.split(',').apply(all_snp)].copy()
		triallelic_snp = multiallelic_snp[multiallelic_snp.ALT.str.split(',').apply(len) == 2].copy().reset_index(drop=True)
		quadriallelic_snp = multiallelic_snp[multiallelic_snp.ALT.str.split(',').apply(len) == 3].copy()
		allele_counts = pandas.DataFrame(triallelic_snp[triallelic_snp.columns[9:]].apply(allele_count, axis='columns').tolist(), columns='ref_allele_count ref_allele_freq primary_AA_count primary_AA_freq secondary_AA_count secondary_AA_freq'.split())
		for column in allele_counts.columns:
			triallelic_snp[column] = allele_counts[column]

		# building dataframe and working on it
		raw_df = pandas.DataFrame.from_dict(allele_count_dict, orient='index')
		het_df = pandas.DataFrame.from_dict(het_dict, orient='index')
		homo_ref_df = pandas.DataFrame.from_dict(homo_ref_dict, orient='index')
		homo_alt_df = pandas.DataFrame.from_dict(homo_alt_dict, orient='index')
		multi_alts_df = pandas.DataFrame.from_dict(multi_alts, orient='index') 
		multi_rows_df = pandas.DataFrame.from_dict(multi_rows, orient='index') 
		raw_df['alt_allele_count'] = raw_df[0]
		raw_df['homo_ref_count'] = homo_ref_df[0]
		raw_df['homo_alt_count'] = homo_alt_df[0]
		raw_df['het_count'] = het_df[0]
		raw_df['multi_alts'] = multi_alts_df[0]
		raw_df['multi_rows'] = multi_rows_df[0]
		del raw_df[0]
		new_df = raw_df.reset_index()
		splitted_index_df = new_df['index'].str.split('|', expand=True)
		for i, element in enumerate('chrom pos ref alt'.split()):
			new_df[element] = splitted_index_df[i]
		del new_df['index']
		new_df['chrom'] =  new_df['chrom'].str[3:]
		new_df['pos'] = new_df['pos'].apply(int)
		x_chrom = new_df['chrom'] == 'X'
		y_chrom = new_df['chrom'] == 'Y'
		new_df.loc[x_chrom, 'chrom'] = '23'
		new_df.loc[y_chrom, 'chrom'] = '24'
		new_df['chrom'] = new_df['chrom'].apply(int)
		final_df = new_df.sort_values(by='chrom pos ref alt'.split())
		final_df.loc[x_chrom,'chrom'] = 'X'
		final_df.loc[y_chrom, 'chrom'] = 'Y'
		final_df.reset_index(drop=True, inplace=True)
		final_df['chrom'] = 'chr'+final_df['chrom'].apply(str)
		final_df['NEXOMES'] = final_df['homo_alt_count'] + final_df['het_count']
		final_df['chrom pos ref alt alt_allele_count multi_alts multi_rows NEXOMES'.split()].to_csv(args.output_file+'.ALL_NORMALIZED_VARIANTS.tsv', sep='\t', index=None, header='#CHROM POS REF ALT NTIMES multi_alts multi_rows NEXOMES'.split())

		mono_allelic_df = final_df.loc[(~final_df.multi_alts) & (~final_df.multi_rows)]
		problematic_allele_df = final_df.loc[(final_df.multi_alts) | (final_df.multi_rows)]
		problematic_allele_df['chrom pos ref alt alt_allele_count multi_alts multi_rows'.split()].to_csv(args.output_file+'.PROBLEMATIC_NORMALIZED_VARIANTS.tsv', sep='\t', index=None, header='#CHROM POS REF ALT NTIMES multi_alts multi_rows'.split())
		
		#embed()
		final_df = mono_allelic_df.copy()
		final_df.reset_index(drop=True, inplace=True)
		final_df['alt_allele_frequency'] = final_df['alt_allele_count']/allele_number	# q
		final_df['ref_allele_frequency'] = 1 - final_df['alt_allele_frequency'] 	# p

		# expected genotype frequencies starting from allele freq
		final_df['expected_het_frequency'] = 2 * final_df['ref_allele_frequency'] * final_df['alt_allele_frequency']	# 2pq
		final_df['expected_homo_alt_frequency'] = final_df['alt_allele_frequency']**2	# q2	
		final_df['expected_homo_ref_frequency'] = final_df['ref_allele_frequency']**2	# p2

		# real genotype frequencies derived from our data
		final_df['het_frequency'] = final_df['het_count']/sample_number
		final_df['homo_ref_frequency'] = final_df['homo_ref_count']/sample_number
		final_df['homo_alt_frequency'] = final_df['homo_alt_count']/sample_number

		# difference between real and expected freq
		for genotype in 'het homo_alt homo_ref'.split():
			final_df[genotype + '_difference'] = final_df[genotype + '_frequency'] - final_df['expected_' + genotype + '_frequency']

		df_blocks = split4parallelization(final_df, os.cpu_count())
		# pvalue computation for each genotype using two-sided binomial test
		for genotype in 'het homo_alt homo_ref'.split():
			print('\rBinomial test on {} genotype...'.format(genotype), end='', file=sys.stderr)
			starting_binomial_time = time.time()
			with multiprocessing.Pool(os.cpu_count()) as pool:
				new_blocks = [(block, genotype) for block in df_blocks]
				results = pool.starmap(apply_binomial_test, new_blocks)
				final_df[genotype+"_pvalue"] = pandas.concat(results)
				#results_het = pool.map(apply_binomial_test_for_het, df_blocks)
				#final_df['het_pvalue'] = pandas.concat(results_het)
				#results_homo_alt = pool.map(apply_binomial_test_for_homo_alt, df_blocks)
				#final_df['homo_alt_pvalue'] = pandas.concat(results_homo_alt)
				#results_homo_ref = pool.map(apply_binomial_test_for_homo_ref, df_blocks)
				#final_df['homo_ref_pvalue'] = pandas.concat(results_homo_ref)
			final_binomial_time = '{0:0.3f}'.format(time.time() - starting_binomial_time)
			print('\rBinomial tests on {} genotype... Done ({}s).'.format(genotype, final_binomial_time), file=sys.stderr)

		# pvalue computation for het genotype using one-sided binomial test --> greater than
		with multiprocessing.Pool(os.cpu_count()) as pool:
			results_het = pool.map(apply_binomial_test_for_het, df_blocks)
			final_df['het_one_tail_pvalue'] = pandas.concat(results_het)
		
		# Benjamini/Hochberg FDR pvalue correction
		starting_correction_time = time.time()
		print("\rP-value correction for multiple tests using Benjamini/Hochberg FDR...", end='', file=sys.stderr)
		for column in 'het homo_alt homo_ref het_one_tail'.split():
			final_df[column+'_corrected_pvalue'] = pandas.Series(multipletests(final_df[column+'_pvalue'], alpha=0.01, method='fdr_bh')[1]) # 0 --> array with True and False, 1 --> array with corrected pvalue
		final_correction_time = '{0:0.3f}'.format(time.time() - starting_correction_time)
		print("\rP-value correction for multiple tests using Benjamini/Hochberg FDR...  Done ({}s).".format(final_correction_time), file=sys.stderr)

		sorted_columns = 'chrom pos ref alt alt_allele_count alt_allele_frequency het_count het_frequency expected_het_frequency het_pvalue het_corrected_pvalue het_one_tail_pvalue het_one_tail_corrected_pvalue homo_alt_count homo_alt_frequency expected_homo_alt_frequency homo_alt_pvalue homo_alt_corrected_pvalue homo_ref_count homo_ref_frequency expected_homo_ref_frequency homo_ref_pvalue homo_ref_corrected_pvalue'.split()
		final_df[sorted_columns].to_csv(args.output_file+'.MONOALLELIC_NORMALIZED_VARIANTS.tsv', sep='\t', index=None, header='#CHROM POS REF ALT'.split()+sorted_columns[4:])
		# computing significant variants
		significant_at_least_for_one = final_df.loc[(final_df.het_corrected_pvalue<0.01) | (final_df.homo_ref_corrected_pvalue<0.01) | (final_df.homo_alt_corrected_pvalue<0.01)].copy()
		significant_for_all_df = final_df.loc[(final_df.het_corrected_pvalue<0.01) & (final_df.homo_ref_corrected_pvalue<0.01) & (final_df.homo_alt_corrected_pvalue<0.01)].copy()
		significant_for_het = final_df.loc[final_df.het_one_tail_corrected_pvalue<0.01].copy()
		#significant_df[(significant_df.het_freq> significant_df.expected_het_frequency)]['het_freq expected_het_frequency'.split()]
		#significant_at_least_for_one.chrom = 'chr' + significant_at_least_for_one['chrom'].apply(str)
		significant_at_least_for_one[sorted_columns].to_csv(args.output_file+'.UNBALANCED_GENOTYPE.tsv', sep='\t', index=None, header='#CHROM POS REF ALT'.split()+sorted_columns[4:])

		#significant_for_het.chrom = 'chr' + significant_for_het['chrom'].apply(str)
		significant_for_het[sorted_columns].to_csv(args.output_file+'.UNBALANCED_HET_GENOTYPE.tsv', sep='\t', index=None, header='#CHROM POS REF ALT'.split()+sorted_columns[4:])
		sys.exit()

		starting_plotting_time = time.time()
		print("\rPlotting figures...", end='', file=sys.stderr)
		# pvalue figures
		het_pvalue_kde = final_df['het_pvalue'].plot(kind='kde', xlim=(0,1))
		het_pvalue_kde.get_figure().savefig('het_pvalue_kde.png')
		plt.clf()
		print("\rPlotting figures... 1/7 Done.", end='', file=sys.stderr)
		homo_alt_pvalue_kde = final_df['homo_alt_pvalue'].plot(kind='kde', xlim=(0,1))
		homo_alt_pvalue_kde.get_figure().savefig('homo_alt_pvalue_kde.png')
		plt.clf()
		print("\rPlotting figures... 2/7 Done.", end='', file=sys.stderr)
		homo_ref_pvalue_kde = final_df['homo_ref_pvalue'].plot(kind='kde', xlim=(0,1))
		homo_ref_pvalue_kde.get_figure().savefig('homo_ref_pvalue_kde.png')
		plt.clf()
		all_pvalue_kde = final_df['homo_ref_pvalue homo_alt_pvalue het_pvalue'.split()].plot(kind='kde', xlim=(0,1), legend=True)
		all_pvalue_kde.get_figure().savefig('all_pvalue.png')
		plt.clf()
		print("\rPlotting figures... 3/7 Done.", end='', file=sys.stderr)
		min_interval = 0.02
		het_pvalue_hist = final_df['het_pvalue'].plot(kind='hist', xlim=(0,1), bins=[i*min_interval for i in range(int(1/min_interval)+1)])
		het_pvalue_hist.get_figure().savefig('het_pvalue_hist.png')

		# difference figures
		diff_het = final_df['het_difference'].plot(kind='kde', xlim=[-0.2, 0.2])
		diff_het.get_figure().savefig('het_difference_kde.png')
		plt.clf()
		print("\rPlotting figures... 4/7 Done.", end='', file=sys.stderr)
		diff_homo_ref = final_df['homo_ref_difference'].plot(kind='kde', xlim=[-0.2, 0.2])
		diff_homo_ref.get_figure().savefig('homo_ref_difference_kde.png')
		plt.clf()
		print("\rPlotting figures... 5/7 Done.", end='', file=sys.stderr)
		diff_homo_alt = final_df['homo_alt_difference'].plot(kind='kde', xlim=[-0.2, 0.2])
		diff_homo_alt.get_figure().savefig('homo_alt_difference_kde.png')
		plt.clf()
		print("\rPlotting figures... 6/7 Done.", end='', file=sys.stderr)
		all_diff = final_df['homo_ref_difference homo_alt_difference het_difference'.split()].plot(kind='kde', xlim=[-0.2, 0.2], legend=True)
		all_diff.get_figure().savefig('all_difference_kde.png')
		final_correction_time = '{0:0.3f}'.format(time.time() - starting_plotting_time)
		print("\rPlotting figures... 7/7 Done ({}s).".format(final_correction_time), file=sys.stderr)

		sys.exit()
		embed()
		
		#final_df['het_pvalue'] = final_df['het_count expected_het_frequency'.split()].apply(binomial_test, args=(['two-sided']), axis='columns')
		#final_df['homo_alt_pvalue'] = final_df['homo_alt_count expected_homo_alt_frequency'.split()].apply(binomial_test, args=(['two-sided']), axis='columns')
		final_df['chrom'] = 'chr' + final_df['chrom'].apply(str)
		final_df['chrom pos ref alt alt_allele_count alt_allele_frequency'.split()].to_csv(args.output_file, sep='\t', index=None, header='#CHROM POS REF ALT ALLELE_COUNT ALLELE_FREQ'.split())  
		

		# more than 1 alt allele
		#final_df.set_index('chrom pos ref'.split()).loc[final_df.groupby('chrom pos ref'.split())['homo_alt_pvalue'].count() > 1]
