#/usr/bin python3

import sys
import re
import argparse
from IPython import embed
import pandas

wholeAAPatternExtended = re.compile('p.\w{3}\d+\w{3}')
AAandPOSpattern = re.compile('p.(\w{3})(\d+)(\w{3})')
onlyAApattern = re.compile('p.(\w{3})\d+(\w{3})')
three_letters_code = {'G':'Gly', 'A':'Ala', 'V':'Val', 'L':'Leu', 'I':'Ile', 'M':'Met', 'F':'Phe', 'W':'Trp', 'P':'Pro', 'S':'Ser', 'T':'Thr', 'C':'Cys', 'Y':'Tyr', 'N':'Asn', 'Q':'Gln', 'D':'Asp', 'E':'Glu', 'K':'Lys', 'R':'Arg', 'H':'His', '?':'?', '*':'Ter'}

def parse_uniprot(uniprot_file):
	'''Returns three dictionary derived from parsing human protein variation of UniProt:
		1) gene_name|p_hgvs as key and type of variant as value
		2) rs_code|p_hgvs as key and type of variant as value
		3) rs_code as key and a list with p_hgvs and type of variant as value'''
	with open(uniprot_file, 'r') as unifile:
		need_to_read = False
		uniprot_dict = {}
		dbsnp_uniprot_dict = {}
		dbsnp_codes = {}
		for line in unifile:
			fields = line.split()
			if len(fields) < 3 or not fields[2].startswith('VAR'):
				continue
			uniprot_dict[fields[0]+"|"+fields[3]] = fields[4]
			dbsnp_uniprot_dict[fields[5]+"|"+fields[3]] = fields[4]
			dbsnp_codes[fields[5]] = [fields[3], fields[4]]
	unifile.close()
	return uniprot_dict, dbsnp_uniprot_dict, dbsnp_codes 

def load_uniprot_isoforms(sequence_file):
	'''Returns a pandas dataframe containing two columns: protein_name and sequence, indexed by gene_name'''
	df = pandas.read_table(sequence_file, sep='\t', names='protein_name gene_name sequence'.split())
	df_filtered = df[~df['protein_name'].isnull()].copy()	# remove rows without gene and protein name
	return df_filtered.set_index('gene_name')

def load_uniprot_sequence_length(sequence_file):
	'''Returns a pandas dataframe containing three columns: protein_name, protein_length and sequence, indexed by gene_name'''
	first_df = pandas.read_table(sequence_file, sep='\t')
	df = first_df[['Entry','Gene names', 'Length', 'Sequence']].copy()
	df.rename(index=str, columns={'Entry':'protein_name','Gene names':'gene_name', 'Length':'length', 'Sequence':'sequence'}, inplace=True)
	df_filtered = df[~df['gene_name'].isnull()].copy()
	return df_filtered.set_index('gene_name')




def get_protein_predictions(line, va):
	'''Returns the protein prediction gene_name|p_hgvs starting from the whole line of vcf and considering only missense snp'''
	fields = line.split('\t')
	index_correction = 0 if va == 'snpeff' else 1
	annotations = fields[7].split(';')[-1]
	protein_predictions = []
	protein_length = []
	if len(fields[3]) == len(fields[4]) and ('missense' in annotations or 'stop_gain' in annotations or 'stop_lost' in annotations):	# get hgvs protein code only for missense or nonsense SNP variants
		for annotation in annotations.split(','):
			columns = annotation.split('|')
			try:
				protein_predictions.append(columns[3]+"|"+re.findall(wholeAAPatternExtended, columns[10+index_correction])[0]) 
			except:
				protein_predictions.append(columns[3]+"|")
			protein_length.append(columns[13+index_correction].split('/')[-1])
	return protein_predictions, protein_length

def get_dbsnp_prediction(rs_code, p):
	'''Return a string with rs_code joined to hgvs protein code'''
	pred_elements = p.split('|')
	try:
		protein_pred = pred_elements[1]
	except:
		protein_pred = ''
	return rs_code + '|' + protein_pred

def inverted_prediction(p):
	'''Returns a string with gene_name joined to hgvs protein code with inverted aa (p.AA2posAA1)'''
	pred_elements = p.split('|')
	gene_name = pred_elements[0]
	try:
		first_aa, pos, second_aa = re.findall(AAandPOSpattern, pred_elements[1])[0]
		return gene_name + '|p.' + second_aa + pos + first_aa
	except:
		return gene_name + '|'

def checkAAidentity(rs_code, dbsnp_codes, p):
	''' Returns True if the couple of AAs is the same in the two pattern (AA1s1 = AA1s2 and AA2s1 = AA2s2) independently from the position on the transcripts.'''
	pred_elements = p.split('|')
	try:
#		if re.findall(onlyAApattern, dbsnp_codes[rs_code][0])[0] == re.findall(onlyAApattern, pred_elements[1])[0]:
#			print(re.findall(onlyAApattern, dbsnp_codes[rs_code][0])[0], re.findall(onlyAApattern, pred_elements[1])[0])
		return re.findall(onlyAApattern, dbsnp_codes[rs_code][0])[0] == re.findall(onlyAApattern, pred_elements[1])[0]
	except:
		return False

def checkAApos(dbsnp_p, p):
	dbsnp_AApos = dbsnp_p[:-3]
	try:
		ann_AApos = p.split('|')[1][:-3]
	except:
		ann_AApos = ''
	return dbsnp_AApos == ann_AApos

def checkAAintoSeq(seq, POS, SECOND_AA, protein_length):
	if len(seq) >= int(POS):
		return (three_letters_code[seq[int(POS)-1]] == SECOND_AA) and (len(seq)==int(protein_length))
	else:
		return False

def analyze_snpeff_vcf(vcf_input_file, vcf_output, uniprot_dict, dbsnp_dict, dbsnp_codes, uniprot_seq):
	'''Function working on SnpEff annotated file to establish if protein variation are present within UniProt curated DB'''
	uniprot_seq_index = uniprot_seq.index
	with open(vcf_input_file, 'r') as vcf_file:
		tot_counter = 0
		unsolved = 0
		completely_solved = 0
		probable_solved = 0
		possible_solved = 0
		variant_in_reference = 0 
		for line in vcf_file:
			line = line.rstrip('\n')
			if line.startswith('#'):
				print(line, file=vcf_output)
				continue
			#print(line)
			else:
				columns = line.split('\t')
				rs_code = columns[2]
				tot_counter += 1 
				protein_predictions, protein_length = get_protein_predictions(line, 'snpeff')
				if protein_predictions:
					c = 1
					for i, prediction in enumerate(protein_predictions):
						gene_name, hgvs_protein = prediction.split('|')
						if hgvs_protein == '':	#variant with no protein prediction
							if c == len(protein_predictions):
								unsolved += 1
								print('\t'.join([line, 'not_found']), file=vcf_output)
							c+=1
							continue
						#print(prediction,hgvs_protein)
						first_aa, pos, second_aa = re.findall(AAandPOSpattern, hgvs_protein)[0]
						# prediction --> GENE_NAME | p.AA1posAA2
						if prediction in uniprot_dict:
							print('\t'.join([line, uniprot_dict[prediction]]), file=vcf_output)
							completely_solved += 1
							break
						elif get_dbsnp_prediction(rs_code, prediction) in dbsnp_dict:
							print('\t'.join([line, dbsnp_dict[get_dbsnp_prediction(rs_code, prediction)]]), file=vcf_output)
							completely_solved += 1
							break
						elif inverted_prediction(prediction) in uniprot_dict:
							print('\t'.join([line, 'variant_is_in_reference']), file=vcf_output)
							variant_in_reference += 1
							#completely_solved += 1
							break
						elif rs_code in dbsnp_codes and checkAAidentity(rs_code, dbsnp_codes, prediction):
							print('\t'.join([line, 'probable_' + dbsnp_codes[rs_code][-1]]),  file=vcf_output)
							probable_solved += 1
							break
						elif rs_code in dbsnp_codes and checkAApos(dbsnp_codes[rs_code][0], prediction): 
							print('\t'.join([line, 'possible_' + dbsnp_codes[rs_code][-1]]),  file=vcf_output)
							possible_solved += 1
							break
						elif gene_name in uniprot_seq_index	and type(uniprot_seq.loc[gene_name,'sequence'])!=str and args.isoforms:
							#if uniprot_seq.loc[gene_name].sequence.apply(lambda x:three_letters_code[x[int(pos)-1]]==second_aa).any():
							#if uniprot_seq.loc[gene_name].sequence.apply(checkAAintoSeq, args=(pos, second_aa)).any():
							if uniprot_seq.loc[gene_name].sequence.apply(checkAAintoSeq, args=(pos, second_aa, protein_length[i])).any() and first_aa!=second_aa:
									#print(prediction, uniprot_seq.loc[gene_name, 'sequence'].apply(len),protein_length[i])
									print('\t'.join([line, 'variant_is_in_reference']), file=vcf_output)
									variant_in_reference += 1
									#completely_solved += 1
									#probable_solved += 1
									break
							elif c == len(protein_predictions):
								unsolved += 1
								print('\t'.join([line, 'not_found']), file=vcf_output)
						elif gene_name in uniprot_seq_index	and type(uniprot_seq.loc[gene_name,'sequence'])==str and args.isoforms:
							#if checkAAintoSeq(uniprot_seq.loc[gene_name].sequence, pos, second_aa):
							if checkAAintoSeq(uniprot_seq.loc[gene_name].sequence, pos, second_aa,  protein_length[i]) and first_aa!=second_aa:
									#print(uniprot_seq.loc[gene_name, 'sequence'], prediction)
									#print(prediction, len(uniprot_seq.loc[gene_name, 'sequence']), protein_length[i])
									print('\t'.join([line, 'variant_is_in_reference']), file=vcf_output)
									#completely_solved += 1
									variant_in_reference += 1
									#probable_solved += 1
									break
							elif c == len(protein_predictions):
								unsolved += 1
								print('\t'.join([line, 'not_found']), file=vcf_output)
#						elif gene_name in uniprot_seq_index and args.sequence_length and type(uniprot_seq.loc[gene_name,'sequence'])==str:  
#							#print(prediction, len(uniprot_seq.loc[gene_name, 'sequence']), protein_length[i])
#							if checkAAintoSeq(uniprot_seq.loc[gene_name].sequence, pos, second_aa, protein_length[i]) and first_aa!=second_aa:# and int(uniprot_seq.loc[gene_name, 'length']) == int(protein_length[i]):
#									print(prediction, len(uniprot_seq.loc[gene_name, 'sequence']), protein_length[i])
#									#print(uniprot_seq.loc[gene_name, 'sequence'], prediction)
#									print('\t'.join([line, 'variant_is_in_reference']), file=vcf_output)
#									#completely_solved += 1
#									variant_in_reference += 1
#									#probable_solved += 1
#									break
#							elif c == len(protein_predictions):
#								unsolved += 1
#								print('\t'.join([line, 'not_found']), file=vcf_output)
#							
								

						elif c == len(protein_predictions):	# if the script comes inside here it means that SnpEff prediction is not inside common UniProt protein variation
							unsolved += 1
							#print(protein_predictions, line.split('\t')[2], file=sys.stdout)
							print('\t'.join([line, 'not_found']), file=vcf_output)
						c += 1
				else:
					unsolved += 1
					print('\t'.join([line, 'unproper']), file=vcf_output)
		print(str(variant_in_reference) + ' have been found as reference AA in the primary protein sequence, corresponding to ' + str(variant_in_reference/tot_counter * 100) + '% of the total variants')
		print(str(completely_solved) + ' have been completely solved as common polymorphism, corresponding to ' + str(completely_solved/tot_counter * 100) + '% of the total variants')
		print(str(probable_solved) + ' have been probably solved, corresponding to ' + str(probable_solved/tot_counter * 100) + '% of the total variants')
		print(str(possible_solved) + ' have been possibly solved, corresponding to ' + str(possible_solved/tot_counter * 100) + '% of the total variants')
		print(str(unsolved) + ' have not been solved, corresponding to ' + str(unsolved/tot_counter * 100) + '% of the total variants')
	#embed()
	

def analyze_vep_vcf(vcf_input_file, vcf_output, uniprot_dict, dbsnp_dict, dbsnp_codes, uniprot_seq):
	'''Function working on SnpEff annotated file to establish if protein variation are present within UniProt curated DB'''
	uniprot_seq_index = uniprot_seq.index
	with open(vcf_input_file, 'r') as vcf_file:
		tot_counter = 0
		unsolved = 0
		completely_solved = 0
		probable_solved = 0
		possible_solved = 0
		variant_in_reference = 0 
		for line in vcf_file:
			line = line.rstrip('\n')
			#print(line)
			if line.startswith('#'):
				print(line, file=vcf_output)
				continue
			else:
				columns = line.split('\t')
				rs_code = columns[2]
				tot_counter += 1 
				protein_predictions, protein_length = get_protein_predictions(line, 'vep')
				if protein_predictions:
					c = 1
					for i, prediction in enumerate(protein_predictions):
						gene_name, hgvs_protein = prediction.split('|')
						if hgvs_protein == '':	#variant with no protein prediction
							if c == len(protein_predictions):
								unsolved += 1
								print('\t'.join([line, 'not_found']), file=vcf_output)
							c+=1
							continue
						first_aa, pos, second_aa = re.findall(AAandPOSpattern, hgvs_protein)[0]
						# prediction --> GENE_NAME | p.AA1posAA2
						if prediction in uniprot_dict:
							print('\t'.join([line, uniprot_dict[prediction]]), file=vcf_output)
							completely_solved += 1
							break
						elif get_dbsnp_prediction(rs_code, prediction) in dbsnp_dict:
							print('\t'.join([line, dbsnp_dict[get_dbsnp_prediction(rs_code, prediction)]]), file=vcf_output)
							completely_solved += 1
							break
						elif inverted_prediction(prediction) in uniprot_dict:
							print('\t'.join([line, 'variant_is_in_reference']), file=vcf_output)
							variant_in_reference += 1
							break
						elif rs_code in dbsnp_codes and checkAAidentity(rs_code, dbsnp_codes, prediction):
							print('\t'.join([line, 'probable_' + dbsnp_codes[rs_code][-1]]),  file=vcf_output)
							probable_solved += 1
							break
						elif rs_code in dbsnp_codes and checkAApos(dbsnp_codes[rs_code][0], prediction): 
							print('\t'.join([line, 'possible_' + dbsnp_codes[rs_code][-1]]),  file=vcf_output)
							possible_solved += 1
							break
						elif gene_name in uniprot_seq_index	and type(uniprot_seq.loc[gene_name,'sequence'])!=str and args.isoforms:
							#if uniprot_seq.loc[gene_name].sequence.apply(lambda x:three_letters_code[x[int(pos)-1]]==second_aa).any():
							#if uniprot_seq.loc[gene_name].sequence.apply(checkAAintoSeq, args=(pos, second_aa)).any():
							if uniprot_seq.loc[gene_name].sequence.apply(checkAAintoSeq, args=(pos, second_aa, protein_length[i])).any() and first_aa!=second_aa:
									#print(prediction, uniprot_seq.loc[gene_name, 'sequence'].apply(len),protein_length[i])
									print('\t'.join([line, 'variant_is_in_reference']), file=vcf_output)
									variant_in_reference += 1
									#completely_solved += 1
									#probable_solved += 1
									break
							elif c == len(protein_predictions):
								unsolved += 1
								print('\t'.join([line, 'not_found']), file=vcf_output)
						elif gene_name in uniprot_seq_index	and type(uniprot_seq.loc[gene_name,'sequence'])==str and args.isoforms:
							#if checkAAintoSeq(uniprot_seq.loc[gene_name].sequence, pos, second_aa):
							if checkAAintoSeq(uniprot_seq.loc[gene_name].sequence, pos, second_aa,  protein_length[i]) and first_aa!=second_aa:
									#print(uniprot_seq.loc[gene_name, 'sequence'], prediction)
									print('\t'.join([line, 'variant_is_in_reference']), file=vcf_output)
									#completely_solved += 1
									variant_in_reference += 1
									#probable_solved += 1
									break
							elif c == len(protein_predictions):
								unsolved += 1
								print('\t'.join([line, 'not_found']), file=vcf_output)
								

						elif c == len(protein_predictions):	# if the script comes inside here it means that SnpEff prediction is not inside common UniProt protein variation
							unsolved += 1
							#print(protein_predictions, line.split('\t')[2], file=sys.stdout)
							print('\t'.join([line, 'not_found']), file=vcf_output)
						c += 1
				else:
					unsolved += 1
					print('\t'.join([line, 'unproper']), file=vcf_output)
		print(str(variant_in_reference) + ' have been found as reference AA in the primary protein sequence, corresponding to ' + str(variant_in_reference/tot_counter * 100) + '% of the total variants')
		print(str(completely_solved) + ' have been completely solved as common polymorphism, corresponding to ' + str(completely_solved/tot_counter * 100) + '% of the total variants')
		#print(str(completely_solved) + ' have been completely solved, corresponding to ' + str(completely_solved/tot_counter * 100) + '% of the total variants')
		print(str(probable_solved) + ' have been probably solved, corresponding to ' + str(probable_solved/tot_counter * 100) + '% of the total variants')
		print(str(possible_solved) + ' have been possibly solved, corresponding to ' + str(possible_solved/tot_counter * 100) + '% of the total variants')
		print(str(unsolved) + ' have not been solved, corresponding to ' + str(unsolved/tot_counter * 100) + '% of the total variants')

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="A tool to check if missense variant in an annotated VCF are included in UniProt common protein variantions.")
	parser.add_argument("-U", "--uniprot", metavar='UNIPROT', required=True, help="Uniprot file containing common protein variation")
	parser.add_argument("-V", "--variants", metavar='VARS', required=True, help="Annotated VCF file")
	parser.add_argument("-I", "--isoforms", metavar='ISO', required=True, help="Uniprot converted file containing protein isoform sequence")
	#parser.add_argument("-s", "--sequence-length", metavar='SL', help="Uniprot converted file containing protein length and sequence")
	parser.add_argument("-A", "--annotator", metavar="ANN", required=True, choices=["snpeff", "vep"], help="Parameter defining which variant annotator has been used") 
	args = parser.parse_args()
#	if args.isoforms and args.sequence_length:
#		sys.exit('Options --isoforms [-i] and --sequence-length [-s] are mutually exclusive. Please select only one!')
	uniprot_dict, dbsnp_dict, dbsnp_codes = parse_uniprot(args.uniprot)
#	if args.isoforms:
	uniprot_seq = load_uniprot_isoforms(args.isoforms)
#	elif args.sequence_length:
#		uniprot_seq = load_uniprot_sequence_length(args.sequence_length)
	output_file = args.variants.split('.vcf')[0] + '.UniProtChecked.vcf'
	vcf_output = open(output_file, 'w') 
	if args.annotator == 'snpeff':
		analyze_snpeff_vcf(args.variants, vcf_output, uniprot_dict, dbsnp_dict, dbsnp_codes, uniprot_seq)
	elif args.annotator == 'vep':
		analyze_vep_vcf(args.variants, vcf_output, uniprot_dict, dbsnp_dict, dbsnp_codes, uniprot_seq)
		#sys.exit('not yet implemented')
	
