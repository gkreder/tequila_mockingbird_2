import pickle
from Bio.Seq import Seq


def barcodes_to_alleles(filename):

	# barcode_count_dict = {}

	# for barcode in test_data:
	# 	# print line
	# 	# barcode = line.rstrip()
	# 	if barcode in barcode_count_dict:
	# 		barcode_count_dict[barcode] += 1
	# 	else:
	# 		barcode_count_dict[barcode] = 1

	# print barcode_count_dict

	# -----------------------------------------------------------------
	# with file input
	# -----------------------------------------------------------------
	barcode_count_dict = {}

	with open(filename, 'rb') as f:
		lines = f.readlines()
	f.close()
	total_count = 0
	num_of_N = 0
	mapped_count = 0

	for line in lines:
		# print line
		barcode = line.rstrip()
		barcode = barcode[0 : 18]
		barcode = Seq(barcode)
		barcode = barcode.reverse_complement()
		barcode = str(barcode)
		total_count += 1
		if barcode in barcode_count_dict:
			barcode_count_dict[barcode] += 1
		else:
			barcode_count_dict[barcode] = 1

		


	# -----------------------------------------------------------------




	# -----------------------------------------------------------------

	# # barcode_allele_dict = ???

	allele_count_dict = {}
	# barcode_allele_dict = pickle.load(open('allele_dic_with_WT.pkl', 'rb'))
	barcode_allele_dict = pickle.load(open('allele_dic.pkl', 'rb'))
	# print barcode_count_dict
	# print barcode_allele_dict
	for barcode in barcode_count_dict:
		barcode = line[0 : 18]
		barcode = Seq(barcode)
		barcode = barcode.reverse_complement()
		barcode = str(barcode)
		if 'N' in barcode:
			num_of_N += 1
			continue
		# allele_count_dict[barcode_allele_dict(barcode)] = barcode_count_dict[barcode]
		elif barcode in barcode_allele_dict.keys():
			mapped_count += 1
			allele_count_dict[barcode_allele_dict[barcode]] = barcode_count_dict[barcode]


	return allele_count_dict,total_count,num_of_N,mapped_count






# ----------------------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------------------


filename_list = ['t0_r1_seqs.txt',
				 't1_r1_seqs.txt',
				 't2_r1_seqs.txt',
				 't0_r2_seqs.txt',
				 't1_r2_seqs.txt',
				 't2_r2_seqs.txt',
				 't0_control_seqs.txt',
				 't1_control_seqs.txt',
				 't2_control_seqs.txt']

for filename in filename_list:
	allele_count_dict,total_count,num_of_N,mapped_count = barcodes_to_alleles(filename)
	# print allele_count_dict
	filename_save = filename.replace('_seqs.txt', '_allele_count_3.pkl')
	# print filename, filename_save

	pickle.dump(allele_count_dict, open(filename_save, "wb" ))

	# print allele_count_dict

	# print total counts, number_of_N, mapped_count
	print filename, total_count, num_of_N, mapped_count

	# print allele_count_dict