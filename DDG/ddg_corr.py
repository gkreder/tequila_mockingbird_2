import csv



ddg_affinities = []
with open('uby_OTU.tsv', 'rb') as csvfile:
	spamreader = csv.reader(csvfile, delimiter='\t')
	for row in spamreader:
		print row[3]
	if row[3] != 'None':
		ddg_affinities.append(float(row[3]))