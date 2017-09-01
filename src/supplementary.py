#!/usr/bin/python
# coding=utf-8
#v2.0
import random
from parse_input_data import *
random.seed(1)

def generate_data(mixture, alex_signature, sigs, chromatin_bool):
	"""
	generate random data corresponding to each variants file
	for each variants file create a random matrix
	chr | pos | vaf | reference | alt | tri | mut_type
	return: mut_type, se, trans, sense, chromatin, p_ce, VAF_list
	"""
	variant_parser = VariantParser()
	# total random mutations
	total = 5000

	# trans = 0/1
	trans = np.asarray([[random.randrange(0,2,1)] for _ in range (total)])

	# chromatin = 0/1
	if chromatin_bool:
		chromatin = np.asarray([[random.randrange(0,2,1)] for _ in range (total)])
	else:
		chromatin = np.asarray([-1]*total).reshape(total, 1)

	sense = []
	for i in trans:
		if i == 0:
			sense.append([-1])
		else:
			sense.append([random.randrange(0,2,1)])
	sense = np.asarray(sense)

	VAF_list = [random.uniform(0, 1)*2 for _ in range (total)]

	se = []
	p_ce = []
	mut_type = []

	for i in VAF_list:
		mut_class = random.randrange(0,96,1)
		## get signature vector from time point file
		# vaf_values = mixture[0, :]
		# # print VAF
		# idx = (np.abs(vaf_values - i)).argmin()
		# # print idx
		# signature_vector = mixture[:, idx][1:]
		signature_vector = mixture
		se.append(signature_vector)
		pce = variant_parser._calculate_pce(signature_vector,
                                      alex_signature, sigs, mut_class)
		p_ce.append(pce)
		mut_type.append(mut_class)

	mut_type_one_hot = get_one_hot_encoding(mut_type)

	mut_type = add_colname(np.asarray(mut_type_one_hot), "Mut_type")

	se = add_colname(np.asarray(se), "Exposure")
	trans_region = add_colname(trans, "Transcribed")
	sense = add_colname(sense, "Strand")

	VAF_list = np.asarray(VAF_list).reshape(np.asarray(VAF_list).shape[0], 1)
	VAF_list = add_colname(np.asarray(VAF_list), "VAF")

	chromatin = add_colname(np.asarray(chromatin), "Chromatin")

	p_ce = np.asarray(p_ce).reshape(np.asarray(p_ce).shape[0], 1)
	p_ce = add_colname(np.asarray(p_ce), "p_ce")

	return combine_column([mut_type, se, trans_region, sense, chromatin,p_ce, VAF_list])

# if __name__ == "__main__":
# 	generate_data()

