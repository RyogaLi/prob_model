#!/usr/bin/python
# coding=utf-8
# v2.0
from read_files import *
from supplementary import *
from prob_model import *
from conf import *
import os
import fnmatch
import subprocess

import argparse
import logging.config

def load_pickle(filename):
	with open(filename, 'rb') as pkl_file:
		filecontent = pickle.load(pkl_file)
	return filecontent

def load_npy(filename):
	filecontent = np.load(filename)
	return filecontent

def save_as_matrix(filename):
	"""
	save file as matrix
	:param filename: 
	:return: 
	"""
	matrix = []
	with open(filename, "r") as mix:
		for line in mix:
			line = line.strip().split(",")
			line = [i.strip("\"\"\n") for i in line]
			matrix.append(line)
	return np.asarray(matrix)

def find(pattern, path):
	"""
	Find file in path that with filename matches pattern
	:param pattern:
	:param path:
	:return: file path
	"""
	result = None
	for root, dirs, files in os.walk(path):
		for name in files:
			if fnmatch.fnmatch(name, pattern):
				result = os.path.join(root, name)
	return result

def write_output_to_file(filename, data):
	# write content in data by line to filename
	with open(filename, "w") as output:
		for item in data:
			output.write("%s\n" % item)

def read_tumour_spreadsheet(input_file):
	"""
	read the spread sheet and save it into numpy arrray
	:param input_file: csv file contains tumour names and chromatin file names
	:return: matrix
	"""
	spreadsheet = []
	with open(input_file, "r") as spreads:
		for line in spreads:
			if "sample_name" not in line:
				spreadsheet.append(line.split())
	return spreadsheet


def get_line(title, matrix):
	"""
	get a line by title(first element)
	:param title: 
	:return: list contains the title and rest of the line 
	"""
	return matrix[np.where(matrix[:, 0] == title)]

def three_fold_validation():
	# three fold validation
	for variants_file in vcf_list[group]:
		# select only breast cancer
		bc = []
		with open("/home/ryogali/dev/prob_model/BRCA_files.txt") as files:
			for line in files:
				line = line.strip("\"\"\n")
				bc.append(line)

		if variants_file.endswith(".vcf"):
			tumour_name = variants_file.split(".")[0]
		else:
			main_logger.debug("Only VCF files are accepted.")
			continue

		if tumour_name not in bc: # not a brca tumor
			print("not bc")
			continue
		# # find mixture directory
		# mixture = None
		# mixture_path = os.path.join(mixture_dir, tumour_name)
		# # check if tumour exist in mixture path
		# if os.path.exists(mixture_path):
		# 	mixture = find("mixtures.csv", mixture_path)
		# # check if mixture.csv found
		# if mixture == None:
		# 	continue
		#
		# main_logger.info("Tumor name: %s", tumour_name)
		# # get signatures and time points
		# sigs, mixture = read_mix_file(mixture)
		# main_logger.info("%s signatures analyzed", str(len(sigs)))

		# to do get new mixture here
		sigs = []
		mixture = []
		# convert file to matrix
		mixture_matrix = save_as_matrix(mixture_overall)
		# print(mixture_matrix)
		# print(tumour_name)
		# select tumor name from the matrix
		tumor_sig = mixture_matrix[(mixture_matrix[:,0]==tumour_name)|(mixture_matrix[:,0]=="")]
		# print(tumor_sig)
		if tumor_sig.shape[0]<2:
			continue
		# select where signatures != 0
		for i in range(len(tumor_sig[1])):
			# print(i)
			if tumor_sig[1][i] != "0":
				# print(tumor_sig[1][i])
				sigs.append(tumor_sig[0][i])
				mixture.append(tumor_sig[1][i])
		for i in range(len(sigs)):
			sigs[i] = "Signature " + sigs[i]

		mixture = [float(i) for i in mixture[1:]]
		# print(sigs)
		# print(mixture)

		##################################################################
		##################################################################

		vf = os.path.join(input_dir, variants_file)
		variants_parser = VariantsFileParser(vf, chromatin_file, mRNA_file, hg19_file, trinuc, mixture, alex_signature_file, sigs)

		# get input data to the model
		# n = n fold validation, 1/n as train data and 2/n as test data
		# low support data are those mutations that has FILTER = LOWSUPPORT
		test, train, low_support = variants_parser._get_input_data(3)

		# get low support data feature
		low_support_data = variants_parser._get_features(low_support)

		# get random data feature
		if chromatin_file:
			random_data = generate_data(mixture, alex_signature_file, sigs, True)
		else:
			random_data = generate_data(mixture, alex_signature_file, sigs, False)

		for i in range(len(train)):
			# get features from train data
			train_data = variants_parser._get_features(train[i])
			# get features from test data
			test_data = variants_parser._get_features(test[i])

			# train the model
			train_matrix = ProbModel(train_data)
			train_matrix._fit()

			# predict probabilities for train data
			if chromatin_file:
				train_pa, train_pt, train_ps = \
                                        train_matrix._predict_proba(train_matrix._mut,
                                                                    train_matrix._tr_X,
                                                                    train_matrix._strand_X,
                                                                    train_matrix._strand)
				test_matrix = ProbModel(test_data)
				test_pa, test_pt, test_ps = \
                                        train_matrix._predict_proba(test_matrix._mut,
                                                                    test_matrix._tr_X,
                                                                    test_matrix._strand_X,
                                                                    test_matrix._strand)
				# predict probabilities for low_sup data
				low_support_matrix = ProbModel(low_support_data)
				lowsup_pa, lowsup_pt, lowsup_ps = train_matrix._predict_proba(low_support_matrix._mut,
                                                                  low_support_matrix._tr_X,
                                                                  low_support_matrix._strand_X,
                                                                  low_support_matrix._strand)
				# predict probabilities for random data
				random_matrix = ProbModel(random_data)
				random_pa, random_pt, random_ps = train_matrix._predict_proba(random_matrix._mut,
                                                                        random_matrix._tr_X,
                                                                        random_matrix._strand_X,
                                                                        random_matrix._strand)
			else:
				train_pa, train_pt, train_ps = train_matrix._predict_proba(train_matrix._mut,
                                                            train_matrix._mut,
                                                               train_matrix._strand_X,train_matrix._strand)

				# predict probabilities for test data
				test_matrix = ProbModel(test_data)
				test_pa, test_pt, test_ps = train_matrix._predict_proba(test_matrix._mut,
                                                            test_matrix._mut ,
                                                            test_matrix._strand_X,
                                                           test_matrix._strand)
				# predict probabilities for low_sup data
				low_support_matrix = ProbModel(low_support_data)
				lowsup_pa, lowsup_pt, lowsup_ps = train_matrix._predict_proba(low_support_matrix._mut,
                                                            low_support_matrix._mut,
                                                                  low_support_matrix._strand_X,
                                                                 low_support_matrix._strand)
				# predict probabilities for random data
				random_matrix = ProbModel(random_data)
				random_pa, random_pt, random_ps = train_matrix._predict_proba(random_matrix._mut,
                                                            random_matrix._mut,
                                                                  random_matrix._strand_X,
                                                                 random_matrix._strand)

			# write the probabilities to file
			write_output_to_file(os.path.join(train_prob_dir,
                                     tumour_name)+ "."+str(i)+".train.txt",
                                     train_matrix._calculate_proba(train_pa, train_pt, train_ps))
			write_output_to_file(os.path.join(test_prob_dir,
                                     tumour_name)+"."+str(i)+".test.txt",
                                     test_matrix._calculate_proba(test_pa, test_pt, test_ps))
			write_output_to_file(os.path.join(lowsup_prob_dir,
                                     tumour_name)+"."+str(i)+".lowsup.txt",
                                     low_support_matrix._calculate_proba(lowsup_pa, lowsup_pt, lowsup_ps))
			write_output_to_file(os.path.join(random_prob_dir,
                                     tumour_name)+"."+str(i)+".random.txt",
                                     random_matrix._calculate_proba(random_pa, random_pt, random_ps))
			main_logger.info("DONE-%s-%s",tumour_name, str(i))

def single_file_main():

	# get mixture matrix (one file for all the tumours)
	mixture_matrix = save_as_matrix(mixture_overall)
	vf = os.path.join(vcf_file_path, vcf_file)



def validated_file_main():
	# get tumour name
	tumour_name = vcf_file.split(".")[0]
	# get mixture matrix (one file for all the tumours)
	mixture_matrix = save_as_matrix(mixture_overall)
	vf = os.path.join(vcf_file_path, vcf_file)
	# parse input file
	validated_parser = ValidatedVCFParser(vf)
	validated_variants = validated_parser._parse()

	# find overlaps between these variants and training data
	# remove overlaps from validated matrix


if __name__ == '__main__':
	# get logger
	logging.config.fileConfig("src/logging.conf")
	main_logger = logging.getLogger("main")

	# all mixture in one file
	mixture_overall = "./data/overall_exposures.sigsBeta2.csv"

	# logging configration
	parser = argparse.ArgumentParser(description='Predict mutation probabilities')
	parser.add_argument('-f', '--file', help='To run the model on single vcf file', required=False)
	parser.add_argument("-v", "--validated", help="To run the model on validated (deep sequenced) tumours", required=False)


	##################################################################
	# this argument is used for spreading the program on different nodes
	parser.add_argument('--group', default=-1,type=int, required=False)
	args = parser.parse_args()

	##################################################################

	# create directory for different class of variants
	# train_prob_dir = os.path.join(output_dir, "train_prob/")
	# test_prob_dir = os.path.join(output_dir, "test_prob/")
	# random_prob_dir = os.path.join(output_dir, "random_prob/")
	# lowsup_prob_dir = os.path.join(output_dir, "lowsup_prob/")

	main_logger.info("output files will be saved into: %s", output_dir)

	##################################################################
	# Load feature files
	feature_data = "./data/"
	main_logger.info("Loading required feature files...")
	try:
		mRNA_file = load_pickle(os.path.join(feature_data, "mRNA.pickle"))
		#main_logger.info("mRNA loaded")
		trinuc = load_pickle(os.path.join(feature_data,"trinucleotide.pickle"))
		#main_logger.info("trinuc loaded")
		alex_signature_file = load_npy(os.path.join(feature_data,"signature.npy"))
		#main_logger.info("alex_signature loaded")
		hg19_file = load_pickle(os.path.join(feature_data,"hg.pickle"))
		#main_logger.info("hg file loaded")
	except Exception as error:
		main_logger.exception("Please provide valid compiled feature data files.")
		exit()
	##################################################################

	all_vcf = os.listdir(vcf_file_path)

	# following code is for paralizing jobs
	if args.group is not None:
		group = args.group
		GROUP_SIZE = 147 #change this based on cores
		vcf_list = [all_vcf[i:i + GROUP_SIZE] for i in xrange(0, len(all_vcf), GROUP_SIZE)]
	else:
		group = 0
		vcf_list = [all_vcf]
	# end

	##################################################################
	# read spread sheet and run the model
	spreadsheet = read_tumour_spreadsheet(tumour_type)
	# for each tumour in the spreadsheet
	# get the tumour name and chromatin profile for this tumour
	for line in spreadsheet:
		tumour_id = line[0]
		vcf_file = subprocess.check_output('find '+vcf_file_path+' -name '+tumour_id+"*",shell=True)
		main_logger.info("Processing tumour: %s", vcf_file)
		main_logger.info("Tumour type: %s", line[1])
		if line[2] == "N/A":
			chromatin_file = 0
			main_logger.info("No chromatin profile specified for this tumour")
		else:
			chromatin_file = os.path.join(chromatin_path, line[2]+".bed")
			main_logger.info("chromatin profile: %s", chromatin_file)
			chromatin_dict = read_chromatin(chromatin_file)
		# run model on this vcf file
		if args.validated == "True":
			validated_file_main()
		else:
			single_file_main()
	##################################################################




