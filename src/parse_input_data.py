#!/usr/bin/python
# coding=utf-8
#v2.1
from __future__ import division
from sklearn.model_selection import KFold
from prob_model import *
import vcf
import math
np.random.seed(1)

def get_one_hot_encoding(matrix,max):
	"""
	get an one hot encoding matrix with corresponding position = 1
	matrix is zero indexed
	"""
	output = []
	for i in matrix:
		cla = int(i)
		one_hot = [0] * max
		one_hot[cla] = 1
		output.append(one_hot)
	output = np.asarray(output)
	return output

def add_colname(matrix, colname):
	"""
	add colname to matrix
	:param matrix:
	:param colname:
	:return:
	"""
	colnames = np.repeat(colname, matrix.shape[1])
	colnames = colnames.reshape(1, matrix.shape[1])
	matrix = np.concatenate((colnames, matrix))
	return matrix

class VariantsFileParser(object):

	def __init__(self, filename, chromatin_dict, mrna_dict, hg, trinuc_file, time_points, alex_sig, signatures):
		self._filename = filename
		self._chromatin_dict = chromatin_dict
		self._mrna_dict = mrna_dict
		self._hg = hg
		self._trinuc_file = trinuc_file
		self._time_points = time_points # mixture file
		self._alex_sig = alex_sig
		self._signatures = signatures

	def _save_as_matrix(self):
		"""
		create a matrix for the vcf data
		:return:
		"""
		matrix = []
		low_sup_matrix = []
		vcf_reader = vcf.Reader(open(self._filename))
		for record in vcf_reader:
			if "VAF" in record.INFO.keys():
				record.INFO["VAF"] = float(record.INFO["VAF"]) * 2
			else:
				continue

			if record.FILTER == ["LOWSUPPORT"]:
				low_sup_matrix.append(record)
			else:
				matrix.append(record)
		matrix = np.asarray(matrix)
		low_sup_matrix = np.asarray(low_sup_matrix)
		return matrix, low_sup_matrix

	def _nfold_shuffle(self, matrix, n):
		"""
		separate data based on n
		"""
		train =[]
		test = []
		kf = KFold(n_splits=n, shuffle=True, random_state=10)
		for train_index, test_index in kf.split(matrix):
			train.append(matrix[train_index])
			test.append(matrix[test_index])
		return train, test

	def _get_input_data(self, n=0, validate=False):
		"""
		read variants file and separate train, test, lowsupport data
		:param n: n-fold cross validation
		:return: train, test, lowsupport
		"""
		data, low_support = self._save_as_matrix()
		# print data.shape
		if validate != False:
			# germline = []
			notseen = []
			passed = []
			train = []
			for record in data:
				chromosome = "chr" + str(record.CHROM)
				position = record.POS-1
				validated = validate[np.where((validate[:, 0] == chromosome) * (validate[:, 1] == str(position)))]
				if validated.shape[0] !=0:
					if validated[0][2] == "PASS":
						passed.append(record)
					else:
						notseen.append(record)
				else:
					train.append(record)
			return np.asarray(train), np.asarray(passed), np.asarray(notseen)
		else:
			# you can choose n fold validation
			if n != 0:
				train, test = self._nfold_shuffle(data, n)
			# if n == 0, train and test data = 2:1 (DEFAULT)
			else:
				np.random.shuffle(data)
				total = data.shape[0]
				sep = math.floor(total*0.67)
				train, test = data[:sep,], data[sep:, ]
			return train, test, low_support

	def _get_features(self, input_data):
		"""
		Get all the features from input data
		:param input_data: 
		:return:
		"""
		chromatin = []
		trans_region = []
		sense = []
		mut_type = []
		VAF_list = []
		se = []
		p_ce = []
		for record in input_data:
			variant_parser = VariantParser(record)
			chromosome = "chr" + str(record.CHROM)
			position = record.POS - 1
			if chromosome != "chrY":
				if self._chromatin_dict:
					find_region = variant_parser._find_variant_in_region(self._chromatin_dict[chromosome])
					if find_region:
						chromatin.append([1])
					else:
						chromatin.append([0])
				else:
					chromatin.append([-1])

				trans_interval = variant_parser._find_variant_in_region(self._mrna_dict[chromosome])
				# trans.append(variant_parser._find_variant_in_region(mrna_dict[chromosome]))
				if trans_interval:
					trans_region.append([1])
					if trans_interval[2] == "+":
						if record.REF == "G" or record.REF == "A":  # the mutation is on antisense strand
							sense.append([1])
						else:  # the mutation is on sense strand
							sense.append([0])
					else:
						if record.REF == "T" or record.REF == "C":
							sense.append([1])
						else:
							sense.append([0])
				else:
					trans_region.append([0])
					sense.append([-1])

				trinucleotide_class = variant_parser._get_trinucleotide(self._hg)
				int_type = self._trinuc_file[trinucleotide_class]

				mut_type.append(int_type)

				VAF = float(record.INFO["VAF"])
				VAF_list.append(VAF)


				## get signature vector from time point file
				# vaf_values = self._time_points[0, :]
				# idx = (np.abs(vaf_values - VAF)).argmin()
				# print idx
				# signature_vector = self._time_points[:, idx][1:]
				# print signature_vector

				## get signature from overall
				signature_vector = self._time_points
				# print(signature_vector)
				se.append(signature_vector)

				pce = variant_parser._calculate_pce(signature_vector, self._alex_sig, self._signatures, int_type)

				p_ce.append(pce)

		mut_type_one_hot = get_one_hot_encoding(mut_type, 96)
		mut_type = add_colname(np.asarray(mut_type_one_hot), "Mut_type")
		se = add_colname(np.asarray(se), "Exposure")
		trans_region = add_colname(np.asarray(trans_region), "Transcribed")
		sense = add_colname(np.asarray(sense), "Strand")
		VAF_list = np.asarray(VAF_list).reshape(np.asarray(VAF_list).shape[0], 1)
		VAF_list = add_colname(np.asarray(VAF_list), "VAF")

		chromatin = add_colname(np.asarray(chromatin), "Chromatin")

		p_ce = np.asarray(p_ce).reshape(np.asarray(p_ce).shape[0],1)
		p_ce = add_colname(np.asarray(p_ce), "p_ce")


		return combine_column([mut_type, se, trans_region, sense, chromatin, p_ce, VAF_list])



class ValidatedVCFParser(object):

	def __init__(self, filename):
		self._filename = filename

	def _parse(self):
		"""
		save validated data into matrix
		:return: a matrix with chr | pos | filter
		"""
		validated_matrix = []
		vcf_reader = vcf.Reader(open(self._filename))
		for record in vcf_reader:
			chromosome = record.CHROM
			position = record.POS - 1
			filter = record.FILTER
			if filter == []:
				filter = "PASS"
			else:
				filter = filter[0]
			validated_matrix.append([chromosome, position, filter])
		# print validated_matrix
		validated_matrix = np.asarray(validated_matrix)
		return validated_matrix


class VariantParser(object):

	def __init__(self, variant=None):
		self._variant = variant

	def _get_trinucleotide(self, hg19):
		"""
		get trinucleotide content from hg19 file
		:param hg19: load from hg19.pickle
		:return:
		"""
		reference = str(self._variant.REF)
		alt = str(self._variant.ALT[0])
		position = self._variant.POS - 1
		mut_pair = {"G": "C", "T": "A", "A": "T", "C": "G"}
		sequence = hg19[self._variant.CHROM] # sequence of corresponding chromosome
		if reference == "G" or reference == "A":
			reference = mut_pair[reference]
			alt = mut_pair[alt]
			tri = mut_pair[sequence[position - 1].upper()] + \
				  mut_pair[sequence[position].upper()] + mut_pair[sequence[position + 1].upper()]
		else:
			tri = sequence[position - 1].upper() + sequence[position].upper() \
				  + sequence[position + 1].upper()

		return (reference, alt, tri)

	def _find_variant_in_region(self, input_list):
		"""
		Use binary serach to find the range
		:param input_list: list of tuples : (start, end)
		:return:
		"""
		position = self._variant.POS - 1
		input_list.sort()
		start = 0
		end = len(input_list) - 1
		while end >= start:
			mid = (end + start) // 2
			if position in range(input_list[mid][0], input_list[mid][1] + 1):
				return input_list[mid]
			elif position > input_list[mid][1]:
				start = mid + 1
			elif position < input_list[mid][0]:
				end = mid - 1
		return None

	def _calculate_pce(self, exposure_vector, alex_signature, sigs, mut_type):
		"""
		calculate p_ce for a given mutation

		:param exposure_vector: exposure at each time point
		:param alex_signature: alex signature file, rows corresponding to mut_type and cols=signatures
		:param sigs: active signature in a given tumor
		:param mut_type: 96 types of mutations (index: 0-95)
		:return:
		"""

		idx = np.where(np.in1d(alex_signature, sigs))
		select_signature_cols = alex_signature[:,2:][mut_type+1,idx].squeeze()
		exposure_vector = [float(i) for i in exposure_vector]
		sum = np.dot(np.asarray(select_signature_cols,dtype=np.float), exposure_vector)
		return sum

