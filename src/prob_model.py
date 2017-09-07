#!/usr/bin/python
# coding=utf-8
#v2.0

import numpy as np
from sklearn.linear_model import LogisticRegression
import logging

logger = logging.getLogger("prob_model")

def lr(X, y):
	logistic_regression = LogisticRegression(penalty="l2")
	logistic_regression.fit(X, y)
	return logistic_regression

def convert_onehot(one_hot_matrix):
	"""
	convert one hot encoding matrix back to vector
	:param one_hot_matrix:
	:return:
	"""
	vector = []
	for i in one_hot_matrix:
		vector.append(int(np.where(i == 1)[0]))
	return vector

def combine_column(l):
	"""
	combine all the matrix in l
	Ex.
	> a = np.array([[1,1],[2,2],[3,3],[4,4],[5,5]])
	> b = np.array([["a"],["a"],["a"],["a"],["a"]])
	> l = [a,b]
	> add_column(l)
	> [['1' '1' 'a']
	> ['2' '2' 'a']
	> ['3' '3' 'a']
	> ['4' '4' 'a']
	> ['5' '5' 'a']]
	"""
	out = l[0]
	l = l[1:]
	for i in l:
		out = np.append(out, i, 1)
	return out

def write_output_to_file(filename, data):
	with open(filename, "w") as output:
		for item in data:
			output.write("%s\n" % item)

def get_col(colname, matrix):
	"""
	return columns with colname in matrix
	"""
	cols = matrix[:, np.where(matrix[0,] == colname)]
	return cols[1:,:].astype(float).squeeze()

class ProbModel(object):

	def __init__(self, input_matrix):
		self._input_matrix = input_matrix # whole matrix
		self._mut_type = get_col("Mut_type", input_matrix) # mut type (one hot encoded)
		self._exposure = get_col("Exposure", input_matrix) # exposure vector
		self._mut = combine_column([self._mut_type, self._exposure]) # combine mut_type and exposure vector for later use
		self._transregion = get_col("Transcribed", input_matrix) # transcribed
		self._strand =get_col("Strand", input_matrix) #strand information
		self._chromatin = get_col("Chromatin", input_matrix) #chromatin accessibility
		self._p_ce = get_col("p_ce", input_matrix) # calculated p_ce vector
		self._vaf = get_col("VAF", input_matrix) # vaf: VAF to actural VAF * 2
		# combine mut and chromatin for later use
		self._tr_X = combine_column([self._mut, self._chromatin.reshape(self._chromatin.shape[0], 1)])
		# combine mut and strand for later use
		self._strand_X = combine_column([self._mut, self._transregion.reshape(self._transregion.shape[0], 1)])

	def _fit_chromatin(self):
		self._ch_lr = lr(self._mut, self._chromatin)

	def _fit_transregion(self):
		if -1 in self._chromatin:
			self._trans_lr = lr(self._mut, self._transregion)
		else:
			self._trans_lr = lr(self._tr_X, self._transregion)

	def _fit_strand(self):
		# select rows not equal to -1
		# get index of row != -1 for features and labels
		idx = np.where(self._strand != -1)[0]
		strand_X = self._strand_X[idx]
		labels = self._strand[idx]
		# fit the model
		self._strand_lr = lr(strand_X, labels)

	def _fit(self):
		if -1 not in self._chromatin:
			self._fit_chromatin()
		self._fit_strand()
		self._fit_transregion()

	def _predict_proba(self, mut, tr_X, strand_X, strand_label):
		if -1 in self._chromatin:
			p_a = np.asarray([1]*self._mut.shape[0]).reshape(self._mut.shape[0],1)
		else:
			p_a = self._ch_lr.predict_proba(mut)

		p_t = self._trans_lr.predict_proba(tr_X)
		idx = np.where(strand_label != -1)[0]
		p_s_predicted = self._strand_lr.predict_proba(strand_X[idx])
		#create empty array
		p_s = np.ones(shape=(strand_label.shape[0], 2), dtype=float)
		# insert predicted prob into p_s
		p_s[idx]=p_s_predicted.astype(float)

		return p_a, p_t, p_s

	def _calculate_proba(self,p_a, p_t, p_s):
		mut_prob = []
		mutation_type = convert_onehot(self._mut_type)
		mut_prob.append("mut\tVAF\tprob\tp_ai\tp_si\tp_ti\tp_ce")
		for i in range(len(self._p_ce)):
			p_ti = p_t[i][int(self._transregion[i])]

			if self._strand[i] == [-1]:
				p_si = 1
				#p_si = p_s[i][0]
			elif self._strand[i] == [0]:
				#print(p_s.shape)
				p_si = p_s[i][0]
			else:
				p_si = p_s[i][1]

			if -1 not in self._chromatin:
				p_ai = p_a[i][int(self._chromatin[i])]
			else: p_ai = 1
			mut_prob.append(str(mutation_type[i]) + "\t" + str(self._vaf[i]) + "\t" +
							str(p_ai * p_ti * p_si * self._p_ce[i]) + "\t" +
							str(p_ai) + "\t" + str(p_si) + "\t" +
							str(p_ti) + "\t" + str(self._p_ce[i]))
		return mut_prob

#
# if __name__ == "__main__":
# 	# test file
# 	#
	# self._filename = filename
	#
	# filecontent = np.load(self._filename)
	#
	# train = TrainFromMatrixFile(train_file)
	# train._fit()
	# p_a, p_t, p_s = train._predict_proba()
	# write_output_to_file("./train.txt", train._calculate_proba(p_a, p_t, p_s))
	#
	# test = TrainFromMatrixFile(test_file)
	# test_pa, test_pt, test_ps = train._predict_proba(test._mut, test._tr_X, test._strand_X)
	# write_output_to_file("./test.txt", test._calculate_proba(test_pa, test_pt, test_ps))
	#
	# low_support = TrainFromMatrixFile(lowsup_file)
	# low_support_pa, low_support_pt, low_support_ps = train._predict_proba(low_support._mut, low_support._tr_X, low_support._strand_X)
	# write_output_to_file("./low_sup.txt", low_support._calculate_proba(low_support_pa, low_support_pt, low_support_ps))
	# random = TrainFromMatrixFile(random_file)
	# random_pa, random_pt, random_ps = train._predict_proba(random._mut, random._tr_X, random._strand_X)
	# write_output_to_file("./random.txt", random._calculate_proba(random_pa, random_pt, random_ps))

