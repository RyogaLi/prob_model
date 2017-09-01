#!/usr/bin/python
# coding=utf-8
from __future__ import division
import numpy as np
import argparse
import random as rd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import itertools
import math
import random as rd

# rd.seed(1)
# np.random.seed(1)


class SignatureFileParser(object):
	def __init__(self, filename):
		self._filename = filename

	def _parse(self):
		matrix = []
		with open(self._filename, "r") as input:
			signatures = input.readline().strip().split(",")
			matrix.append(signatures)
			for line in input:
				line = line.split(",")
				matrix.append(line)
		signatures = np.asarray(matrix)
		return signatures


class RandomData(object):
	def __init__(self, change_points, time_points, signatures):
		self._change_points = change_points
		self._time_points = time_points
		self._signatures = signatures

	def _generate_e_i(self, distribution="random"):
		"""
		generate e_i based on time points and signatures
		:return exposure
		>> t * n: t time points * n signatures matrix
		"""
		# generate constant value

		constant = np.random.dirichlet(np.ones(len(self._signatures)), size=1)
		exposure = constant.reshape(constant.shape[1], 1)
		if distribution == "constant":
			exposure = np.repeat(exposure, self._time_points, axis=1)
		else:
			# generate list of change points
			i = True
			while i:
				cp = rd.sample(range(3, self._time_points - 3), self._change_points)
				if 1 in [j - i for i, j in zip(constant[:-1], constant[1:])]:
					continue
				else:
					i = False

			exp = np.empty((len(self._signatures), 1))
			cp.insert(0, 0)
			cp.insert(-1, self._time_points)
			cp.sort()
			idx = 0
			while idx < len(cp)-1:
				constant = np.random.dirichlet(np.ones(len(self._signatures)), size=1)
				exposure = constant.reshape(constant.shape[1], 1)
				exp = np.concatenate((exp, np.repeat(exposure, cp[idx+1]-cp[idx], axis=1)), axis=1)
				idx += 1
			exposure = exp[:,1:]
		return exposure

	def _calculate_probability(self, signature_matrix, exposure):
		"""
		Calculate sum_i e_i S_i
		:param signature_matrix:
		:return: 96*t matrix containing mutation probabilities
		"""
		time_points = []

		selected_signatures = signature_matrix[1:, np.where(np.intersect1d(signature_matrix[0, :], self._signatures))]

		selected_signatures = selected_signatures.reshape(selected_signatures.shape[0], len(self._signatures)).astype(float)
		# print selected_signatures
		for i in exposure.T:
			time_points.append(np.sum(i * selected_signatures, axis=1))
		# print time_points
		# print np.asarray(time_points)
		return np.asarray(time_points)

	def _sample_mut_count(self, mutation_prob, mut_types):
		mut_count = []
		for col in mutation_prob:
			# print col
			mut_count.append(np.random.multinomial(100, col))
		return np.concatenate((mut_types[1:,:], np.asarray(mut_count, dtype=int).T), axis=1)

	def _sample_from_one_exposure(self, exposure,signature_matrix):

		# select signatures
		selected_signatures = signature_matrix[1:, np.where(np.intersect1d(signature_matrix[0, :], self._signatures))]
		selected_signatures = selected_signatures.reshape(selected_signatures.shape[0], len(self._signatures)).astype(
			float)
		mutations = []
		for i in exposure.T:
			mut_counts = []
			# for each exposure at this time point
			for index in range(len(i)):
				# get distribution for each signature
				distribution = i[index]*selected_signatures[:,index]
				# sample 100 * i[index] mutation for this signature
				print(int(round(100*i[index],0)))
				mut_vector = np.random.multinomial(int(round(100*i[index],0)), distribution)
				mut_counts.append(mut_vector)
				print(mut_vector)
			# sum the vecotr to get mut count at this time point (96*# of sig)
			mut_counts = np.asarray(mut_counts)
			mut_counts = np.sum(mut_counts,axis=0)
			# create a matrix for all the mut counts (96*t)
			mutations.append(mut_counts)
			break
		return np.asarray(mutations)

	def _plot(self, exposure, file_name):
		# plt.figure()
		for i in range(len(exposure)):
			plt.title("# of time points: "+ str(self._time_points))
			plt.ylabel("exposure")
			x = range(self._time_points)
			y = exposure[i]
			plt.plot(x, y, label=self._signatures[i], linestyle='--')
			plt.ylim(0,1)
			plt.xlim(0,self._time_points-1)

			plt.xticks(x, range(self._time_points))

		plt.legend(self._signatures, loc="upper right", ncol =2, fontsize="small")
		plt.grid()
		plt.savefig(file_name)
		plt.close()

def generate_signature_combinations(signature_list):
	"""

	:param signature_list:
	:return:
	"""
	# print signature_list
	combinations = list(itertools.combinations(signature_list, 2))
	return combinations
	# return signature_list

if __name__=="__main__":

	# input srguments
	parser = argparse.ArgumentParser(description="Generate random data")
	parser.add_argument("-dis", help="choose a distribution from ['constant', 'random'], default=random", default="random")
	parser.add_argument("-timepoint", help="number of time points", default=20)

	args = parser.parse_args()

	time_points = args.timepoint

	# signature file
	sig_file = "./PCAWG_signature_patterns_beta.csv"

	# save the file as matrix
	sig = SignatureFileParser(sig_file)
	signature_matrix = sig._parse()

	all_sigs = signature_matrix[0,2:] # signature names
	# select sig 1 and 5
	# idx = np.where(np.logical_or(all_sigs=='Signature 1', all_sigs=='Signature 5'))

	# select only sig 1
	# idx = np.where(all_sigs == 'Signature 1')

	# select only sig 5
	idx = np.where(all_sigs == 'Signature 5')

	# delete selected to generate all combindations
	all_sigs = np.delete(all_sigs, idx)
	combinations = generate_signature_combinations(all_sigs)

	# signatures = [['Signature 1', 'Signature 5'] + list(_) for _ in combinations]
	# signatures = [['Signature 1'] + list(_) for _ in combinations]
	signatures = [['Signature 5'] + list(_) for _ in combinations]

	for sig_comb in signatures:
		change_points = rd.randint(1, 3)
		# output file name
		# file_name = "./simulated_data/"+"_".join(sig_comb[2:])+ "_" + str(change_points)
		file_name = "./simulated_data/" + "_".join(sig_comb[1:]) + "_" + str(change_points)

		random = RandomData(change_points, time_points, sig_comb)
		exposure = random._generate_e_i(args.dis)

		sig_names = np.asarray(sig_comb).T.reshape(3,1)


		ex = exposure
		ex = np.append(sig_names, ex, axis=1)
		np.savetxt(file_name+".exposure.csv", ex,delimiter=",", fmt="%s")
		random._plot(exposure, file_name+".png")

		# mut_prob = random._calculate_probability(signature_matrix[:,2:], exposure)
		# mut_counts = random._sample_mut_count(mut_prob, signature_matrix[:,:2])
		#
		mut_counts = random._sample_from_one_exposure(exposure, signature_matrix[:,2:])
		# check sum and add/minus one
		sum = np.sum(mut_counts, axis=1)
		for i in range(len(sum)):
			if sum[i] == 99:
				# find a random col and add one
				col = rd.randint(0,95)
				mut_counts[i,:][col] = mut_counts[i,:][col] + 1
			elif sum[i] == 101:
				# find a count >= 1 and minus one
				mut_counts[i,:][np.where(mut_counts[i,:]>0)[0][0]] -= 1
		mut_counts = np.concatenate((signature_matrix[:,:2][1:,:], np.asarray(mut_counts, dtype=int).T), axis=1)
		np.savetxt(file_name+"mut_counts.csv",mut_counts, delimiter=",", fmt="%s")
		break