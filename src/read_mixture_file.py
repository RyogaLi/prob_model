#!/usr/bin/python
# coding=utf-8
#v2.0
import os
import fnmatch

def find(pattern, path):
	result = None
	for root, dirs, files in os.walk(path):
		for name in files:
			if fnmatch.fnmatch(name, pattern):
				result = os.path.join(root, name)
	return result

def read_mixture(mixture_dir, output_file):
	"""
	Go through every sub dir in the mixture dir
	:param mixture_dir:
	:return:
	"""
	os.chdir(mixture_dir)
	with open(output_file, "wa") as output:
		for directory in os.listdir(os.getcwd()):
			tumour_name = directory.split(".")[0]
			if os.path.isdir(tumour_name):
				mixture = find("mixtures.csv", directory)
				if mixture:
					with open(mixture, "r") as mix:
						timepoints = mix.readline().split(",")[1:]
						output.write(tumour_name+"\t"+ "\t".join(timepoints))


if __name__ == "__main__":
	mixture_dir_vaf = "/mnt/raisin/yulia/pwgs/samples/vaf/results_onlyKnownSignatures_vaf/BRCA"
	mixture_dir_phi = "/mnt/raisin/yulia/pwgs/samples/psub/results_onlyKnownSignatures_psub/BRCA"
	output_vaf = "/home/ryogali/vaf.txt"
	output_phi = "/home/ryogali/phi.txt"
	# print read_mixture(mixture_dir_vaf, output_vaf)