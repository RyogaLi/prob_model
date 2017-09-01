#!/usr/bin/python
# coding=utf-8
#v2.0
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import os
import fnmatch

def convert_file(filename):
	matrix = []
	with open(filename, "r") as inputfile:
		for line in inputfile:
			line = line.split("=")
			line[1] = float(line[1].strip("/n"))
			matrix.append(line)
	matrix = np.asarray(matrix)
	return matrix

def compare_vaf(VAF, Phi):
	VAF_matrix = convert_file(VAF)
	phi_matrix = convert_file(Phi)

	# sort by vaf

	vaf_y = [float(i) for i in VAF_matrix[VAF_matrix[:,1].argsort()][:,1].tolist()]
	# vaf_x = range(1,138)
	# vaf_y = np.arange(1,(len(vaf_x)/100)+0.01,0.01).tolist()
	vaf_x = range(0,len(vaf_y))
	print len(vaf_x)
	print len(vaf_y)
	phi_y = [float(i) for i in phi_matrix[phi_matrix[:,1].argsort()][:,1].tolist()]
	# phi_y = np.arange(1,len(phi_x)+0.01,0.01).tolist()
	phi_x = range(0,len(phi_y))

	plt.plot(vaf_x, vaf_y, label='label here', color="red")
	plt.plot(phi_x, phi_y, label='label here', color="blue")
	plt.title('vaf and phi values')
	plt.ylabel('')
	plt.xlabel('')
	plt.savefig(VAF+".png") 
	plt.close()

def find(pattern, path):
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result = os.path.join(root, name)
    return result

if __name__ == '__main__':
	VAF_dir = "/Users/RyogaLi/Desktop/phi/ssmvaf/"
	phi_dir = "/Users/RyogaLi/Desktop/phi/consprelim.sample-5000.psub/"
	for phi_file in os.listdir(phi_dir):
		tumour_name = phi_file.split(".")[0]

		phi = phi_dir+phi_file
		vaf = find(tumour_name+"*", VAF_dir)

		compare_vaf(vaf, phi)
