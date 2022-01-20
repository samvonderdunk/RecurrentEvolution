#!/usr/bin/python

import math, sys, os
import numpy as np

if len(sys.argv) < 4:
	print("Input error.\nUsage: ./cluster_bootstraps.py [output.fates] [originals.fcl] [[bootstraps.fcl]]")
else:
	jaccard_output = sys.argv[1]		#Output for the merged clusters
	real_set = sys.argv[2]				#Clusters from the original data
	comparison_sets = sys.argv[3:]		#Clusters from each of the bootstraps

# Read clusters from original data
with open(real_set, 'r') as fin:
	IndieClusters = []
	for line in fin:
		if len(line)>3:
			species=line.split("\t")
			cluster = []
			for s in species:
				if s != "\n":
					sp = s.strip(" ")
					cluster.append(sp.strip("\n"))
			IndieClusters.append(cluster)

# Read clusters from all bootstrapped data
IndieClusters_comp = [[] for n in range(len(comparison_sets))]
empty_cluster_files = 0
for n in range(len(comparison_sets)):
	with open(comparison_sets[n], 'r') as fin:
		for line in fin:
			if len(line)>3:
				species=line.split("\t")
				cluster = []
				for s in species:
					if s != "\n":
						sp = s.strip(" ")
						cluster.append(sp.strip("\n"))
				IndieClusters_comp[n].append(cluster)

# Find all unique clusters in the ensemble of bootstraps and the original data
AllUniqueClusters = [sorted(cluster) for cluster in IndieClusters]
for bootstrap in IndieClusters_comp:
	for cluster in bootstrap:
		if sorted(cluster) not in AllUniqueClusters:
			AllUniqueClusters.append(sorted(cluster))

# Define Jaccard overlap between two sets
def Jaccard(set_A, set_B):
	return len(set_A.intersection(set_B))/float(len(set_A.union(set_B)))

# Calculate the jaccard congruence of each unique cluster over all bootstrap clusterings.
# In other words, for each unique cluster find its best matching cluster in every bootstrap clustering, and average these jaccard scores.
# Last line makes sure to only include values higher than zero, because a zero indicates that the genes of the unique clusters were absent from the particular clustering file. Generally this means that a clustering file was completely empty, in which case the first line of the output file should inform the user (i.e. "#Bootstraps = 99").
IndieJaccards = []
for cluster in AllUniqueClusters:
	IndieJacs = []
	for bootstrap in IndieClusters_comp:
		max_jaccard = 0.0
		for try_comparison_cluster in bootstrap:		#Compare the unique cluster to each bootstrap clustering
			jac = Jaccard(set(cluster), set(try_comparison_cluster))
			if jac > max_jaccard:
				max_jaccard = jac
		IndieJacs.append(max_jaccard)
	max_jaccard = 0.0
	for try_comparison_cluster in IndieClusters:		#Compare the unique cluster to the original clustering
		jac = Jaccard(set(cluster), set(try_comparison_cluster))
		if jac > max_jaccard:
			max_jaccard = jac
	IndieJacs.append(max_jaccard)
	del IndieJacs[IndieJacs.index(1.0)]
	IndieJaccards.append(round(np.average([x for x in IndieJacs if x != 0.0]),3))

# Simple heuristic to pick the unique clusters starting with the one with the highest jaccard congruence, and then picking the next best cluster that does not have any genes that are already in the first cluster, etcetera.
JaccardSorted_Clusters = [x for _,x in sorted(zip(IndieJaccards, AllUniqueClusters))]
JaccardSorted_Clusters.reverse()
ChooseIndies = []
ChooseJaccard = []
UsedSpecies = []
for cluster in JaccardSorted_Clusters:
	species_double = 0
	for species in cluster:
		if species in UsedSpecies:
			species_double = 1
			break
	if species_double == 0:
		ChooseIndies.append(cluster)
		ChooseJaccard.append(sorted(IndieJaccards)[len(JaccardSorted_Clusters)-JaccardSorted_Clusters.index(cluster)-1])
		UsedSpecies.extend([species for species in cluster])

# Below are commented out some outputs that can give the user an idea of the congruence of the clusters throughout bootstraps

# Print all clusters found anywhere along with there averaged jaccard coefficients
# Note that sometimes the species can be the same, but that two paralogs of a species have swapped places in clusters
"""
for x,y in zip(JaccardSorted_Clusters[::2], reversed(sorted(IndieJaccards)[::2])):
 	print str(y)+"\t"+str(x)

# Print for each cluster found only in the real data, its members and its jaccard coefficient
for x,y in zip(JaccardSorted_Clusters, reversed(sorted(IndieJaccards))):
 	if x in IndieClusters:
 		print str(y)+"\t"+str(x)

# Print in a greedy way the best clusters according to their jaccard coefficients. This is the approach we use in the pipeline.
for x,y in zip(ChooseIndies, ChooseJaccard):
 	print str(y)+"\t"+str(x)
"""

with open(jaccard_output, 'w') as fout:
	fout.write("#Bootstraps = "+str(len(comparison_sets)-empty_cluster_files)+"\n")
	for x,y in zip(ChooseIndies, ChooseJaccard):
		fout.write(str(y)+"\t"+"\t".join([z for z in x])+"\n")
