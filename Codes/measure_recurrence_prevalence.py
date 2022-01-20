#!/usr/bin/python

import math,sys,os
import numpy as np
from decimal import *

#Decimal precision
getcontext().prec = 7

if len(sys.argv) < 4:
	print("Input error.\nUsage: ./measure_recurrence_prevalence.py [groups.indies] [groups.fates] [output.score]")
	sys.exit(1)
else:
	indie_clusters = sys.argv[1]
	fate_clusters = sys.argv[2]
	score_output = sys.argv[3]

def GROUPVALUE(N):
	return 1.0

def FATESCORE(list_in):
	fatescore = 0.0
	save_scores = []
	for i in set(list_in):
		fatescore += GROUPVALUE(list_in.count(i))
		save_scores.append(GROUPVALUE(list_in.count(i)))
	return fatescore

#Import cluster data
with open(indie_clusters, 'r') as fin:
	IndieClusters = []
	next(fin)
	for line in fin:
		if len(line)>3:
			species=line.split("\t")
			cluster = []
			if species[0][0] == "0" or species[0][0] == "1":
				start_list = 1
			else:
				start_list = 0
			for s in species[start_list:]:
				if s != "\n":
					cluster.append(s[:4])
			IndieClusters.append(cluster)

with open(fate_clusters, 'r') as fin:
	FateClusters = []
	next(fin)
	for line in fin:
		if len(line)>3:
			fate=line.split("\t")
			cluster = []
			if fate[0][0] == "0" or fate[0][0] == "1":
				start_list = 1
			else:
				start_list = 0
			for f in fate[start_list:]:
				if f != "\n":
					c = f.strip(" ")
					cluster.append(c.strip("\n"))
			FateClusters.append(cluster)

if len(FateClusters) < 2:
	with open(score_output, 'a') as fout:
		fout.write("P = 0.0\n")
		sys.exit(1)

#Check whether fate clustering is symmetric
spoofFateClusters = []
for Fate in FateClusters:
	spoofFateClusters.append([F[:4] for F in Fate])
for i in range(0, len(spoofFateClusters)):
	spFate1 = spoofFateClusters[i]
	validated = 0
	for j in range(0, len(spoofFateClusters)):
		if j != i:
			spFate2 = spoofFateClusters[j]
			if sorted(spFate1) == sorted(spFate2):
				validated = 1
	if validated == 0:
		print("Warning: fate clusters are not symmetric.\n")
		with open(score_output, 'a') as fout:
			fout.write("P = NA\n")
		sys.exit(1)

#Calculate score
lfatescores = []
for Fate in FateClusters:
	Indie_IDs = []
	for Gene in Fate:
		for i in range(0, len(IndieClusters)):
			if Gene[:4] in IndieClusters[i]:
				Indie_IDs.append(i)
	#Now check that only complete indies are in my fate (not just the number of unique Indie_IDs)
	Indie_IDs_complete = Indie_IDs
	for ID in set(Indie_IDs):
		number_ids = len([x for x in Indie_IDs if x==ID])
		if number_ids > len(IndieClusters[ID]):
			print("Error: indie does not contain so many species.\n")
			with open(score_output, 'a') as fout:
				fout.write("P = NA\n")
			sys.exit(1)
		elif number_ids < len(IndieClusters[ID]):
			Indie_IDs_complete = [x for x in Indie_IDs_complete if x!=ID]
	if len(Indie_IDs_complete) > 0:
		lfatescores.append(FATESCORE(Indie_IDs_complete))
	else:
		lfatescores.append(0.0)

TotalScore = max(lfatescores)

Indie_Splits = 0
for Indie in IndieClusters:
	Fate_IDs = []
	for Species in Indie:
		for i in range(0, len(FateClusters)):
			if Species in [f[:4] for f in FateClusters[i]]:
				Fate_IDs.append(i)
	if len(set(Fate_IDs)) > 2:
		Indie_Splits += (len(set(Fate_IDs))-2)/2

count_fates = 0
with open(score_output, 'a') as fout:
	for s, F in zip(lfatescores, FateClusters):
		if s == TotalScore:
			count_fates += 1
			fout.write("Fate "+str(count_fates)+": "+"\t".join(F)+"\n")
	fout.write("P = "+str(TotalScore)+"\n")
