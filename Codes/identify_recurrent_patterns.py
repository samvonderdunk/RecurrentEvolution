#!/usr/bin/python
#############################################################################################################################################################################################################
import math,sys,os
import numpy as np
import matplotlib.pyplot as plt
from decimal import *
from random import shuffle
import subprocess
getcontext().prec = 7	#Decimal precision, don't remember the reason for setting this
#############################################################################################################################################################################################################

#############################################################################################################################################################################################################
#My functions
def GROUPVALUE(N):
	return 1.0

def FATESCORE(list_in):
	fatescore = 0.0
	save_scores = []
	for i in set(list_in):
		fatescore += GROUPVALUE(list_in.count(i))
		save_scores.append(GROUPVALUE(list_in.count(i)))
	return fatescore

def Shannon_Entropy(list_in):	#Information Content
	H = 0
	if sum(list_in) != 1.0:
		temp_list = [None] * len(list_in)
		for x in range(0,len(list_in)):
			temp_list[x] = list_in[x]/float(sum(list_in))
		list_in = [round(x,3) for x in temp_list]
	for x in range(0,len(list_in)):
		if list_in[x] != 0.0:
			H += list_in[x]*math.log(list_in[x]/(1./len(list_in)),2)
	return -H
#############################################################################################################################################################################################################

#############################################################################################################################################################################################################
#Receive input-arguments
if len(sys.argv) < 4:
	print("Input error.\nUsage: ./identify_recurrent_patterns.py [alignment.tcg] [groups.indies] [groups.fates] [poscol.jalview] [groups.jalview] [alignment.aln] [seqlogo.png]")
	sys.exit(1)
else:
	alignment = sys.argv[1]
	indie_clusters = sys.argv[2]
	fate_clusters = sys.argv[3]

position_coloring = sys.argv[4]
groups_file = sys.argv[5]
alignment_out = sys.argv[6]
save_logo = sys.argv[7]
#############################################################################################################################################################################################################

#############################################################################################################################################################################################################
#############################################################################################################################################################################################################
#Put labels and sequences in dictionary (only the ones that are in the highest-scoring fate clusters, see below)
seqd = dict()
with open(alignment, 'r') as fin:
	for line in fin:
		if ">" in line:
			current_key = line[1:-1]
			current_key = current_key.partition(" ")[0]
			current_key = current_key.partition("/")[0]
		else:
			if current_key in seqd.keys():
				seqd.update({current_key:seqd[current_key]+line[:-1]})
			else:
				seqd.update({current_key:line[:-1]})

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

#############################################################################################################################################################################################################

#############################################################################################################################################################################################################
#Determine what the best-scoring convergent fate is
lfatescores = []
lfate_indies = []
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
		elif number_ids < len(IndieClusters[ID]):
			Indie_IDs_complete = [x for x in Indie_IDs_complete if x!=ID]
	lfatescores.append(FATESCORE(Indie_IDs_complete))
	lfate_indies.append(Indie_IDs_complete)

count = 0
MaxClusters = []
for fs in lfatescores:
	if fs == max(lfatescores):
		MaxClusters.append(FateClusters[count])
		MaxClusters_Indies = lfate_indies[count]
	count += 1

Included_IndieClusters = []
for i in range(0,len(IndieClusters)):
	if i in set(MaxClusters_Indies):
		Included_IndieClusters.append(IndieClusters[i])

#Remove genes that are not in the highest-scoring fates
name_dict = dict()
for key in list(seqd.keys()):
	if key not in MaxClusters[0] and key not in MaxClusters[1]:
		sequence = seqd[key]
		del seqd[key]
	elif key in MaxClusters[0]:
		sequence = seqd[key]
		del seqd[key]
		seqd.update({key[:4]+"_A":sequence})
		name_dict.update({key[:4]+"_A":key})
	elif key in MaxClusters[1]:
		sequence = seqd[key]
		del seqd[key]
		seqd.update({key[:4]+"_B":sequence})
		name_dict.update({key[:4]+"_B":key})
#############################################################################################################################################################################################################

#############################################################################################################################################################################################################
#Making an alignment-file of the duplicates that show fate recurrence, duplicates ordered as fate-A, fate-B, etc. and species grouped (sorted) per indie.
with open(alignment_out,'a') as fout:
	for Indie in Included_IndieClusters:
		names = [x for x in seqd.keys() if x[:4] in Indie]
		for seq_name in sorted(names):
			fout.write(">"+name_dict[seq_name]+"\n"+seqd[seq_name]+"\n")
#############################################################################################################################################################################################################

#############################################################################################################################################################################################################
#Making the grouping file for jalview
with open(groups_file, 'a') as fout:
	fout.write("JALVIEW_ANNOTATION\n")
	count_seqs = 0
	for Indie in Included_IndieClusters:
		number_of_seqs = len([x for x in seqd.keys() if x[:4] in Indie])
		fout.write("SEQUENCE_GROUP\tGROUP_"+str(Included_IndieClusters.index(Indie))+"\t1\t"+str(len(sequence))+"\t"+str(count_seqs+1)+"-"+str(count_seqs+number_of_seqs)+"\n")
		count_seqs += number_of_seqs
	for Indie in Included_IndieClusters:
		fout.write("PROPERTIES\tGROUP_"+str(Included_IndieClusters.index(Indie))+"\toutlineColour=gray\tdisplayBoxes=true\n")
#############################################################################################################################################################################################################

#############################################################################################################################################################################################################
#Blending with white means that the 'zeros' in the rgb-arrays take on these values: 0, 26, 51, 77, 102, 128, 153, 179, 204, 230, 255 (https://meyerweb.com/eric/tools/color-blend/#FFFFFF:00FFFF:9:rgbd)
with open(position_coloring, 'a') as fout:

	#From cyan to green or blue
	fout.write("cFate_m1.0_0.9_AxxA\t0,255,255|0,0,255|absolute|0.0|1.0\n")
	fout.write("cFate_m1.0_0.9_xAAx\t0,255,255|0,255,0|absolute|0.0|1.0\n")
	fout.write("cFate_m0.9_0.8_AxxA\t26,255,255|26,26,255|absolute|0.0|0.9\n")
	fout.write("cFate_m0.9_0.8_xAAx\t26,255,255|26,255,26|absolute|0.0|0.9\n")
	fout.write("cFate_m0.8_0.7_AxxA\t51,255,255|51,51,255|absolute|0.0|0.8\n")
	fout.write("cFate_m0.8_0.7_xAAx\t51,255,255|51,255,51|absolute|0.0|0.8\n")
	fout.write("cFate_m0.7_0.6_AxxA\t77,255,255|77,77,255|absolute|0.0|0.7\n")
	fout.write("cFate_m0.7_0.6_xAAx\t77,255,255|77,255,77|absolute|0.0|0.7\n")
	fout.write("cFate_m0.6_0.5_AxxA\t102,255,255|102,102,255|absolute|0.0|0.6\n")
	fout.write("cFate_m0.6_0.5_xAAx\t102,255,255|102,255,102|absolute|0.0|0.6\n")
	fout.write("cFate_m0.5_0.4_AxxA\t128,255,255|128,128,255|absolute|0.0|0.5\n")
	fout.write("cFate_m0.5_0.4_xAAx\t128,255,255|128,255,128|absolute|0.0|0.5\n")
	fout.write("cFate_m0.4_0.3_AxxA\t153,255,255|153,153,255|absolute|0.0|0.4\n")
	fout.write("cFate_m0.4_0.3_xAAx\t153,255,255|153,255,153|absolute|0.0|0.4\n")
	fout.write("cFate_m0.3_0.2_AxxA\t179,255,255|179,179,255|absolute|0.0|0.3\n")
	fout.write("cFate_m0.3_0.2_xAAx\t179,255,255|179,255,179|absolute|0.0|0.3\n")
	fout.write("cFate_m0.2_0.1_AxxA\t204,255,255|204,204,255|absolute|0.0|0.2\n")
	fout.write("cFate_m0.2_0.1_xAAx\t204,255,255|204,255,204|absolute|0.0|0.2\n")
	fout.write("cFate_m0.1_0.0_AxxA\t255,255,255|230,230,255|absolute|0.0|0.1\n")
	fout.write("cFate_m0.1_0.0_xAAx\t255,255,255|230,255,230|absolute|0.0|0.1\n")

	#From red to magenta or yellow (check visibility)
	fout.write("Fate_p0.0_0.1_AxAx\t255,255,255|255,230,255|absolute|0.0|0.1\n")
	fout.write("Fate_p0.0_0.1_xAxA\t255,255,255|255,255,230|absolute|0.0|0.1\n")
	fout.write("Fate_p0.1_0.2_AxAx\t255,204,204|255,204,255|absolute|0.0|0.2\n")
	fout.write("Fate_p0.1_0.2_xAxA\t255,204,204|255,255,204|absolute|0.0|0.2\n")
	fout.write("Fate_p0.2_0.3_AxAx\t255,179,179|255,179,255|absolute|0.0|0.3\n")
	fout.write("Fate_p0.2_0.3_xAxA\t255,179,179|255,255,179|absolute|0.0|0.3\n")
	fout.write("Fate_p0.3_0.4_AxAx\t255,153,153|255,153,255|absolute|0.0|0.4\n")
	fout.write("Fate_p0.3_0.4_xAxA\t255,153,153|255,255,153|absolute|0.0|0.4\n")
	fout.write("Fate_p0.4_0.5_AxAx\t255,128,128|255,128,255|absolute|0.0|0.5\n")
	fout.write("Fate_p0.4_0.5_xAxA\t255,128,128|255,255,128|absolute|0.0|0.5\n")
	fout.write("Fate_p0.5_0.6_AxAx\t255,102,102|255,102,255|absolute|0.0|0.6\n")
	fout.write("Fate_p0.5_0.6_xAxA\t255,102,102|255,255,102|absolute|0.0|0.6\n")
	fout.write("Fate_p0.6_0.7_AxAx\t255,77,77|255,77,255|absolute|0.0|0.7\n")
	fout.write("Fate_p0.6_0.7_xAxA\t255,77,77|255,255,77|absolute|0.0|0.7\n")
	fout.write("Fate_p0.7_0.8_AxAx\t255,51,51|255,51,255|absolute|0.0|0.8\n")
	fout.write("Fate_p0.7_0.8_xAxA\t255,51,51|255,255,51|absolute|0.0|0.8\n")
	fout.write("Fate_p0.8_0.9_AxAx\t255,26,26|255,26,255|absolute|0.0|0.9\n")
	fout.write("Fate_p0.8_0.9_xAxA\t255,26,26|255,255,26|absolute|0.0|0.9\n")
	fout.write("Fate_p0.9_1.0_AxAx\t255,0,0|255,0,255|absolute|0.0|1.0\n")
	fout.write("Fate_p0.9_1.0_xAxA\t255,0,0|255,255,0|absolute|0.0|1.0\n")
#############################################################################################################################################################################################################

#############################################################################################################################################################################################################
#############################################################################################################################################################################################################
lindie_averaged_position_score = [[None for a in range(7)] for b in range(len(sequence) * len(Included_IndieClusters))]
for indie_x in range(0,len(Included_IndieClusters)):
	countx = 0
	lcolumns_position_score = [[None for a in range(7)] for b in range(len(sequence) * len(Included_IndieClusters[indie_x]))]
	#Start doing species X vs species Y.
	for x in sorted(seqd.keys()):
		indie_counter_x = 0
		for Indie in Included_IndieClusters:
			if x[:4] in Indie:
				break
			indie_counter_x += 1
		if indie_counter_x != indie_x:
			continue

		lindie_position_score = [[None for a in range(7)] for b in range(len(sequence) * (len(Included_IndieClusters)-1))]
		#What indie does the focal species belong to
		indie_counter = 0
		for Indie in Included_IndieClusters:
			if x[:4] in Indie:
				IndieX = indie_counter
			indie_counter += 1
		countx += 1
		#We can skip every other X, because we do comparisons per species not per duplicate.
		if countx % 2 == 0:
			continue
		else:
			county = 0
			compare_indie_counter = 0
			for IndieY in range(0,len(Included_IndieClusters)):
				lposition_score = [None] * len(sequence) * len(Included_IndieClusters[IndieY])	#Note that sequence should still be one of the sequences as it was set in the loop that adjusts the dictionary based on the highest-scoring fates, see above
				indie_member_counter = 0
				if IndieY == IndieX:
					continue
				else:
					for y in sorted(seqd.keys()):
						#Treat other species from the same duplication special (not or differently)
						indie_counter = 0
						for Indie in Included_IndieClusters:
							if y[:4] in Indie:
								break
							indie_counter += 1
						county += 1
						#Do only the species that belong to the indie that we are now comparing to the focal species
						if indie_counter != IndieY:
							continue
						#Same for Y; skip redundant comparisons.
						elif county % 2 == 0:
							continue
						else:
							#Now find the duplicate of our current x and y.
							x1_label = x
							x1 = seqd[x]
							x2_label = [n for n in seqd.keys() if (n[:4]==x[:4] and n!=x)][0]
							x2 = seqd[x2_label]
							y1_label = y
							y1 = seqd[y]
							y2_label = [n for n in seqd.keys() if (n[:4]==y[:4] and n!=y)][0]
							y2 = seqd[y2_label]

							position_counter = 0
							for px1, px2, px3, px4 in zip(x1,x2,y1,y2):
								if px1 == '-':
									if px2 == '-':
										if px3 == '-':
											if px4 == '-':									# ---- #
												Fate = 'OOOO'
											elif px4 != '-':								# ---A #
												Fate = 'OOOA'

										elif px3 != '-':
											if px4 == '-':									# --A- #
												Fate = 'OOAO'
											elif px4 != '-':
												if px3 == px4:								# --AA #
													Fate = 'OOAA'
												elif px3 != px4:							# --AB #
													Fate = 'OOAB'
									elif px2 != '-':
										if px3 == '-':
											if px4 == '-':									# -A-- #
												Fate = 'OAOO'
											elif px4 != '-':
												if px2 == px4:								# -A-A #
													Fate = 'OAOA'
												elif px2 != px4:							# -A-B #
													Fate = 'OAOB'
										elif px3 != '-':
											if px4 == '-':
												if px2 == px3:								# -AA- #
													Fate = 'OAAO'
												elif px2 != px3:							# -AB- #
													Fate = 'OABO'
											elif px4 != '-':							## -xxx ##
												if px2 == px3:
													if px2 == px4:							# -AAA #
														Fate = 'OAAA'
													elif px2 != px4:						# -AAB #
														Fate = 'OAAB'
												elif px2 != px3:
													if px2 == px4:							# -ABA #
														Fate = 'OABA'
													elif px2 != px4:
														if px3 == px4:						# -ABB #
															Fate = 'OABB'
														elif px3 != px4:					# -ABC #
															Fate = 'OABC'
								elif px1 != '-':
									if px2 == '-':
										if px3 == '-':
											if px4 == '-':
												Fate = 'AOOO'								# A--- #
											elif px4 != '-':
												if px1 == px4:								# A--A #
													Fate = 'AOOA'
												elif px1 != px4:							# A--B #
													Fate = 'AOOB'
										elif px3 != '-':
											if px4 == '-':
												if px1 == px3:								# A-A- #
													Fate = 'AOAO'
												elif px1 != px3:							# A-B- #
													Fate = 'AOBO'
											elif px4 != '-':							## x-xx ##
												if px1 == px3:
													if px1 == px4:							# A-AA #
														Fate = 'AOAA'
													elif px1 != px4:						# A-AB #
														Fate = 'AOAB'
												elif px1 != px3:
													if px1 == px4:							# A-BA #
														Fate = 'AOBA'
													elif px1 != px4:
														if px3 == px4:						# A-BB #
															Fate = 'AOBB'
														elif px3 != px4:					# A-BC #
															Fate = 'AOBC'
									elif px2 != '-':
										if px3 == '-':
											if px4 == '-':
												if px1 == px2:								# AA-- #
													Fate = 'AAOO'
												elif px1 != px2:							# AB-- #
													Fate = 'ABOO'
											elif px4 != '-':							## xx-x ##
												if px1 == px2:
													if px1 == px4:							# AA-A #
														Fate = 'AAOA'
													elif px1 != px4:						# AA-B #
														Fate = 'AAOB'
												elif px1 != px2:
													if px1 == px4:							# AB-A #
														Fate = 'ABOA'
													elif px1 != px4:
														if px2 == px4:						# AB-B #
															Fate = 'ABOB'
														elif px2 != px4:					# AB-C #
															Fate = 'ABOC'
										elif px3 != '-':
											if px4 == '-':								## xxx- ##
												if px1 == px2:
													if px1 == px3:							# AAA- #
														Fate = 'AAAO'
													elif px1 != px3:						# AAB- #
														Fate = 'AABO'
												elif px1 != px2:
													if px1 == px3:							# ABA- #
														Fate = 'ABAO'
													elif px1 != px3:
														if px2 == px3:						# ABB- #
															Fate = 'ABBO'
														elif px2 != px3:					# ABC- #
															Fate = 'ABCO'
											elif px4 != '-':							## xxxx ##
												if px1 == px2:
													if px1 == px3:
														if px1 == px4:						# AAAA #
															Fate = 'AAAA'
														elif px1 != px4:					# AAAB #
															Fate = 'AAAB'
													elif px1 != px3:
														if px1 == px4:						# AABA #
															Fate = 'AABA'
														elif px1 != px4:
															if px3 == px4:					# AABB #
																Fate = 'AABB'
															elif px3 != px4:				# AABC #
																Fate = 'AABC'
												elif px1 != px2:
													if px1 == px3:
														if px1 == px4:						# ABAA #
															Fate = 'ABAA'
														elif px1 != px4:
															if px2 == px4:					# ABAB #
																Fate = 'ABAB'
															elif px2 != px4:				# ABAC #
																Fate = 'ABAC'
													elif px1 != px3:
														if px1 == px4:
															if px2 == px3:					# ABBA #
																Fate = 'ABBA'
															if px2 != px3:					# ABCA #
																Fate = 'ABCA'
														elif px1 != px4:
															if px2 == px3:
																if px2 == px4:				# ABBB #
																	Fate = 'ABBB'
																elif px2 != px4:			# ABBC #
																	Fate = 'ABBC'
															elif px2 != px3:
																if px2 == px4:				# ABCB #
																	Fate = 'ABCB'
																elif px2 != px4:
																	if px3 == px4:			# ABCC #
																		Fate = 'ABCC'
																	elif px3 != px4:		# ABCD #
																		Fate = 'ABCD'
								Fate_AxAx = ['AOAO','AOAB','ABAO','ABAC','AOBO']
								Fate_xAxA = ['OAOA','OABA','ABOB','ABCB','OAOB']
								Fate_AxxA = ['AOOA','AOBA','ABOA','ABCA','AOOB']
								Fate_xAAx = ['OAAO','OAAB','ABBO','ABBC','OABO']
								if Fate == 'ABAB':
									lposition_score[indie_member_counter*len(sequence)+position_counter] = 5
								elif Fate == 'ABBA':
									lposition_score[indie_member_counter*len(sequence)+position_counter] = 6
								elif Fate in Fate_AxAx:
									lposition_score[indie_member_counter*len(sequence)+position_counter] = 1
								elif Fate in Fate_xAxA:
									lposition_score[indie_member_counter*len(sequence)+position_counter] = 2
								elif Fate in Fate_AxxA:
									lposition_score[indie_member_counter*len(sequence)+position_counter] = 3
								elif Fate in Fate_xAAx:
									lposition_score[indie_member_counter*len(sequence)+position_counter] = 4
								else:
									lposition_score[indie_member_counter*len(sequence)+position_counter] = 0

								position_counter += 1
							indie_member_counter += 1

				for i in range(0,len(sequence)):
					for j in range(0,7):
						lindie_position_score[compare_indie_counter*len(sequence)+i][j] = lposition_score[i::len(sequence)].count(j)/float(len(lposition_score[i::len(sequence)]))
				compare_indie_counter += 1

		for i in range(0,len(sequence)):
			for j in range(0,7):
				sum_this = 0
				for k in range(0,len(Included_IndieClusters)-1):
					sum_this += lindie_position_score[k*len(sequence)+i][j]
				lcolumns_position_score[(int((countx+1)/2)-1)*len(sequence)+i][j] = sum_this/float(len(Included_IndieClusters)-1)
#############################################################################################################################################################################################################

#############################################################################################################################################################################################################
		count_non_gaps = 0
		for i in range(0,len(sequence)):

			pRest = round(lcolumns_position_score[(int((countx+1)/2)-1)*len(sequence)+i][0],3)
			pAxAx = round(lcolumns_position_score[(int((countx+1)/2)-1)*len(sequence)+i][1],3)
			pxAxA = round(lcolumns_position_score[(int((countx+1)/2)-1)*len(sequence)+i][2],3)
			pAxxA = round(lcolumns_position_score[(int((countx+1)/2)-1)*len(sequence)+i][3],3)
			pxAAx = round(lcolumns_position_score[(int((countx+1)/2)-1)*len(sequence)+i][4],3)
			pABAB = round(lcolumns_position_score[(int((countx+1)/2)-1)*len(sequence)+i][5],3)
			pABBA = round(lcolumns_position_score[(int((countx+1)/2)-1)*len(sequence)+i][6],3)
			Fate = pABAB + pAxAx + pxAxA - pABBA - pAxxA - pxAAx

			if x1[i] != '-':
				count_non_gaps += 1
				if Fate > 0:
					if Fate == 1.0:
						Left_Border = "0.9"
						Right_Border = "1.0"
					else:
						Left_Border = str(float(str(round(Fate,3))[:3]))	#e.g. 0.4 if Fate = 0.4325
						Right_Border = str(float(str(round(Fate,3))[:3])+0.1)	#e.g. 0.5 if Fate = 0.4325
					if pAxAx >= pxAxA:
						with open(position_coloring, 'a') as fout:
							fout.write("Q\t"+name_dict[x1_label]+"\t-1\t"+str(count_non_gaps)+"\t"+str(count_non_gaps)+"\tFate_p"+Left_Border+"_"+Right_Border+"_AxAx\t"+str(min(abs(pAxAx-pxAxA),abs(Fate)))+"\n")
					elif pAxAx < pxAxA:
						with open(position_coloring, 'a') as fout:
							fout.write("Q\t"+name_dict[x1_label]+"\t-1\t"+str(count_non_gaps)+"\t"+str(count_non_gaps)+"\tFate_p"+Left_Border+"_"+Right_Border+"_xAxA\t"+str(min(abs(pAxAx-pxAxA),abs(Fate)))+"\n")

				elif Fate <= 0:
					Left_Border = str(float(str(round(Fate,3))[1:4])+0.1)	#e.g. 0.5 if Fate = -0.4325
					Right_Border = str(float(str(round(Fate,3))[1:4]))	#e.g. 0.4 if Fate = -0.4325
					if pAxxA >= pxAAx:
						with open(position_coloring, 'a') as fout:
							fout.write("Q\t"+name_dict[x1_label]+"\t-1\t"+str(count_non_gaps)+"\t"+str(count_non_gaps)+"\tcFate_m"+Left_Border+"_"+Right_Border+"_AxxA\t"+str(min(abs(pAxxA-pxAAx),abs(Fate)))+"\n")
					elif pAxxA < pxAAx:
						with open(position_coloring, 'a') as fout:
							fout.write("Q\t"+name_dict[x1_label]+"\t-1\t"+str(count_non_gaps)+"\t"+str(count_non_gaps)+"\tcFate_m"+Left_Border+"_"+Right_Border+"_xAAx\t"+str(min(abs(pAxxA-pxAAx),abs(Fate)))+"\n")
#############################################################################################################################################################################################################
		count_non_gaps = 0
		for i in range(0,len(sequence)):

			pRest = round(lcolumns_position_score[(int((countx+1)/2)-1)*len(sequence)+i][0],3)
			pAxAx = round(lcolumns_position_score[(int((countx+1)/2)-1)*len(sequence)+i][1],3)
			pxAxA = round(lcolumns_position_score[(int((countx+1)/2)-1)*len(sequence)+i][2],3)
			pAxxA = round(lcolumns_position_score[(int((countx+1)/2)-1)*len(sequence)+i][3],3)
			pxAAx = round(lcolumns_position_score[(int((countx+1)/2)-1)*len(sequence)+i][4],3)
			pABAB = round(lcolumns_position_score[(int((countx+1)/2)-1)*len(sequence)+i][5],3)
			pABBA = round(lcolumns_position_score[(int((countx+1)/2)-1)*len(sequence)+i][6],3)
			Fate = round(pABAB + pAxAx + pxAxA - pABBA - pAxxA - pxAAx, 3)

			if x2[i] != '-':
				count_non_gaps += 1
				if Fate > 0:
					if Fate == 1.0:
						Left_Border = "0.9"
						Right_Border = "1.0"
					else:
						Left_Border = str(float(str(round(Fate,3))[:3]))	#e.g. 0.4 if Fate = 0.4325
						Right_Border = str(float(str(round(Fate,3))[:3])+0.1)	#e.g. 0.5 if Fate = 0.4325
					if pAxAx >= pxAxA:
						with open(position_coloring, 'a') as fout:
							fout.write("Q\t"+name_dict[x2_label]+"\t-1\t"+str(count_non_gaps)+"\t"+str(count_non_gaps)+"\tFate_p"+Left_Border+"_"+Right_Border+"_AxAx\t"+str(min(abs(pAxAx-pxAxA),abs(Fate)))+"\n")
					elif pAxAx < pxAxA:
						with open(position_coloring, 'a') as fout:
							fout.write("Q\t"+name_dict[x2_label]+"\t-1\t"+str(count_non_gaps)+"\t"+str(count_non_gaps)+"\tFate_p"+Left_Border+"_"+Right_Border+"_xAxA\t"+str(min(abs(pAxAx-pxAxA),abs(Fate)))+"\n")

				elif Fate <= 0:
					Left_Border = str(float(str(round(Fate,3))[1:4])+0.1)	#e.g. 0.5 if Fate = -0.4325
					Right_Border = str(float(str(round(Fate,3))[1:4]))	#e.g. 0.4 if Fate = -0.4325
					if pAxxA >= pxAAx:
						with open(position_coloring, 'a') as fout:
							fout.write("Q\t"+name_dict[x2_label]+"\t-1\t"+str(count_non_gaps)+"\t"+str(count_non_gaps)+"\tcFate_m"+Left_Border+"_"+Right_Border+"_AxxA\t"+str(min(abs(pAxxA-pxAAx),abs(Fate)))+"\n")
					elif pAxxA < pxAAx:
						with open(position_coloring, 'a') as fout:
							fout.write("Q\t"+name_dict[x2_label]+"\t-1\t"+str(count_non_gaps)+"\t"+str(count_non_gaps)+"\tcFate_m"+Left_Border+"_"+Right_Border+"_xAAx\t"+str(min(abs(pAxxA-pxAAx),abs(Fate)))+"\n")
#############################################################################################################################################################################################################

#############################################################################################################################################################################################################
	#Average position scores per indie
	for i in range(0,len(sequence)):
		for j in range(0,7):
			sum_this = 0
			for k in range(0,len(Included_IndieClusters[indie_x])):
				sum_this += lcolumns_position_score[k*len(sequence)+i][j]
			lindie_averaged_position_score[indie_x*len(sequence)+i][j] = round(sum_this/float(len(Included_IndieClusters[indie_x])),3)
#############################################################################################################################################################################################################

#############################################################################################################################################################################################################
linfo_cont = [[None for a in range(len(sequence))] for b in range(7)]
for i in range(0,len(sequence)):
	p = [None] * 7
	for j in range(0,7):
		sum_this = 0
		for k in range(0,len(Included_IndieClusters)):
			sum_this += lindie_averaged_position_score[k*len(sequence)+i][j]
		p[j] = round(sum_this/float(len(Included_IndieClusters)),3)
	for l in range(0,7):
		linfo_cont[l][i] = round(p[l] * -Shannon_Entropy(p), 3)

lplot_cont = [[] for a in range(7)]
lplot_gaps = []
lplot_sims = []
lplot_conservation = []
lx_axis = []
lx_axis_referenced = []
thresholds = [0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55]
conservation_thresholds = [0.0,0.01,0.05,0.1,0.2]
count_axax = [[0 for b in range(5)] for a in range(10)]
count_xaxa = [[0 for b in range(5)] for a in range(10)]
reference_seq_count = 0
for i in range(0,len(sequence)):
	reference = seqd[list(seqd.keys())[0]]
	if reference[i] != '-':
		reference_seq_count += 1
	gaps_in_duplicates = 0
	for key in seqd.keys():
		dup_seq = seqd[key]
		if dup_seq[i] == '-':
			gaps_in_duplicates += 1
	if float(gaps_in_duplicates)/len(seqd.keys()) > 0.5:
		pass
	else:
		lx_axis.append(i+1)
		lx_axis_referenced.append(reference_seq_count)
		for j in range(7):
			lplot_cont[j].append(linfo_cont[j][i]/2.807)


fig, ax = plt.subplots(1,1,figsize=(15,5))
p1 = ax.bar(range(len(lx_axis)), lplot_cont[2], 1, color='yellow', alpha=1.0, linewidth=0, align="edge")
p7 = ax.bar(range(len(lx_axis)), lplot_cont[5], 1, color='red', alpha=1.0, linewidth=0, bottom=lplot_cont[2], align="edge")
p2 = ax.bar(range(len(lx_axis)), lplot_cont[1], 1, color='magenta', alpha=1.0, linewidth=0, bottom=[x1+x2 for x1, x2 in zip(lplot_cont[2],lplot_cont[5])], align="edge")
p3 = ax.bar(range(len(lx_axis)), lplot_cont[3], 1, color='blue', alpha=1.0, linewidth=0, bottom=[x1+x2+x3 for x1, x2, x3 in zip(lplot_cont[2],lplot_cont[5],lplot_cont[1])], align="edge")
p8 = ax.bar(range(len(lx_axis)), lplot_cont[6], 1, color='cyan', alpha=1.0, linewidth=0, bottom=[x1+x2+x3+x4 for x1, x2, x3, x4 in zip(lplot_cont[2],lplot_cont[5],lplot_cont[1],lplot_cont[3])], align="edge")
p4 = ax.bar(range(len(lx_axis)), lplot_cont[4], 1, color='green', alpha=1.0, linewidth=0, bottom=[x1+x2+x3+x4+x5 for x1, x2, x3, x4, x5 in zip(lplot_cont[2],lplot_cont[5],lplot_cont[1],lplot_cont[3],lplot_cont[6])], align="edge")

ax.set_ylim([0.0,1.0])
plt.yticks([0.0,0.25,0.5,0.75])
ax.set_xlim([0,len(lx_axis)])
ax.yaxis.grid(True)
plt.savefig(save_logo)
#############################################################################################################################################################################################################
