#!/usr/bin/python

import math, sys, os
import numpy as np
import random

if len(sys.argv) < 3:
	print("Error in input.\nUsage: ./generate_bootstrap_alignments.py [alignment.mafft] [bootstrap nr]")
	sys.exit(1)
else:
	alignment_in = sys.argv[1]
	bootstrap_number = int(sys.argv[2])

#Put labels and sequences from alignment_in in dictionary
seqd = dict()
with open(alignment_in, 'r') as fin:
	for line in fin:
		if ">" in line:
			current_key = line[1:-1]
			current_key = current_key.partition(" ")[0]
			current_key = current_key.partition("/")[0]
			if current_key in seqd.keys():   #Two genes have the same name!
				print("Error: Two genes appear to have the same name: "+current_key)
				sys.exit(1)
		else:
			if current_key in seqd.keys():
				seqd.update({current_key:seqd[current_key]+line[:-1]})
			else:
				seqd.update({current_key:line[:-1]})

#Make a bootstrapped dictionary, randomly drawing columns from the alignment (with replacement)
seq_length = len(seqd[list(seqd.keys())[0]])
random.seed(1000+bootstrap_number)
seqd_bootstrap = dict()
for p in range(seq_length):
	draw = random.randint(0,seq_length-1)
	for key in seqd.keys():
		if p == 0:
			seqd_bootstrap.update({key:seqd[key][draw]})
		else:
			seqd_bootstrap.update({key:seqd_bootstrap[key]+seqd[key][draw]})

#Output each bootstrapped dictionary
with open(alignment_in+"_bs:"+str(bootstrap_number), 'a') as fout:
	for bootkey in sorted(seqd_bootstrap.keys()):
		fout.write(">"+bootkey+"\n"+seqd_bootstrap[bootkey]+"\n")
