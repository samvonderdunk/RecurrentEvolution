#!/usr/bin/python

import math,sys,os

if len(sys.argv) < 3:
 	print("Input error.\nUsage: ./extract_two_copy_genes.py [alignment.mafft] [duplicates.tcg]")
 	sys.exit(1)
else:
	filename_in = sys.argv[1]
	duplicates_output = sys.argv[2]

#Dictionaries which will contain the sequences.
#The idea is that if you encounter a species the first time, the gene goes into seq_dict. Once you encounter the same species again, both sequences are put in the dupl_dict and deleted from seq_dict. If we encounter the same species more than twice, we will put the sequences in multicopy_dict and remove them from dupl_dict.
seq_dict = dict()
dupl_dict = dict()
multicopy_dict = dict()
with open(filename_in, 'r') as fin:
	for line in fin:
		if '>' in line:
			if line[1:5] not in [x[1:5] for x in seq_dict.keys()]:
				if line[1:5] not in [x[1:5] for x in dupl_dict.keys()]:
					# Species is not present in any list:
					if line[1:5] not in [x[1:5] for x in multicopy_dict.keys()]:
						current_key = line
						target_list = 1
					# Species is present in the multiple-copy list:
					else:
						previous_tags = len([key for key in multicopy_dict.keys() if key[-4:-2] == ".c"])
						if previous_tags == 0:
							multicopy_key = line
						else:
							multicopy_key = line[:-1]+".c"+str(previous_tags)+"\n"
						target_list = 3

				# Species is present in the two-copy list:
				else:
					del_key = [key for key in dupl_dict.keys() if line[1:5] in key[1:5]]
					multicopy_dict.update({del_key[0]:dupl_dict[del_key[0]]})
					multicopy_dict.update({del_key[1]:dupl_dict[del_key[1]]})
					if del_key[1][-3:-1] == ".c":
						multicopy_key = line[:-1]+".c3\n"
					else:
						multicopy_key = line

					del dupl_dict[del_key[0]]
					del dupl_dict[del_key[1]]
					target_list = 3

			# Species is present in the single-copy list:
			else:
				del_key = [key for key in seq_dict.keys() if line[1:5] in key[1:5]]
				if del_key[0]  == line:
					dupl_dict.update({del_key[0][:-1]+".c1\n":seq_dict[del_key[0]]})
					dupl_key = line[:-1]+".c2\n"
				else:
					dupl_dict.update({del_key[0]:seq_dict[del_key[0]]})
					dupl_key = line

				del seq_dict[del_key[0]]#[1:5]]
				target_list = 2
		else:
			if target_list == 1:
				if current_key[1:5] not in [x[1:5] for x in seq_dict.keys()]:	# When you add a new entry to the dictionary
					seq_dict.update({current_key:line})
				else:	# When the key is already there, you append the sequence entry (in case the sequence spans multiple lines)
					seq_dict.update({current_key:seq_dict[current_key]+line})
			elif target_list == 2:
				if dupl_key not in dupl_dict.keys():	# Same as above
					dupl_dict.update({dupl_key:line})
				else:
					dupl_dict.update({dupl_key:dupl_dict[dupl_key]+line})
			else:
				if multicopy_key not in multicopy_dict.keys():	# Same as above
					multicopy_dict.update({multicopy_key:line})
				else:
					multicopy_dict.update({multicopy_key:multicopy_dict[multicopy_key]+line})

with open(duplicates_output, 'a') as fout:
	for k, s in sorted(dupl_dict.items()):
		fout.write(k+s)
