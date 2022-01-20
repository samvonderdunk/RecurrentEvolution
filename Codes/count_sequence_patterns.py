#!/usr/bin/python

import math, sys, os

if len(sys.argv) < 3:
	print("Input error.\nUsage: ./count_sequence_patterns.py [duplicates.tcg] [counts.pts]")
	sys.exit(1)
else:
	alignment = sys.argv[1]
	data_output = sys.argv[2]

#Put labels and sequences in dictionary
seqd = dict()
with open(alignment, 'r') as fin:
	for line in fin:
		if ">" in line:
			current_key = line[1:-1]
			current_key = current_key.partition(" ")[0]
		else:
			if current_key in seqd.keys():
				seqd.update({current_key:seqd[current_key]+line[:-1]})
			else:
				seqd.update({current_key:line[:-1]})

countx = 0
#Start doing species X vs species Y.
for x in sorted(seqd.keys()):
	countx += 1
	#We can skip every other X, because we do comparisons per species not per duplicate.
	if countx % 2 == 0:
		continue
	else:
		county = 0
		for y in sorted(seqd.keys()):
			county += 1
			#Same for Y; skip redundant comparisons.
			if county % 2 == 0:
				continue
			else:
				#X and Y cannot be the same species.
				if y[:4] == x[:4]:
					continue

				#Now find the duplicate of our current x and y.
				else:
					x1_label = x
					x1 = seqd[x]
					x2_label = [n for n in seqd.keys() if (n[:4]==x[:4] and n!=x)][0]
					x2 = seqd[x2_label]
					y1_label = y
					y1 = seqd[y]
					y2_label = [n for n in seqd.keys() if (n[:4]==y[:4] and n!=y)][0]
					y2 = seqd[y2_label]

				#We're counting all possible patterns
				OOOO = 0
				OOOA = 0
				OOAO = 0
				OOAA = 0
				OOAB = 0
				OAOO = 0
				OAOA = 0
				OAOB = 0
				OAAO = 0
				OABO = 0
				OAAA = 0
				OAAB = 0
				OABA = 0
				OABB = 0
				OABC = 0
				AOOO = 0
				AOOA = 0
				AOOB = 0
				AOAO = 0
				AOBO = 0
				AOAA = 0
				AOAB = 0
				AOBA = 0
				AOBB = 0
				AOBC = 0
				AAOO = 0
				ABOO = 0
				AAOA = 0
				AAOB = 0
				ABOA = 0
				ABOB = 0
				ABOC = 0
				AAAO = 0
				AABO = 0
				ABAO = 0
				ABBO = 0
				ABCO = 0
				AAAA = 0
				AAAB = 0
				AABA = 0
				AABB = 0
				AABC = 0
				ABAA = 0
				ABAB = 0
				ABAC = 0
				ABBA = 0
				ABCA = 0
				ABBB = 0
				ABBC = 0
				ABCB = 0
				ABCC = 0
				ABCD = 0


				for px1, px2, px3, px4 in zip(x1,x2,y1,y2):
					if px1 == '-':
						if px2 == '-':
							if px3 == '-':
								if px4 == '-':									# ---- #
									OOOO += 1
								elif px4 != '-':								# ---A #
									OOOA += 1

							elif px3 != '-':
								if px4 == '-':									# --A- #
									OOAO += 1
								elif px4 != '-':
									if px3 == px4:								# --AA #
										OOAA += 1
									elif px3 != px4:							# --AB #
										OOAB += 1
						elif px2 != '-':
							if px3 == '-':
								if px4 == '-':									# -A-- #
									OAOO += 1
								elif px4 != '-':
									if px2 == px4:								# -A-A #
										OAOA += 1
									elif px2 != px4:							# -A-B #
										OAOB += 1
							elif px3 != '-':
								if px4 == '-':
									if px2 == px3:								# -AA- #
										OAAO += 1
									elif px2 != px3:							# -AB- #
										OABO += 1
								elif px4 != '-':							## -xxx ##
									if px2 == px3:
										if px2 == px4:							# -AAA #
											OAAA += 1
										elif px2 != px4:						# -AAB #
											OAAB += 1
									elif px2 != px3:
										if px2 == px4:							# -ABA #
											OABA += 1
										elif px2 != px4:
											if px3 == px4:						# -ABB #
												OABB += 1
											elif px3 != px4:					# -ABC #
												OABC += 1
					elif px1 != '-':
						if px2 == '-':
							if px3 == '-':
								if px4 == '-':
									AOOO += 1								# A--- #
								elif px4 != '-':
									if px1 == px4:								# A--A #
										AOOA += 1
									elif px1 != px4:							# A--B #
										AOOB += 1
							elif px3 != '-':
								if px4 == '-':
									if px1 == px3:								# A-A- #
										AOAO += 1
									elif px1 != px3:							# A-B- #
										AOBO += 1
								elif px4 != '-':							## x-xx ##
									if px1 == px3:
										if px1 == px4:							# A-AA #
											AOAA += 1
										elif px1 != px4:						# A-AB #
											AOAB += 1
									elif px1 != px3:
										if px1 == px4:							# A-BA #
											AOBA += 1
										elif px1 != px4:
											if px3 == px4:						# A-BB #
												AOBB += 1
											elif px3 != px4:					# A-BC #
												AOBC += 1
						elif px2 != '-':
							if px3 == '-':
								if px4 == '-':
									if px1 == px2:								# AA-- #
										AAOO += 1
									elif px1 != px2:							# AB-- #
										ABOO += 1
								elif px4 != '-':							## xx-x ##
									if px1 == px2:
										if px1 == px4:							# AA-A #
											AAOA += 1
										elif px1 != px4:						# AA-B #
											AAOB += 1
									elif px1 != px2:
										if px1 == px4:							# AB-A #
											ABOA += 1
										elif px1 != px4:
											if px2 == px4:						# AB-B #
												ABOB += 1
											elif px2 != px4:					# AB-C #
												ABOC += 1
							elif px3 != '-':
								if px4 == '-':								## xxx- ##
									if px1 == px2:
										if px1 == px3:							# AAA- #
											AAAO += 1
										elif px1 != px3:						# AAB- #
											AABO += 1
									elif px1 != px2:
										if px1 == px3:							# ABA- #
											ABAO += 1
										elif px1 != px3:
											if px2 == px3:						# ABB- #
												ABBO += 1
											elif px2 != px3:					# ABC- #
												ABCO += 1
								elif px4 != '-':							## xxxx ##
									if px1 == px2:
										if px1 == px3:
											if px1 == px4:						# AAAA #
												AAAA += 1
											elif px1 != px4:					# AAAB #
												AAAB += 1
										elif px1 != px3:
											if px1 == px4:						# AABA #
												AABA += 1
											elif px1 != px4:
												if px3 == px4:					# AABB #
													AABB += 1
												elif px3 != px4:				# AABC #
													AABC += 1
									elif px1 != px2:
										if px1 == px3:
											if px1 == px4:						# ABAA #
												ABAA += 1
											elif px1 != px4:
												if px2 == px4:					# ABAB #
													ABAB += 1
												elif px2 != px4:				# ABAC #
													ABAC += 1
										elif px1 != px3:
											if px1 == px4:
												if px2 == px3:					# ABBA #
													ABBA += 1
												if px2 != px3:					# ABCA #
													ABCA += 1
											elif px1 != px4:
												if px2 == px3:
													if px2 == px4:				# ABBB #
														ABBB += 1
													elif px2 != px4:			# ABBC #
														ABBC += 1
												elif px2 != px3:
													if px2 == px4:				# ABCB #
														ABCB += 1
													elif px2 != px4:
														if px3 == px4:			# ABCC #
															ABCC += 1
														elif px3 != px4:		# ABCD #
															ABCD += 1

				with open(data_output,'a') as fout:
					fout.write(x1_label+" "+x2_label+" "+y1_label+" "+y2_label+" "+str(OOOO)+" "+str(OOOA)+" "+str(OOAO)+" "+str(OOAA)+" "+str(OOAB)+" "+str(OAOO)+" "+str(OAOA)+" "+str(OAOB)+" "+str(OAAO)+" "+str(OABO)+" "+str(OAAA)+" "+str(OAAB)+" "+str(OABA)+" "+str(OABB)+" "+str(OABC)+" "+str(AOOO)+" "+str(AOOA)+" "+str(AOOB)+" "+str(AOAO)+" "+str(AOBO)+" "+str(AOAA)+" "+str(AOAB)+" "+str(AOBA)+" "+str(AOBB)+" "+str(AOBC)+" "+str(AAOO)+" "+str(ABOO)+" "+str(AAOA)+" "+str(AAOB)+" "+str(ABOA)+" "+str(ABOB)+" "+str(ABOC)+" "+str(AAAO)+" "+str(AABO)+" "+str(ABAO)+" "+str(ABBO)+" "+str(ABCO)+" "+str(AAAA)+" "+str(AAAB)+" "+str(AABA)+" "+str(AABB)+" "+str(AABC)+" "+str(ABAA)+" "+str(ABAB)+" "+str(ABAC)+" "+str(ABBA)+" "+str(ABCA)+" "+str(ABBB)+" "+str(ABBC)+" "+str(ABCB)+" "+str(ABCC)+" "+str(ABCD)+"\n")
