#!/usr/bin/env python

import os
import sys

def cr_sqrt(x):
	# returns the square root of a list of variables [x]
	if x == 0.0:
		return 0.0
	else:
		M = 1.0
		xN = x 
		while xN >= 2.0:
			xN = 0.25*xN
			M = 2.0*M
		while xN < 0.5:
			xN = 4.0*xN
			M = 0.5*M
		A = xN
		B = 1.0-xN
		while 1==1:
			A = A*(1.0+0.5*B)
			B = 0.25*(3.0+B)*B*B
			if B < 1.0E-15:
				return A*M


def cr_sum(x):
	# returns the sum of a list by removing 'na' values (-9)
	x2 = [ i for i in x if i!=-9 ]
	if len(x2) != 0:
		return(sum(x2))
	else:
		return(-9)


def cr_mean(x):
	# returns the mean of a list
	# removes the 'na' values from the list (-9)
	x2 = [ i for i in x if i!=-9 ] # -9 = 'na' value
	nElement = len(x2)
	if nElement == 0:
		return(-9)
	else:
		return(sum(x2)/(1.0 * nElement))

def cr_std(x, exp_X):
	# returns the standard variation of a list
	# removes the 'na' values from the list (-9)
	x2 = [ i for i in x if i!=-9 ] # -9 = 'na' value
	nElement = len(x2)
	if nElement == 0:
		return(-9)
	else:
		if sum(x2) == 0:
			return(0)
		else:
			A = sum( [ i**2 for i in x2 ] )
			A = A/(1.0 * nElement)
			return(cr_sqrt(A-exp_X**2))


def cr_pearsonR(x, y):
	# computes correlation between arrays x and y
	## removes 'na' values (-9)
	x2 = [ x[i] for i in range(len(x)) if x[i]!=-9 and y[i]!=-9 ]
	y2 = [ y[i] for i in range(len(y)) if x[i]!=-9 and y[i]!=-9 ]
	
	sumXi = 0.0
	sumYi = 0.0
	sumSquareX = 0.0
	sumSquareY = 0.0
	sumXiYi = 0.0
	
	nX = len(x2)

	if nX!=0:
		for i in range(nX):
			sumXi += x2[i];
			sumYi += y2[i];
			sumSquareX += x2[i]*x2[i];
			sumSquareY += y2[i]*y2[i];
			sumXiYi += x2[i]*y2[i];
		
		numerator = nX*sumXiYi - sumXi * sumYi
		denom1 = cr_sqrt(nX*sumSquareX - sumXi*sumXi)
		denom2 = cr_sqrt(nX*sumSquareY - sumYi*sumYi)
		if denom1 == 0 or denom2 == 0:
			return(0)
		else:
			return(numerator/(denom1*denom2))
	else:
		return(-9)


def compFreq(sequences, segsites):
	# returns derived allele frequency for a number of 'segsites' positions
	nDerAll = []
	nInd = len(sequences)
	nPair = nInd*(nInd-1)/2.0
	for i in range(segsites):
		nDerAll.append(0)
		for j in sequences:
			if j[i] == "1":
				nDerAll[i] += 1.0
	pi = [ i * (nInd-i) / nPair for i in nDerAll ]
	freq = [ i/(1.0 * nInd) for i in nDerAll ]
	res = {}
	res['nDer'] = nDerAll # list of the number of occurence of the derived allele at each position
	res['pi_SNPs'] = pi # list of pi measured at each position
	res['pi'] = sum(pi) # total pi = sum of the pi = total number of differences over all positions + over pairwise comparisons between species
	res['freq'] = freq # frequencies of the derived allele at each position
	res['nSNPs'] = len([ i for i in freq if i>0.0 and i<1.0 ]) # number of SNPs (with frequencies in ]0, 1[) in the alignement
	return(res)


def sites(freqA, freqB, segsites):
	# test whether the SNP is a polymorphism specific to species A (sxA), species B (sxB), is found in both species (Ss) or differentialy fixed among species (Sf)
	sxA, sxB, ss, sf, same = 0, 0, 0, 0, 0
	successive_ss = [0]
	previous_site = ""
	for i in range(segsites):
		if freqA[i] == 0:
			if freqB[i] == 0:
				same += 1
			if freqB[i] == 1:
				sf += 1
				if previous_site == "ss" or previous_site == "":
					successive_ss.append(0)
				previous_site = "sf"
			if freqB[i] < 1:
				sxB += 1
				if previous_site == "ss" or previous_site == "":
					successive_ss.append(0)
				previous_site = "sxB"
			continue
		if freqA[i] == 1:
			if freqB[i] == 1:
				same += 1
			if freqB[i] == 0:
				sf += 1
				if previous_site == "ss" or previous_site == "":
					successive_ss.append(0)
				previous_site = "sf"
			if freqB[i] > 0:
				sxB += 1
				if previous_site == "ss" or previous_site == "":
					successive_ss.append(0)
				previous_site = "sxB"
			continue
		else:
			if freqB[i] == 0 or freqB[i] == 1:
				sxA += 1
				if previous_site == "ss" or previous_site == "":
					successive_ss.append(0)
				previous_site = "sxA"
			else:
				ss += 1
				if previous_site == "ss" or previous_site == "":
					successive_ss[len(successive_ss)-1] += 1
				previous_site = "ss"
			continue
	res = {'sxA':sxA, 'sxB':sxB, 'sf':sf, 'ss':ss, 'same':same, 'successive_ss':max(successive_ss)}
	return(res)

def piTot(nDerA, nDerB, nSamA, nSamB, segsites):
	# returns pi tot for the pooled populations from two vectors of allele count per locus
	piT = []
	nTot = nSamA + nSamB
	for i in range(segsites):
		if nDerA==nSamA and nDerB==nSamB:
			piT.append(-9)
		if nDerA==0 and nDerB==0:
			piT.append(-9)
		else:
			tmp = nDerA[i] + nDerB[i]
			tmp = (nTot - tmp) * tmp / (nTot*(nTot-1)/2.0) # "n ancesral" x "n derived" / C(2,k)
			piT.append(tmp)
	return(piT)


def Fst(piA, piB, piT):
	# returns Fst as: 1 - mean_pi_over_populations / pi_in_the_whole_alignment
	if piT==0:
		res = -9
	else:
		res = 1.0 - (piA+piB)/(2.0*piT)
	return(res)



# read information about loci from the bpfile 
if os.path.isfile("bpfile") == False:
	sys.exit("\n\t\033[1;31;40mERROR: bpfile was not found\n\033[0m\n")

infile = open("bpfile", "r")


tmp = infile.readline() # first empty line
L = [ float(i) for i in infile.readline().strip().split("\t") ]
nSamA = [ int(i) for i in infile.readline().strip().split("\t") ]
nSamB = [ int(i) for i in infile.readline().strip().split("\t") ]

nLoci = len(L) 
infile.close()

a1_spA, a1_spB, a2_spA, a2_spB= [], [], [], []
for nsam in nSamA:
	a1_spA.append(sum([ 1.0/i for i in range(1, nsam) ]))
	a2_spA.append(sum([ 1.0/(i**2) for i in range(1, nsam) ]))
for nsam in nSamB:
	a1_spB.append(sum([ 1.0/i for i in range(1, nsam) ]))
	a2_spB.append(sum([ 1.0/(i**2) for i in range(1, nsam) ]))

## ms' output file
outfile = open("ABCstat.txt", "w")
# header
res = "dataset\tbialsites_avg\tbialsites_std\t"
res += "piA_avg\tpiA_std\t"
res += "piB_avg\tpiB_std\t"
res += "FST_AB_avg\tFST_AB_std\n"

outfile.write(res)

test = 0 
nSim_cnt = 0 # count the number of treated multilocus simulations
nLoci_cnt = 0 # count the number of treated loci within a simulation


# READ THE MS's OUTPUTFILE
##for line in infile:
for line in sys.stdin: # read the ms's output from the stdin
	line = line.strip()
	if "segsites" in line:
		if nLoci_cnt == 0:
			#ss_sf, noSs_sf, ss_noSf, noSs_noSf = 0, 0, 0, 0
			bialsites = []
			sxA, sxB = [], []

			piA, piB = [], []
			FST_AB = []
		
		nLoci_cnt += 1
		nSam_cnt = 0 # count the number of treated individuals within a locus
		test = 1
		segsites = int(line.split(":")[1])
		bialsites.append(segsites)
		spA, spB, spC, spD = [], [], [], []
		continue
	if test == 1:
		if segsites == 0:
			test = 0
			sxA.append(0)
			sxB.append(0)

			piA.append(0)
			piB.append(0)
			FST_AB.append(-9)
			
			#noSs_noSf += 1
		if segsites != 0:
			if "positions" not in line and line!="\n":
				nSam_cnt += 1
				if nSam_cnt <= nSamA[nLoci_cnt - 1]:
					spA.append(line.strip())
				if nSam_cnt > nSamA[nLoci_cnt - 1] and nSam_cnt <= (nSamA[nLoci_cnt - 1] + nSamB[nLoci_cnt - 1]):
					spB.append(line.strip())
				
				# end of the block of sequences -> start the calculations of various statistics
				if nSam_cnt == (nSamA[nLoci_cnt - 1] + nSamB[nLoci_cnt - 1] ):
					tmpA = compFreq(spA, segsites)
					freqA = tmpA['freq']
					piA.append(tmpA['pi']/L[nLoci_cnt - 1])

					tmpB = compFreq(spB, segsites)
					freqB = tmpB['freq']
					piB.append(tmpB['pi']/L[nLoci_cnt - 1])

					## Sx
					tmpABC = compFreq(spA+spB, segsites)
					freqABC = tmpABC['freq']
					
					# SxA
					tmp_A_B = sites(freqA, freqB, segsites)
					sxA.append(tmp_A_B['sxA']/(1.0*L[nLoci_cnt - 1]))
					# SxB
					tmp_BA = sites(freqB, freqABC, segsites)
					sxB.append(tmp_BA['sxA']/(1.0*L[nLoci_cnt - 1]))

					# Fst
					# vector of piT over segsites
					piT_AB = piTot(tmpA['nDer'], tmpB['nDer'], nSamA[nLoci_cnt - 1], nSamB[nLoci_cnt - 1], segsites)
					
					# mean Fst
					FST_AB.append(Fst(cr_sum(tmpA['pi_SNPs']), cr_sum(tmpB['pi_SNPs']), cr_sum(piT_AB)))
					
				
	# compute average and std over of statistics over loci
	if nLoci_cnt != 0 and len(sxA) == nLoci:
		test = 0
		nSim_cnt += 1
		nLoci_cnt = 0
		
		# statistics
		bialsites_avg = cr_mean(bialsites)
		bialsites_std = cr_std(bialsites, bialsites_avg)
		
		#sxA_avg = cr_mean(sxA)
		#sxA_std = cr_std(sxA, sxA_avg)

		piA_avg = cr_mean(piA)
		piA_std = cr_std(piA, piA_avg)
		piB_avg = cr_mean(piB)
		piB_std = cr_std(piB, piB_avg)
		FST_AB_avg = cr_mean(FST_AB)
		FST_AB_std = cr_std(FST_AB, FST_AB_avg)
		
		res = ""
		res += "{0}\t{1:.5f}\t{2:.5f}\t".format(nSim_cnt-1, bialsites_avg, bialsites_std)
		res += "{0:.5f}\t{1:.5f}\t".format(piA_avg, piA_std)
		res += "{0:.5f}\t{1:.5f}\t".format(piB_avg, piB_std)
		res += "{0:.5f}\t{1:.5f}\t".format(FST_AB_avg, FST_AB_std)

		res += "\n"
		outfile.write(res)
infile.close()
outfile.close()

