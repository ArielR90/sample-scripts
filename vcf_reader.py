#!/usr/bin/python

import sys
import re
import os
import numpy as np

if len(sys.argv) < 2:
	print "#1- working directory\n\
#2- frequency histogram output\n\
#3- allele counts file from bedtools genomecov\n\
#4- frequency by SNP output\n\
#5- depth by SNP location output"
	sys.exit()

wkdir = sys.argv[1]
files = sorted(os.listdir(sys.argv[1]))

# file to write histogram
outf = open(sys.argv[2], 'w+')
header = ['Fish', 'Genotype', 'Tissue', 'Total_SNP', 'Under_0.5_SNP', 
'Under_0.1_SNP', 'Under_0.01_SNP', 'Under_0.005_SNP', 'Total_Indel', 
'Under_0.5_Indel', 'Under_0.1_Indel', 'Under_0.01_Indel', 'Under_0.005_Indel']
outf.write('\t'.join([str(i) for i in header]))
outf.write('\n')

sites = {}
n = 0

# counts file from genomecov
# generate sites dictionary
PA = np.genfromtxt(sys.argv[3], dtype = 'str')
for i in range(len(PA)):
	site = '-'.join(PA[i][:4])
	if site in sites:
		print "duplicate"
	else:
		sites[site] = [PA[i][4], [0.0 for j in range(18)], [0 for j in range(18)]]
		print PA[i]

# populate depth of coverage array for each site in sites
for file in files:
	if file.endswith('cov'):
		cov = np.genfromtxt(wkdir+'/'+file, dtype=None)
		for i in range(len(cov)):
			pos = cov[i][0]+'-'+str(cov[i][1])
			k = [ v for v in sites.keys() if pos in v]
			if k != []:
				sites[k[0]][2][n] = cov[i][2]
		n +=1		

# populate allele frequency array for each site in sites
# count allele frequencies in each  range for each sample
n = 0
for file in files:
	if file.endswith("PELE.vcf"):
		with open(wkdir+'/'+file, "r") as f:
			snps = [0,0,0,0,0]
			indels = [0,0,0,0,0]
			if file.startswith('W'):
				geno = 'WT'
			else:
				geno = 'fancD1'
			labels = file.split('_')
			for line in f:
				line = line.strip('\n')
				if line.startswith('#'):
					pass
				else:
					line=re.split('\t|;|=', line)
					depth = int(line[8])
					af = float(line[10])
					pos = '-'.join(line[:2])+'-'+'-'.join(line[3:5])
					if pos in sites:
						sites[pos][1][n] = float(af)
						sites[pos][2][n] = depth
					else:
						print "new alleles?", pos
					if 'INDEL' in line:
						indels[0] += 1
						if af < 0.5 and af > 0.1:
							indels[1] += 1
						elif af <= 0.1 and af > 0.01:
							indels[2] += 1
						elif af <= 0.01 and af > 0.005:
							indels[3] += 1
						elif af <= 0.005:
							indels[4] += 1
					else:
						snps[0] += 1
						if af < 0.5 and af > 0.1:
							snps[1] += 1
						elif af <= 0.1 and af > 0.01:
							snps[2] += 1
						elif af <= 0.01 and af > 0.005:
							snps[3] += 1
						elif af <= 0.005:
							snps[4] += 1
			out = [labels[0], geno, labels[1], snps[0], snps[1], snps[2], snps[3], snps[4], indels[0], indels[1], indels[2], indels[3], indels[4]]
			outf.write('\t'.join([str(i) for i in out]))
			outf.write('\n')
		n +=1

# write SNP frequency and depth files
# frequency by site output file
outf2 = open(sys.argv[4], 'w+')
# depth by site output
outf3 = open(sys.argv[5], 'w+')
header = ['Chromosome', 'Position', 'Reference', 'Alternate', 
'FD1_G', 'FD1_HK', 'FD1_T', 'FD2_G', 'FD2_HK', 'FD2_T', 'FD3_G', 'FD3_HK', 'FD3_T', 
'WT1_G', 'WT1_HK', 'WT1_T', 'WT2_G', 'WT2_HK', 'WT2_T', 'WT3_G', 'WT3_HK', 'WT3_T']
outf2.write('\t'.join([str(i) for i in header]))
outf2.write('\n')
outf3.write('\t'.join([str(i) for i in header]))
outf3.write('\n')
for site in sites:
	label = site.split('-')
	outf2.write('\t'.join(label)+'\t'+'\t'.join([str(i) for i in sites[site][1]])+'\n')
	outf3.write('\t'.join(label)+'\t'+'\t'.join([str(i) for i in sites[site][2]])+'\n')
