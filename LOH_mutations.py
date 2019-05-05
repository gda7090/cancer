#!/usr/bin/env python
import os 
import sys
import argparse,HTSeq
from argparse import RawTextHelpFormatter
import math

parser = argparse.ArgumentParser(description="LOH module.",formatter_class=RawTextHelpFormatter)
parser.add_argument('-i','--cnv',help="A tab seperated file with these infos:\nchrom, start, end, copynumber, genotype, genotype score(optional)",required=True)
parser.add_argument('-m','--mutation', help="Mutation file", required=True)
parser.add_argument('-t','--mutation-type', help="Type of mutation file, default: maf", choices=['vcf','maf'],default='maf')
parser.add_argument('-o','--out',help="Output result",required=True)
parser.add_argument('-c','--chr',type=int,default=0,help="For CNV, chr index number, default: 0")
parser.add_argument('-s','--start',type=int,default=1,help="For CNV, start index number, default: 1")
parser.add_argument('-e','--end',type=int,default=2,help="For CNV, end index number, default: 2")
parser.add_argument('-n','--cn',type=int,default=4,help="For CNV, copy number index number, default: 4")
parser.add_argument('-g','--genotype',type=int,default=7,help="For CNV, copy number\' genotype(AAB, AAA ...), default: 7")
parser.add_argument('-gs','--genotype-score',type=int,help="Score for genotype filter.", default=None)
parser.add_argument('-gh','--genotype-score-threshold',help="Expression for score filter, such as: lt:score, gt:score, eq:score", default=None)
parser.add_argument('-tsg','--tsgene', help="TSGenes database,TSGene.list", default='TSGene.list')
argv=parser.parse_args()


def gs_filter(score):
	if argv.genotype_score_threshold:
		e,s = argv.genotype_score_threshold.split(':',1)
		assert e in ['lt', 'gt', 'eq']
		s = float(s)
		flag = False
		if e == 'lt' and score < s and score != '-1':
			flag = True
		if e == 'gt' and score > s and score != '-1':
			flag = True
		if e == 'eq' and score == s and score != '-1':
			flag = True

		return flag
	else:
		return True

def func_filter(func):
	myfuncs = ['Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation', 'Splice_Site','In_Frame_Del', \
		'Frame_Shift_Del','In_Frame_Ins', 'Frame_Shift_Ins','Frame_Shift_Ins']
	if func in myfuncs:
		return True
	return False

TSGenes = [each.split()[0] for each in open(argv.tsgene).read().strip().split('\n')]


## read info from cnv file
cnvinfo = HTSeq.GenomicArrayOfSets( "auto", stranded=False )
i = 1
for line in open(argv.cnv):
	if i == 1:
		i += 1
		continue
	array = line.strip().split('\t')
	if array[9] =="germline":
		continue
	genotypes = set(array[argv.genotype])

	info_tmp = (array[argv.chr],array[argv.start],array[argv.end],array[argv.cn],array[argv.genotype])
	## genotype filter
	if not (int(array[argv.cn]) > 0 and len(genotypes) == 1):
		continue
	## score filter
	if argv.genotype_score:
		info_tmp = (array[argv.chr],array[argv.start],array[argv.end],array[argv.cn],array[argv.genotype],array[argv.genotype_score])
		if not gs_filter(float(array[argv.genotype_score])):
			continue
	iv = HTSeq.GenomicInterval(array[argv.chr],int(array[argv.start]),int(array[argv.end]))
	cnvinfo[iv] += info_tmp


if argv.mutation_type == 'maf':
	outfile = open(argv.out,'w')
	head = 'Chrom\tStart\tEnd\tCopynumber\tCopynumber Genotype\t'
	if argv.genotype_score:
		head += 'Genotype Score\t'
	i = 1
	for line in open(argv.mutation):
		if i == 1:
			i += 1
			head += line.strip()+'\tTSGene\n'
			outfile.write(head)
			continue
		array = line.strip().split('\t')
		if not func_filter(array[8]):
			continue
		iv = HTSeq.GenomicPosition(array[4],int(array[5]))
		tsg = ''
		if array[0] in TSGenes:
			tsg = 'TSGene'
		if cnvinfo[iv]:
			outfile.write('\t'.join(list(cnvinfo[iv])[0])+'\t'+line.strip()+'\t%s\n'%(tsg))
	outfile.close()
	
elif argv.mutation_type == 'vcf':
	outfile = open(argv.out,'w')
	head = 'Chrom\tStart\tEnd\tCopynumber\tCopynumber Genotype\t'
	if argv.genotype_score:
		head += 'Genotype Score\t'
	for line in open(argv.mutation):
		if line.startswith('##'):
			continue
		if line.startswith('#'):
			head += line[1:]
			outfile.write(head)
			continue
		array = line.strip().split('\t')
		iv = HTSeq.GenomicPosition(array[0],int(array[1]))
		if cnvinfo[iv]:
			outfile.write('\t'.join(list(cnvinfo[iv])[0])+'\t'+line)
	outfile.close()
