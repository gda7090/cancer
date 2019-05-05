#!/usr/bin/python
import os,sys
import argparse
import numpy as np
import copy
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

parser = argparse.ArgumentParser(description='Landscape.')
parser.add_argument('-m', '--maf', help='Maf file', required=True)
parser.add_argument('-s', '--samplelist', help='Samplelist', default=None)
parser.add_argument('-g', '--gene-pvalue', help='Genes\' pvalue', default=None)
parser.add_argument('-i', '--ipvalue', help='Pvalue column, 0 based, default 1', default=1,type=int)
parser.add_argument('-o', '--outdir', help='Output dir', required=True)
parser.add_argument('-c', '--cnv', help='gene cnv', default=None)
parser.add_argument('-t', '--threshold', help='CNV threshold for plot', type=float, default=3)
parser.add_argument('-p', '--prefix', help='Output file basename', default="Sample")
parser.add_argument('-n', '--gene-number', help='Gene number shown on the figure, default 30', default=30,type=int)
parser.add_argument('-w', '--width', help='Figure width(px), default', default=None,type=float)
parser.add_argument('-e', '--height', help='Figure height(px), default', default=None,type=float)
parser.add_argument('-k', '--scale', help='Figure scale', default=None,type=float)
parser.add_argument('-v', '--coverage', help='Sample coverage file', default=None)
args = parser.parse_args()

if not os.path.exists(args.outdir):
	assert not os.system('mkdir '+args.outdir)

type = {
	'Frame_Shift_Ins'   : ('Shift_indel', 7),
	'Frame_Shift_Del'   : ('Shift_indel', 7),
	'Nonsense_Mutation' : ('Nonsense', 6),
	'Nonstop_Mutation'  : ('Nonstop', 5),
	'In_Frame_Del'      : ('Inframe_indel', 4),
	'In_Frame_Ins'      : ('Inframe_indel', 4),
	'Missense_Mutation' : ('Missense', 3),
	'Splice_Site'       : ('Splice', 2),
	'gain'              : ('Gain', 1),
	'loss'              : ('Loss', 0.5),
	'Silent'            : ('Silent', 0)
}

#colmap = { 7 : ('Shift_indel', (1.0, 0.6274509803921569, 0.47843137254901963, 1.0)),
colmap = { 7 : ('Shift_indel', (0.28627451,0.03137255,0.93333333, 1.0)),
	   6 : ('Nonsense', (0.2549019607843137, 0.4117647058823529, 0.8823529411764706, 1.0)),
	   5 : ('Nonstop', (0.5568627450980392, 0.5568627450980392, 0.2196078431372549, 1.0)),
	   4 : ('Inframe_indel', (0.0, 0.9803921568627451, 0.6039215686274509, 1.0)),
	   3 : ('Missense', (0.9411764705882353, 0.3686274509803922, 0.30980392156862746, 1.0)),
	   2 : ('Splice', (0.6, 0.19607843137254902, 0.8, 1.0)),
	   1 : ('gain', (1.0, 0.8274509803921568, 0.6078431372549019, 1.0)),
	 0.5 : ('loss', (0.49411764705882355, 0.7529411764705882, 0.9333333333333333, 1.0)),
	   0 : ('', (0.9215686274509803, 0.9215686274509803, 0.9215686274509803, 1.0))
}

######### CNV
cnvs = {}
if args.cnv:
	i = 1
	cnv_samples = []
	for line in open(args.cnv):
		array = line.strip().split('\t')
		if i == 1:
			i += 1
			cnv_samples = array[1:]
			continue
		cnvs[array[0]] = {}
		for i,each in enumerate(array[1:]):
			cnvs[array[0]][cnv_samples[i]] = int(each)
#			ratio = 0
#			if int(each) == 0:
#				ratio = -2
#			else:
#				ratio = np.log2(int(each)) - 1
#			cnvs[array[0]][cnv_samples[i]] = ratio
cnv_threshold = args.threshold

matrix = {}
genes_maf = {}
samples_maf = []

mut_counts = {}
i = 1
for line in open(args.maf):
	if i == 1:
		i += 1
		continue
	array = line.strip().split('\t')
	if array[0] not in genes_maf:
		genes_maf[array[0]] = []
		matrix[array[0]] = {}
	if array[15] not in matrix[array[0]]:
		matrix[array[0]][array[15]] = []
	if array[15] not in mut_counts:
		mut_counts[array[15]] = np.array([0,0,0],dtype=np.float) ## NS, S, INDEL
	mutType = type.get(array[8],('RNA',0))

	matrix[array[0]][array[15]].append(mutType)
	if array[15] not in samples_maf:
		samples_maf.append(array[15])
	if array[15] not in genes_maf[array[0]] and mutType[1] > 0:
		genes_maf[array[0]].append(array[15])
	if array[8] in type:
		if array[9] == 'SNP' or array[9] == 'SNV':
			if array[8] == 'Silent':
				mut_counts[array[15]][1] += 1
			else:
				mut_counts[array[15]][0] += 1
		else:
			mut_counts[array[15]][2] += 1

mut_coverage = np.array([mut_counts[each] for each in samples_maf])
if args.coverage:
	coverage = dict([(each.split()[0],int(each.split()[1])) for each in open(args.coverage).read().strip().split('\n')[1:]])
	mut_coverage = np.array([mut_counts[each]/coverage[each]*1000000 for each in samples_maf])
## user defined samples and genes
samples_use = []
genes_use = []

if args.samplelist:
	if ',' in args.samplelist:
		samples_use = args.samplelist.split(',')
	else:
		samples_use = [each.split()[0] for each in open(args.samplelist).read().strip().split('\n')]

N = args.gene_number
genes_maf_count = dict([(each,len(genes_maf[each])) for each in genes_maf])
genes_top = [each[0] for each in sorted(genes_maf_count.items(),key = lambda d:d[1],reverse=True)]

pvalues = []
print args.gene_pvalue
if args.gene_pvalue:
	if ',' in args.gene_pvalue:
		genes_use = args.gene_pvalue.split(',')
	else:
		gp_pair = [(each.split()[0],float(each.split()[args.ipvalue])) for each in open(args.gene_pvalue).read().strip().split('\n') if not each.startswith('#')]
		gp_pair = sorted(gp_pair,key = lambda d:d[1])
		genes_use = [each[0] for each in gp_pair]
		pvalues = [each[1] for each in gp_pair][:N]
#		genes_use = [each.split()[0] for each in open(args.gene_pvalue).read().strip().split('\n')]
#		pvalues = [float(each.split()[1]) for each in open(args.gene_pvalue).read().strip().split('\n') if len(each.split()) > 1]


genes = genes_use or genes_top
genes = genes[:N]
samples = samples_use or samples_maf

assert len(pvalues) == 0 or len(pvalues) == len(genes)

matrix_raw = open(os.path.join(args.outdir,args.prefix+'.raw.matrix.txt'),'w')
matrix_cln = open(os.path.join(args.outdir,args.prefix+'.clean.matrix.txt'),'w')
png_file = os.path.join(args.outdir,args.prefix+'.landscape.png')
pdf_file = os.path.join(args.outdir,args.prefix+'.landscape.pdf')
matrix_raw.write('Gene\t'+'\t'.join(samples)+'\n')
matrix_cln.write('Gene\t'+'\t'.join(samples)+'\n')
data = []
order_genes = []
for gene in genes:
	if gene not in matrix:
		continue
	line1 = gene+'\t'
	line2 = gene+'\t'
	line = []
	
	for each in samples:
		type_tmp = [x[1] for x in matrix[gene].get(each,[('RNA',0)])]
		line1 += ','.join([str(x) for x in type_tmp])+'\t'
		line2 += str(max(type_tmp))+'\t'
		line.append(max(type_tmp))

	if args.cnv:
		cnv_tmp = [cnvs[gene][each] for each in samples]
		for i,each in enumerate(samples):
			if (not cnv_tmp[i] == 0) and line[i] == 0:
				if cnv_tmp[i] == 2:
					continue
				line[i] = 10+cnv_tmp[i]

	data.append(line)
	order_genes.append(gene)
	matrix_raw.write(line1.strip()+'\n')
	matrix_cln.write(line2.strip()+'\n')
matrix_raw.close()
matrix_cln.close()

order_samples = np.array(samples)
order_genes = np.array(order_genes)
data = np.array(data)
e = np.array((data > 0).tolist())
cols = range(data.shape[1])
rows = range(data.shape[0])
npcols = np.arange(data.shape[1])
nprows = np.arange(data.shape[0])
ncol = data.shape[1]
nrow = data.shape[0]
## mutually exclusive sort samples
for x in rows[::-1]:
	index = copy.deepcopy(cols)
	index.sort(key=lambda i : -e[x,i])
	e = e[::,index]
	data = data[::,index]
	order_samples = order_samples[index]
	mut_coverage = mut_coverage[index]

def format(x):
	if x % 1 == 0:
		return str(int(x))
	else:
		return str(x)

snv_matrix = open(os.path.join(args.outdir,args.prefix+'.snv.matrix.txt'),'w')
snv_matrix.write('Gene\t'+'\t'.join(order_samples)+'\n')
for x in rows:
	snv_matrix.write(order_genes[x]+'\t'+'\t'.join([format(d) for d in data[x,::]])+'\n')
snv_matrix.close()

def T2O(two):
	one = []
	for each in two:
		one += list(each)
	return one

event_pos = [(npcols+0.5).tolist() for each in rows]
event_offset = [np.array([each+0.5]*ncol).tolist() for each in rows[::-1]]
event_color = []
for row in rows:
	event_color.append([colmap[data[row,col]][1] for col in cols])

event_pos = [[each] for each in T2O(event_pos)]
event_offset = T2O(event_offset)
event_color = T2O(event_color)

fw=max(ncol*0.23,3)
fh = 5
scale = args.scale or 1
if nrow < 10:
	fh = 2.5
	scale = 0.5
width = args.width or fw
height = args.height or fh

fig = plt.figure(figsize=(width,height))
f1 = fig.add_subplot(111)
plt.sca(f1)
heatmap = f1.eventplot(event_pos,orientation='horizontal', lineoffsets=event_offset, linelengths=0.96, linewidths=5, linestyles='solid',colors=event_color)
print data.shape[1]
print data.shape[0]
f1.set_xbound(lower = 0, upper = data.shape[1])
f1.set_ybound(lower = 0, upper = data.shape[0])
f1.set_yticks(np.arange(data.shape[0])+0.5, minor=False)
f1.set_xticks(np.arange(data.shape[1])+0.5, minor=False)
f1.set_yticklabels(order_genes[::-1], minor=False,size=3)
#f1.xaxis.tick_top()
f1.set_xticklabels(order_samples, minor=False,size=2.2,rotation=90)

for t in f1.yaxis.get_major_ticks():
	t.tick1On = False
	t.tick2On = False
for t in f1.xaxis.get_ticklines():
	t.set_visible(False)

f1.set_aspect('auto',adjustable='box-forced',anchor='C')
divider = make_axes_locatable(f1)
### mutation rate
ns_cov = [each[0] for each in mut_coverage]
s_cov = [each[1] for each in mut_coverage]
indel_cov = [each[2] for each in mut_coverage]

ns_color = (0.06274509803921569, 0.3058823529411765, 0.5450980392156862, 1.0)
s_color = (0.4, 0.803921568627451, 0.6666666666666666, 1.0)
indel_color = (240.0/255,128.0/255,128.0/255,1.0)
top_bar = divider.append_axes("top",0.8*scale,pad=0.1*scale,sharex=f1)
plt.setp(top_bar.get_xticklabels(),visible=False)
top_bar.bar(cols,ns_cov,width=0.88,color=ns_color,label='non-silent',linewidth=0,edgecolor='none',figure=mpl.figure.Figure(frameon=False))
top_bar.bar(cols,s_cov,width=0.88,color=s_color,bottom=ns_cov,label='silent',linewidth=0,edgecolor='none',figure=mpl.figure.Figure(frameon=False))
top_bar.bar(cols,indel_cov,width=0.88,color=indel_color,bottom=[ns_cov[x]+s_cov[x] for x in cols],label='indel',linewidth=0,edgecolor='none',figure=mpl.figure.Figure(frameon=False))
top_bar.set_xbound(lower = 0, upper = ncol)
top_bar.spines['top'].set_visible(False)
top_bar.spines['bottom'].set_visible(False)
top_bar.spines['right'].set_visible(False)
top_bar.tick_params(axis='y', direction='out')
if pvalues:
	top_bar.legend(bbox_to_anchor=(1.0,0.6),loc=2,framealpha=0,fontsize=3,ncol=1)
else:
	top_bar.legend(bbox_to_anchor=(0.5,1.02),loc=0,framealpha=0,fontsize=3,ncol=3)

if args.coverage:
	top_bar.set_ylabel('Mutations/Mb',fontsize=3.5)
else:
	top_bar.set_ylabel('Mutations',fontsize=3.5)
#top_bar.text(-6,15,'a',fontsize=12,weight='bold')
for t in top_bar.xaxis.get_major_ticks():
    t.tick1On = False
    t.tick2On = False
for t in top_bar.get_xticklabels():
    t.set_visible(False)
for t in top_bar.yaxis.get_major_ticks():
    t.tick1line.set_markersize(2)
    t.tick2line.set_markersize(0)
    t.label1.set_fontsize(3)

## gene counts
gene_count = [float(np.array(each>0,dtype=np.int).sum())/ncol*100 for each in data]
xmax = int(max(gene_count)/10+1)*10
left_bar = divider.append_axes("left", 0.8*scale, pad=0.5*scale, sharey=f1)
plt.setp(left_bar.get_yticklabels(),visible=False)
left_bar.barh(rows[::-1],gene_count,height=0.88, color=(0.0, 0.3607843137254902, 0.6, 1.0),linewidth=0)
left_bar.spines['left'].set_visible(False)
left_bar.spines['top'].set_visible(False)
left_bar.spines['right'].set_visible(False)
left_bar.set_ybound(lower = 0, upper = nrow)
left_bar.set_xlim([xmax, 0])
xticks = range(0,xmax+1,10)
if len(xticks)>7:
	xticks = range(0,xmax+1,20)
left_bar.set_xticks(xticks,minor=False)
left_bar.set_xticklabels(range(0,xmax+1,10),size='xx-small',minor=False)
left_bar.grid(True, which='major',axis='x', color='#EDEDED', linestyle='-', linewidth=1.0)
left_bar.set_axisbelow(True)
left_bar.set_xlabel('Percent of mutations (%)',fontsize=3)
for t in left_bar.yaxis.get_major_ticks():
	t.tick1On = False
	t.tick2On = False
for t in left_bar.xaxis.get_major_ticks():
	t.tick1line.set_markersize(2)
	t.tick2line.set_markersize(0)
	t.label1.set_fontsize(3)
for t in left_bar.get_yticklabels():
	t.set_visible(False)


## pvalues
def logpv(pv):
	import math
	pv2 = np.array(pv,dtype=np.float)
	pv2[pv2==0] = min(pv2[pv2>0])/1000
	return [-math.log(v,10) for v in pv2]

if pvalues:
	logPvalue = logpv(pvalues)
	xmax = int(max(logPvalue)+1)
	right_bar = divider.append_axes("right", 0.6*scale, pad=0.1*scale, sharey=f1)
	plt.setp(right_bar.get_yticklabels(),visible=False)
	right_bar.barh(rows[::-1],logPvalue,height=0.88,color='#B7B7B7',linewidth=0)
	right_bar.set_ybound(lower = 0, upper = nrow)
	right_bar.spines['left'].set_visible(False)
	right_bar.spines['top'].set_visible(False)
	right_bar.spines['right'].set_visible(False)
	right_bar.set_xlabel('-log'+'$_{10}$'+'Pvalue',fontsize=3)
	right_bar.set_xlim([0,xmax])
	right_bar.set_xticks(range(0,xmax+1,max(1,xmax/5)),minor=False)
	right_bar.grid(True, which='major',axis='x', color='#EDEDED', linestyle='-', linewidth=1.0)
	right_bar.set_axisbelow(True)
	right_bar.set_xticklabels(range(0,xmax+1,max(xmax/5,1)),size='xx-small',minor=False)
	for t in right_bar.yaxis.get_major_ticks():
		t.tick1On = False
		t.tick2On = False
	for t in right_bar.xaxis.get_major_ticks():
		t.tick1line.set_markersize(2)
		t.tick2line.set_markersize(0)
		t.label1.set_fontsize(3)
	for t in right_bar.get_yticklabels():
		t.set_visible(False)

## legend
legends = [colmap[each] for each in set(T2O(data)) if not each == 0]
bottom = divider.append_axes("bottom", 0.3*scale, pad=0.6*scale, sharex=f1)
bottom._frameon = False
plt.setp(bottom.get_xticklabels()+bottom.get_yticklabels(),visible=False)
colorstart = 2
line = int(0.3*scale/(float(height)/nrow))

for leg_lab,leg_col in legends:
	bottom.scatter(colorstart,line,s=20,c=leg_col,marker='s')
	bottom.text(colorstart+2,line,leg_lab,fontsize=5,verticalalignment='center')
	print colorstart
	print line
	colorstart += 10
	if colorstart+10 >ncol:
		line -= 1
		colorstart = 2
bottom.set_xbound(lower = 0, upper = ncol)
for t in bottom.yaxis.get_major_ticks():
    t.tick1On = False
    t.tick2On = False
for t in bottom.xaxis.get_ticklines():
    t.set_visible(False)
for t in bottom.yaxis.get_ticklines():
    t.set_visible(False)

plt.show()
fig.tight_layout()
fig.savefig(png_file,dpi=500)
fig.savefig(pdf_file)

open(os.path.join(os.path.join(args.outdir,'figure.parameters')),'w').write('--width:\t%f\n--height:\t%f\n--scale:\t%f\n'%(width,height,scale))

