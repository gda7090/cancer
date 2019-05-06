#!/usr/bin/python
import os,sys
import argparse
from rpy2.robjects.packages import importr
import rpy2.robjects as robj
r = robj.r

def hyper(k, m, n, N):
	if n < k:
		return 0
	else:
		return r['signif'](1 - r.phyper(k-1, n, N-n, m)[0])
	#if is_tail:
	#	return r['signif'](1 - r.phyper(k-1, n, N-n, m)[0])
	#else:
#		return r['signif'](r.phyper(k-1, n, N-n, m, is_tail = FALSE)[0])]
def chisq_test(k, m, n, N, **kargs):
	# k	m-k
	# n-k	N-n-m+k
	cmd = 'chisq.test(matrix(c(%d, %d, %d, %d), nc=2))' %(k, n-k, m-k, N-n-m+k)
	res = r(cmd)
	if res[7][0] >0:
		return r['signif'](res[2][0])
	else:
		return r['signif'](1 - res[2][0] / 2)

def improved_chisq_test(k, m, n, N, **kargs):
	# k	m-k
	# n-k	N-n-m+k
	max_e = 10
	s1r,s2r,s1c,s2c,t = m, (N-m), k, (m-k), float(N)
	es = (s1r*s1c/t, s1r*s2c/t, s2r*s1c/t, s2r*s2c/t)
	for e in es:
		if e < max_e:
			return r['signif'](fisher_test(k, m, n, N, **kargs))
def fisher_test(k, m, n, N, **kargs):
	alternative = kargs.get('alternative','greater')
	cmd = 'fisher.test(matrix(c(%d, %d, %d, %d), nc=2), alternative ="%s")' %(k, n-k, m-k, N-n-m+k, alternative)
	return r['signif'](r(cmd)[0][0])

def p_adjust(ps, method ='BH', test_num = None):
	if test_num:
		correct_ps = list(r['p.adjust'](ps, method = method, n = test_num))
		#correct_ps = list(r['p.adjust'](robj.FloatVector(ps), method = method, n = test_num))
	else:
		correct_ps = list(r['p.adjust'](ps, method = method))
		#correct_ps = list(r['p.adjust'](robj.FloatVector(ps), method = method))
	return correct_ps

def symbol():
	genes = {'primary':{},'second':{}}
	for line in open('AnnotationDB/Homo_sapiens.gene_info'):
		array = line.strip().split('\t')
		if line.startswith('#'):
			continue
		genes['primary'][array[2]] = array[2]
		if not array[4] == '-':
			for each in array[4].split('|'):
				genes['second'][each] = array[2]
	return genes

def maf2snv(maf,outdir):
	mut={}
	interest=[]
	for line in open(maf):
		array=line.strip().split('\t')
		if line.startswith('Hugo_Symbol') or array[0] == '.' or array[8]== 'Silent':
			continue
		if array[15] not in mut:
			mut[array[15]]=[]
		for gene in array[0].split(','):
			mut[array[15]].append(array[0])
			interest.append(gene)
	out=open(os.path.join(outdir,'all_mutgene.xls'),'w')
	gene=open(os.path.join(outdir,'all_go_gene'),'w')
	for sample in sorted(mut):
		out.write(sample+'\t'+'\t'.join(sorted(mut[sample]))+'\n')
		gene.write('\t'.join(sorted(mut[sample]))+'\n')
	out.close()
	return set(interest)
##the genelist file has two col , first is sample ,second is gene symbol
def genelist2snv(genelist,outdir) :
	mut={}
	interest=[]
	for line in open(genelist):
		array = line.strip().split("\t")
		if array[0] not in mut :
			mut[array[0]]=[]
		mut[array[0]].append(array[1])
		interest.append(array[1])
	out=open(os.path.join(outdir,'all_mutgene.xls'),'w')
	gene=open(os.path.join(outdir,'all_go_gene'),'w')
	for sample in sorted(mut):
		out.write(sample+'\t'+'\t'.join(sorted(mut[sample]))+'\n')
		gene.write('\t'.join(sorted(mut[sample]))+'\n')
	out.close()
	return set(interest)


parser = argparse.ArgumentParser(description="Pathway enrichment analysis,input file is maf or genelist")
parser.add_argument('--db',help="Pathway database",required=True)
parser.add_argument('--maf',help="the maf file",required=False)
parser.add_argument('--outdir',help="Output dir",required=True)
parser.add_argument('--prefix',help="Prefix for files",required=True)
parser.add_argument('--cut',help="The cutoff of genes for pathway enrich",default=1,type=int)
parser.add_argument('--genelist',help="The input is genelist",required=False)
parser.add_argument('--method',help="Method for enrichment",default='hyper',choices=['hyper','fisher','chisq'])
parser.add_argument('--fdr',help="Qvalue cutoff for enrichment",default=0.05,type=float)
argv = parser.parse_args()


if not os.path.exists(argv.outdir):
	assert not os.system('mkdir '+argv.outdir)

methods = {'hyper':'q_hyper', 'fisher':'q_fisher', 'chisq':'q_chisq'}
genes = symbol()
N_genes = set() ## bg
pathways = {}
for line in open(argv.db):
	if line.startswith('#'):
		continue
	array = line.strip().split('\t')
	if array[2] not in pathways:
		pathways[array[2]] = {'genes':set(),'msig':array[1],'name':array[4],'link':array[5],'dbtype':array[3]}
	gid = genes['primary'].get(array[0],genes['second'].get(array[0],array[0]))
	if not gid in pathways[array[2]]['genes']:
		pathways[array[2]]['genes'].add(gid)
	if gid not in N_genes:
		N_genes.add(gid)

#interest = set([each.split()[0] for each in open(argv.genelist).read().strip().split('\n')])
if argv.genelist :
	interest = genelist2snv(argv.genelist,argv.outdir)
else:
	interest=maf2snv(argv.maf,argv.outdir)
m_genes = N_genes & interest
N = len(N_genes)
m = len(m_genes)

pvalue = {'hyper':[],'chisq':[],'fisher':[]}
pathnames = []
k_n = []
p_genes = []
for each in pathways:
	k_genes = pathways[each]['genes'] & m_genes
	n = len(pathways[each]['genes'])
	k = len(k_genes)
	if k < argv.cut:
		continue
	pathnames.append(each)
	k_n.append([k,n])
	p_genes.append(k_genes)
	if argv.method == 'hyper':
		pvalue['hyper'].append(hyper(k, m, n, N))
	if argv.method == 'chisq':
		pvalue['chisq'].append(chisq_test(k, m, n, N))
	if argv.method == 'fisher':
		pvalue['fisher'].append(fisher_test(k, m, n, N))

qvalue = {}
if argv.method == 'hyper':
	qvalue['q_hyper'] = p_adjust(pvalue['hyper'])
if argv.method == 'fisher':
	qvalue['q_fisher'] = p_adjust(pvalue['fisher'])
if argv.method == 'chisq':
	qvalue['q_chisq']= p_adjust(pvalue['chisq'])

res = []
for i,each in enumerate(pathnames):
	tmp = { 'name':each,argv.method:pvalue[argv.method][i],methods[argv.method]:qvalue[methods[argv.method]][i],
		'k_n':k_n[i], 'genes':p_genes[i]}
	res.append(tmp)

res.sort(lambda x,y : cmp(x[argv.method],y[argv.method]))
out = open(os.path.join(argv.outdir,argv.prefix+'.pathway_enrichment.detailed.xls'),'w')
enrich = open(os.path.join(argv.outdir,argv.prefix+'.pathway_enrichment.xls'),'w')
out.write('Pathway\tPathway_name\tPathway_id\tDatabase\tfg_items\tbg_items\trich_factor\tpvalue\tqvalue\tGenes\tHttp_link\n')
enrich.write('Pathway\tPathway_name\tPathway_id\tDatabase\tfg_items\tbg_items\trich_factor\tpvalue\tqvalue\tGenes\tHttp_link\n')
for each in res:
	pid = each['name']
	out.write('%s\t%s\t%s\t%s\t%d\t%d\t%.2f\t%0.4e\t%0.4e\t%s\t%s\n' % \
	    (pid, pathways[pid]['name'], pathways[pid]['msig'], pathways[pid]['dbtype'], \
		each['k_n'][0], each['k_n'][1], float(each['k_n'][0])/float(each['k_n'][1]),\
	    each[argv.method][0],each[methods[argv.method]],\
		'|'.join(each['genes']), pathways[pid]['link']))
	if each[methods[argv.method]] < 0.05:
		enrich.write('%s\t%s\t%s\t%s\t%d\t%d\t%.2f\t%0.4e\t%0.4e\t%s\t%s\n' % \
		    (pid, pathways[pid]['name'], pathways[pid]['msig'], pathways[pid]['dbtype'], \
			each['k_n'][0], each['k_n'][1],float(each['k_n'][0])/float(each['k_n'][1]), \
		    each[argv.method][0],each[methods[argv.method]],\
		    '|'.join(each['genes']), pathways[pid]['link']))
out.close()
enrich.close()
