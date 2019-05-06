#!Python-2.7.11/bin/python
import os,sys
import argparse
import numpy as np
import copy
sys.path.append('Python-2.7.11/lib/python2.7/site-packages')
from ete2 import Tree,TreeStyle,TextFace,NodeStyle

parser = argparse.ArgumentParser(description='Phylogenetic Tree analysis for Cancer Evolution.')
parser.add_argument('-i', '--maf', help='Maf file', required=True)
parser.add_argument('-g', '--germline', help='Germline mutation file', default=None)
parser.add_argument('-n', '--root-number', help='Root number for construction tree, default None', type=int, default=None)
parser.add_argument('-s', '--samplelist', help='Samplelist', default=None)
parser.add_argument('-o', '--outdir', help='Output dir', required=True)
parser.add_argument('-p', '--patient', help='Patient ID', default="patient")
parser.add_argument('-f', '--func', help='Function filter, default: coding',default='coding',choices=['coding','non-synonymous','exon','gene','all'])
parser.add_argument('-c', '--chr', help='Chr column index for Germline, default=0', default=0,type=int)
parser.add_argument('-t', '--position', help='Position column index for Germline, default=1', default=1,type=int)
parser.add_argument('-r', '--ref', help='Ref-base column index for Germline, default=3', default=3,type=int)
parser.add_argument('-a', '--alt', help='Alt-base column index for Germline, default=4', default=4,type=int)
parser.add_argument('-x', '--gene', help='Gene column index for Germline, default=6', default=6,type=int)
parser.add_argument('-l', '--indel', help='Wether INDEL used for tree construction', action='store_true')
parser.add_argument('-d', '--drivergene', help='Driver Genes for marke, default None', default=None)
args = parser.parse_args()


def filter_funcs(func, type='coding'):
	coding_funcs  = ['Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Silent', 'Splice_Site', \
		'Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins']
	nonsyn_funcs  = ['Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', \
		'Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins']
	exon_funcs    = ['3\'UTR','5\'UTR','Missense_Mutation','Nonsense_Mutation','Nonstop_Mutation','Silent','Splicing_Site',\
		'Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins']
	notgene_funcs = ['3\'Flank','5\'Flank','IGR']
	if type == 'all':
		return True
	if type == 'coding':
		return func in coding_funcs
	if type == 'non-synonymous':
		return func in nonsyn_funcs
	if type == 'exon':
		return func in exon_funcs
	if type == 'gene':
		return func not in notgene_funcs
	return True

if not os.path.exists(args.outdir):
	assert not os.system('mkdir %s'%args.outdir)
## user defined samples and genes
samples_use = []
if args.samplelist:
	if ',' in args.samplelist:
		samples_use = args.samplelist.split(',')
	else:
		samples_use = [each.split()[0] for each in open(args.samplelist).read().strip().split('\n')]


matrix = {}
mutation = {}
ids_maf = []
samples_maf = []

i = 1
for line in open(args.maf):
	if i == 1:
		i += 1
		continue
	array = line.strip().split('\t')
	if array[15] not in samples_use and args.samplelist:
		continue
	## INDEL filter
	if (not args.indel) and (array[9] == 'INS' or  array[9] == 'DEL' or  array[9] == 'INDEL'):
		continue
	## func filter
	if not filter_funcs(array[8], args.func):
		continue
	id = '%s:%s:%s>%s' % (array[4], array[5], array[10], array[12])
	if id not in ids_maf:
		ids_maf.append(id)
		matrix[id] = {}
		mutation[id] = {'ref':array[10], 'gene':array[0]}
	mutation[id][array[15]] = array[12]
	matrix[id][array[15]] = float(array[33])/(int(array[33])+int(array[32])) ## AF
	if array[15] not in samples_maf:
		samples_maf.append(array[15])

samples = samples_use or samples_maf
samples.insert(0,args.patient)

ids_germline = []
if args.germline:
	i = 1
	for line in open(args.germline):
		if line.startswith('##'):
			continue
		if i == 1:
			i += 1
			continue
		array = line.strip().split('\t')
		id = '%s:%s:%s>%s' % (array[args.chr], array[args.position], array[args.ref], array[args.alt])
		if id in mutation:
			continue
		ids_germline.append(id)
		mutation[id] = {'ref':array[args.ref], 'gene':array[args.gene]}
		mutation[id][args.patient] = array[args.alt]
		matrix[id] = {args.patient:1}
#		mutation[id][args.patient] = array[args.alt]
#		matrix[id] = {args.patient:1}
		for each in samples[1:]:
			mutation[id][each] = array[args.alt]
			matrix[id][each] = 0

prefix = os.path.join(args.outdir,args.patient+'.'+args.func)

matrix_cln = open(prefix+'.clean.matrix.txt','w')
matrix_cln.write('ID\t'+'\t'.join(samples)+'\n')
data = []

N = 0
if args.root_number:
	N = args.root_number

id_select = []
for id in ids_maf:
	line = [matrix[id].get(each,0) for each in samples]
	if len([each for each in line if each>0]) == 0:
		continue
	if not (args.germline and args.root_number):
		line[0] = int(len([each for each in line[1:] if each>0]) == len(line[1:]))
	data.append(line)
	matrix_cln.write('\t'.join([id]+[str(each) for each in line])+'\n')
	id_select.append(id)

for id in ids_germline:
	line = [matrix[id].get(each,0) for each in samples]
	data.append(line)
	matrix_cln.write('\t'.join([id]+[str(each) for each in line])+'\n')
	id_select.append(id)
	if args.root_number and N <= 0:
		break
	N = N -1
if args.root_number and N > 0:
	for n in range(N):
		id = 'G:%d'%n
		ref = ['A','T','C','G'][n%4]
		alt = ['A','T','C','G'][(n+1)%4]
		mutation[id] = {'ref':ref, 'gene':'NA'}
		mutation[id][args.patient] = alt
		matrix[id] = {args.patient:1}
		for each in samples[1:]:
			mutation[id][each] = ref
			matrix[id][each] = 0
		data.append([1]+[0]*(len(samples)-1))
		id_select.append(id)
		matrix_cln.write(id+'\t1\t'+'\t'.join(['0']*(len(samples)-1))+'\n')
matrix_cln.close()

## root
order_ids = np.array(id_select)
order_samples = np.array(samples)
d = np.array(data)
f = np.array((d>0).tolist(),dtype=np.int)

## add root
#root_stat = np.array([each.sum()==e.shape[1] for each in e], dtype=np.int)
#f = e.transpose().tolist()
#f.insert(0,root_stat)
#f = np.array(f).transpose()
cols = range(f.shape[1])
rows = range(f.shape[0])

## mutual exclude sort
for x in cols[::-1]:
	index = copy.deepcopy(rows)
	index.sort(key = lambda i : -f[i,x])
	f = f[index, ::]
	order_ids = order_ids[index]

stat_file = open(prefix+'.stat','w')
meg_file = open(prefix+'.nucl.meg','w')
meg_file.write('#mega\n!Title evolution;\n!Format DataType=DNA indel=-\n\n')
stat_file.write('\t%d\t%d\n' % (f.shape[1], f.shape[0]))
sample_mapping = {}

for i in np.arange(f.shape[1]):
	column = [str(each) for each in f[::,i].tolist()]
	sample_mapping['sample%d'%i] = order_samples[i]
	sample = 'sample%d                '%i
	stat_file.write('%s%s\n' %(sample[0:10], ''.join(column)))
	meg_file.write('#'+order_samples[i]+'\n')
	seq = ''.join([mutation[each].get(order_samples[i],mutation[each]['ref']) for each in order_ids])
	meg_file.write('\n'.join([seq[each*60:min((each+1)*60,len(seq))] for each in range(len(seq)/60+1)])+'\n')
#	meg_file.write(''.join([mutation[each].get(order_samples[i],mutation[each]['ref']) for each in order_ids])+'\n')
stat_file.close()
meg_file.close()


## phylip 
open(os.path.join(args.outdir,'pars.params.'+args.func),'w').write('%s\nY\n'%(prefix+'.stat'))
phylip = '''
log=%s/log
cd %s && \\
cp Evolution/phylip-3.695/src/font1 fontfile && \\
if [ -f outtree ];then rm -f outtree; fi && \\
if [ -f outfile ];then rm -f outfile; fi && \\
if [ -f plotfile ];then rm -f plotfile; fi && \\
echo %s ... pars infer tree ... > $log
echo ........................................... >> $log
Evolution/phylip-3.695/exe/pars < pars.params.%s >> $log && \\
mv -f outfile %s.pars_out && \\
mv -f outtree %s.tree
#echo %s ... drawtree ... > $log
#echo ........................................... >> $log
#Evolution/phylip-3.695/exe/drawtree < drawtree.params.%s >> $log && \\
#mv -f plotfile %s.tree.pdf && \\
#convert %s.tree.pdf %s.tree.png
''' % (args.outdir, args.outdir, args.patient, args.func, prefix, prefix, args.patient, args.func, prefix, prefix, prefix)
open(os.path.join(args.outdir,args.func+'.createTree.sh'),'w').write(phylip)
assert not os.system('sh %s' % (os.path.join(args.outdir,args.func+'.createTree.sh')))

t = Tree(prefix+'.tree')
for node in t:
	node.name = sample_mapping.get(node.name,node.name)
t.write(format=1,outfile=prefix+'.tree.nk')

open(os.path.join(args.outdir,'drawtree.params.'+args.func),'w').write('%s\nY\n'%(prefix+'.tree.nk'))
drawtree = '''
cd %s && \
drawtree < drawtree.params.%s && \
mv -f plotfile %s.phylip.tree.pdf && \
convert %s.phylip.tree.pdf %s.phylip.tree.png
''' % (args.outdir, args.func, prefix, prefix, prefix)
open(os.path.join(args.outdir,args.func+'.drawTree.sh'),'w').write(drawtree)
assert not os.system('sh %s' % (os.path.join(args.outdir,args.func+'.drawTree.sh')))


t.name = 'node0'
t.set_outgroup(args.patient)

i = 1
for node in t.traverse("levelorder"):
	if node.name == "":
		node.name = "node%d"%i
		i += 1 

node2labels = t.get_cached_content(store_attr="name")
drivergenes = []
if args.drivergene:
	drivergenes = [each.split()[0] for each in open(args.drivergene).read().strip().split('\n')]
branches = ['Germline','Germline-node1']
pat2branch = {'1'+'0'*(len(order_samples)-1):'Germline', '1'*len(order_samples):'Germline-node1'}
edge_muts = {'Germline':{'id':[],'driver':[]}, 'Germline-node1':{'id':[],'driver':[]}}

for node in t.iter_descendants("preorder"):
	if node.name == 'node0' or node.up.name == 'node0':
		continue
	pattern = ''.join([str(int(each in node2labels[t&node.name])) for each in order_samples])
	branches.append('%s-%s'%(node.up.name,node.name))
	pat2branch[pattern] = '%s-%s'%(node.up.name,node.name)
	edge_muts['%s-%s'%(node.up.name,node.name)] = {'id':[],'driver':[]}

branches.append('Out-of-tree')
pat2branch['0'*len(order_samples)] = 'Out-of-tree'
edge_muts['Out-of-tree'] = {'id':[],'driver':[]}

for i,id in enumerate(order_ids):
	pattern = ''.join(np.array(f[i],dtype=np.str))
	if pattern in pat2branch:
		edge_muts[pat2branch[pattern]]['id'].append(id)
		if mutation[id]['gene'] in drivergenes:
			edge_muts[pat2branch[pattern]]['driver'].append(id)
	else:
		edge_muts['Out-of-tree']['id'].append(id)
		if mutation[id]['gene'] in drivergenes:
			edge_muts['Out-of-tree']['driver'].append(id)


mutation_in_tree = open(prefix+'.tree_mutations.xls','w')
mutation_in_tree.write('Branches\t'+'\t'.join(['%s.Mutation\t%s.Gene'%(each,each) for each in branches])+'\n')
mutation_in_tree.write('Count\t'+'\t'.join(['%d\t%d' %(len(edge_muts[each]['id']),\
	len(set([mutation[x]['gene'] for x in edge_muts[each]['id']]))) for each in branches])+'\n')

driver_in_tree = open(prefix+'.tree_mutations.driver.xls','w')
driver_in_tree.write('Branches\t'+'\t'.join(['%s.Mutation\t%s.Gene'%(each,each) for each in branches])+'\n')
driver_in_tree.write('Count\t'+'\t'.join(['%d\t%d' %(len(edge_muts[each]['driver']),\
	len(set([mutation[x]['gene'] for x in edge_muts[each]['driver']]))) for each in branches])+'\n')

i = 0
while True:
	lines = ['']
	if i == 0:
		lines[0] = 'Mutations'
	flag = False
	for each in branches:
		if len(edge_muts[each]['id']) > i:
			lines.append(edge_muts[each]['id'][i])
			lines.append(mutation[edge_muts[each]['id'][i]]['gene'])
			flag = True
		else:
			lines.append('')
			lines.append('')
	i += 1
	if not flag:
		break
	mutation_in_tree.write('\t'.join(lines)+'\n')
mutation_in_tree.close()
i = 0
while True:
	lines = ['']
	if i == 0:
		lines[0] = 'Mutations'
	flag = False
	for each in branches:
		if len(edge_muts[each]['driver']) > i:
			lines.append(edge_muts[each]['driver'][i])
			lines.append(mutation[edge_muts[each]['driver'][i]]['gene'])
			flag = True
		else:
			lines.append('')
			lines.append('')
	i += 1
	if not flag:
		break
	driver_in_tree.write('\t'.join(lines)+'\n')
driver_in_tree.close()


ts = TreeStyle()
ts.show_leaf_name = True
ts.margin_bottom = 40
ts.margin_top = 40
ts.margin_left = 10
ts.margin_right = 10
ts.show_scale = False
ts.scale=0.3

style1 = NodeStyle()
style1['size'] = 4
style1['fgcolor'] = '#5500ff'
style1['hz_line_color'] = '#55aa00'
style1['vt_line_color'] = '#55aa00'
style1['hz_line_width'] = 2
style1['vt_line_width'] = 2
style2 = NodeStyle()
style2['size'] = 6
style2['fgcolor'] = '#0055ff'
style2['hz_line_color'] = '#55aa00'
style2['vt_line_color'] = '#55aa00'
style2['hz_line_width'] = 2
style2['vt_line_width'] = 2
t0 = t&'node1'
t0.add_face(TextFace('1'),column=0,position='branch-right')
t0.img_style = style1
for node in t0.iter_descendants("preorder"):
	if node.name.startswith('node'):
		(t0&node.name).add_face(TextFace(node.name.replace('node','')),column=0,position='branch-right')
		node.img_style = style1
	else:
		node.img_style = style2

t0.render(prefix+'.ete.tree.png',dpi=300,tree_style=ts,w=2000)
#t0.render(prefix+'.ete.tree.pdf',tree_style=ts,w=2000)
'''
t0 = t&'node1'
t0.render(prefix+'.ete.tree.png',dpi=300,w=1000,h=600)
t0.render(prefix+'.ete.tree.pdf',w=1000,h=600)
'''

Rscript='''
setwd('%s')
dat<-read.table("%s",head=T,row.names=1)
dat[dat>0] <- 1
dat$%s<-as.numeric(sapply(1:nrow(dat),function(x) all(dat[x,]>0)))
d<-data.matrix(dat)
h<-hclust(dist(t(d)))
png("%s",type="cairo-png")
plot(h,hang=-1,xlab="",ylab="",axes=FALSE,sub="")
dev.off()
'''%(args.outdir, prefix+'.clean.matrix.txt', args.patient, args.patient+'.tree.png')
#open(os.path.join(args.outdir,'hclust.R'),'w').write(Rscript)
#assert not os.system('Rscript %s' % (os.path.join(args.outdir,'hclust.R')))
