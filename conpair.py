import os
import sys
import optparse
import threading
import Queue
import time

gatk_jar="gatk/GenomeAnalysisTK-3.1-1-g07a4bf8/GenomeAnalysisTK.jar"
conpair_dir="onpair/Conpair-master"
desc = """Conpair is software to evaluate concordance and contamination estimator for matched tumor-normal pairs. """
parser = optparse.OptionParser(version='%prog version 1.0 8/March/2017', description=desc)
parser.add_option('-i', '--input', help='INPUT: samplelist', action='store')
parser.add_option('-o', '--outdir', help='OUTDIR: the output directory')
parser.add_option('-J', '--java', help='PATH to JAVA [java by default]', default='java', action='store')
parser.add_option('-G', '--gatk', help='GATK JAR [$GATK by default]',default=gatk_jar, action='store')
parser.add_option('-D', '--conpair_dir', help='CONPAIR DIR [$CONPAIR_DIR by default]',default=conpair_dir, action='store')
parser.add_option('-g', '--genome',help='genome: B37,hg19,B38,hg38; default:B37',default='B37',choices=["B37","hg19","B38","hg38"],action='store')
parser.add_option('-m', '--xmx_java', help='Xmx java memory setting [default: 12g]', default='12g', action='store')
parser.add_option('-t', '--thread', type=int, help='Thread number for calculate, default 4', default=4)
parser.add_option('-C', '--min_cov', help='MIN COVERAGE TO CALL GENOTYPE [default: 10]', default=10, type='int', action='store')
parser.add_option('-Q', '--min_mapping_quality', help='MIN MAPPING QUALITY TO CALL GENOTYPE [default: 10]', default=10, type='int', action='store')
parser.add_option('-B', '--min_base_quality', help='MIN BASE QUALITY TO CALL GENOTYPE [default: 20]', default=20, type='int', action='store')
parser.add_option('-H', '--normal_homozygous_markers_only', help='USE ONLY MARKERS THAT ARE HOMOZYGOUS IN THE NORMAL SAMPLE (concordance will not be affected by CNV)', default=False, action='store_true')
parser.add_option('-r', '--rpath', help='R pathway [default:R-3.3.1/bin]', default='R-3.3.1/bin', action='store')
parser.add_option('-T', '--threshold', help='The threshold value of contamination value[default: 0.05]', default='0.05',action='store')
(opts, args) = parser.parse_args()

#thread
exitFlag = 0
class myThread ( threading.Thread ):
	def __init__ (self, threadID, name, q):
		threading.Thread.__init__ (self)
		self.threadID = threadID
		self.name = name
		self.q = q
	def run (self):
		print 'Starting ' + self.name
		doJob(self.name,self.q)
		print 'Exiting: ' + self.name
def doJob(threadName, q):
	while not exitFlag:
		queueLock.acquire()
		if not workQueue.empty():
			data = q.get()
			queueLock.release()
			assert not os.system(data['script'])
			print '%s process %s' % (threadName, data['name'])
		else:
			queueLock.release()
		time.sleep(1)


i=0
samplelist={}
# key: patient; value:list of samples
patient_sample={}
for line in open(opts.input):
	line=line.lstrip('#').strip()
	array=line.split('\t')
	if i==0:
		sample_index='None'
		patient_index='None'
		bam_index='None'
		type_index='None'
		try:
			sample_index=array.index('Sample')
			patient_index=array.index('Patient')
			bam_index=array.index('Bam_path')
			type_index=array.index('Type')
		except ValueError:
			pass
		if sample_index=='None'	or patient_index=='None' or bam_index=='None':
			print 'The input file must contain: Sample,Patient,Bam_Path.'
			exit(0)
		i += 1
		continue
	#dictionary: patient_sample  key:patient  value:sample
	if array[patient_index] not in patient_sample.keys():
		patient_sample[array[patient_index]]=[]
		patient_sample[array[patient_index]].append(array[sample_index])
	else:
		patient_sample[array[patient_index]].append(array[sample_index])
	#dictionary: samplelist
	samplelist[array[sample_index]]={}
	samplelist[array[sample_index]]['patient']=array[patient_index]
	samplelist[array[sample_index]]['bam_path']=array[bam_index]
	if type_index != 'None':
		samplelist[array[sample_index]]['type']=array[type_index]
		
sample_order=[]
for item in patient_sample.keys():
	sample_order += patient_sample[item]

genome_files = {
	'B37':{
		'reference_path':'Database/Genome/human_B37/GRCh37.fasta',
		'markers_file':os.path.join(opts.conpair_dir,'data', 'markers', 'GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt'),
		'markers_bed':os.path.join(opts.conpair_dir,'data', 'markers', 'GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.bed')
	},
	'hg19':{
		'reference_path':'Database/Genome/human_hg19/ucsc.hg19.fasta',
		'markers_file':os.path.join(opts.conpair_dir,'data', 'markers', 'hg19.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt'),
		'markers_bed':os.path.join(opts.conpair_dir,'data', 'markers', 'hg19.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.bed')
	},
	'B38':{
		'reference_path':'Database/Genome/human_B38/human_B38.fa',
		'markers_file':'Conpair/Conpair-master/data/markers/GRCh38.changed.txt',
		'markers_bed':'Conpair/Conpair-master/data/markers/GRCh38.changed.bed'
	},
	'hg38':{
		'reference_path':'Genome/human_hg38/hg38.fasta',
		'markers_file':'Conpair/Conpair-master/data/markers/hg38.changed.txt',
		'markers_bed':'Conpair/Conpair-master/data/markers/hg38.changed.bed'
	}
}
################################
#### calculate pileup
#################################
pileup_dir=os.path.join(opts.outdir,'pileup')
if not os.path.exists(pileup_dir):
	assert not os.system('mkdir '+pileup_dir)

threadList = ['Thread-%d' % each for each in range(0,opts.thread)]
myjobs = []

for each in sample_order:
	script="python %s -D %s -R %s -M %s -G %s -J %s -m %s -B %s -O %s" % (os.path.join(opts.conpair_dir,'scripts','run_gatk_pileup_for_sample.py'),opts.conpair_dir,genome_files[opts.genome]['reference_path'],genome_files[opts.genome]['markers_bed'],opts.gatk,opts.java,opts.xmx_java,samplelist[each]['bam_path'],os.path.join(pileup_dir,each+'.pileup'))
	myjobs.append({'name':each,'script':script})
	open(os.path.join(pileup_dir,'calc-pileup.%s.sh'%each),'w').write(script+'\n')

queueLock = threading.Lock()
workQueue = Queue.Queue(0)
threads = []
threadID = 1
## start threads
for tName in threadList:
	thread = myThread(threadID, tName, workQueue)
	thread.start()
	threads.append(thread)
	threadID += 1
## Queue and job
queueLock.acquire()
for eachjob in myjobs:
	workQueue.put(eachjob)
queueLock.release()

while not workQueue.empty():
	pass
exitFlag = 1

for t in threads:
	t.join()

##################################
##### calculate concordance
##################################
concordance_dir=os.path.join(opts.outdir,'concordance')
if not os.path.exists(concordance_dir):
        assert not os.system('mkdir '+concordance_dir)

concordance={}
for sample1 in sample_order:
	concordance[sample1]={}
	for sample2 in sample_order:
		if sample2==sample1:
			concordance[sample1][sample2]=1.0
		else:
			concor_script = "export CONPAIR_DIR=%s\n" % (conpair_dir)
			concor_script += "export PYTHONPATH=${PYTHONPATH}:%s\n" % (os.path.join(conpair_dir,'modules'))
			concor_script += "python %s -C %s -Q %s -B %s -H %s -N %s -T %s -O %s -M %s" % (os.path.join(opts.conpair_dir,'scripts','verify_concordance.py'),opts.min_cov,opts.min_mapping_quality,opts.min_base_quality,opts.normal_homozygous_markers_only,os.path.join(pileup_dir,sample1+'.pileup'),os.path.join(pileup_dir,sample2+'.pileup'),os.path.join(concordance_dir,'concordance.%s-N.%s-T.txt'%(sample1,sample2)),genome_files[opts.genome]['markers_file'])
			open(os.path.join(concordance_dir,'calc-concordance.%s-N.%s-T.sh'%(sample1,sample2)),'w').write(concor_script+'\n')
			assert not os.system(concor_script)
			concor_value=float(open(os.path.join(concordance_dir,'concordance.%s-N.%s-T.txt'%(sample1,sample2))).readline().replace('Concordance:','').replace('%','').strip())/100
			concordance[sample1][sample2]=concor_value
# generate xls formate
out=open(os.path.join(opts.outdir,'concordance.pairwise.xls'),'w')
out.write('\t'+'\t'.join(sample_order)+'\n')
for sample1 in sample_order:
	line=sample1+'\t'
	for sample2 in sample_order:
		line += '%f\t'% concordance[sample1][sample2]
	line=line.rstrip()
	out.write(line+'\n')
out.close()
# plot
sample_num=len(sample_order)
rscript = '''
library(ggplot2)
setwd("%s")
datas1 <- read.table("%s",header=TRUE,row.names=1,sep="\\t")
datas2 <- data.frame()
for(i in 1:dim(datas1)[1]){
	for(j in 1:dim(datas1)[2]){
		datas3 <- data.frame(X1=c(i),X2=c(j),value=c(datas1[i,j]))
		datas2=rbind(datas2,datas3)
	}
}
sample_num <- %d
base_size <- sample_num*3.2
p <- ggplot(datas2,aes(X2,X1)) + geom_tile(aes(fill = value),colour="red",width=rep(1,dim(datas1)[2]))
p <- p + scale_fill_gradientn(colours=c("white","yellow","red","black"))
p <- p + theme_grey(base_size = base_size) + labs(x="",y="")
p <- p + scale_x_continuous(expand=c(0,0),labels=names(datas1),breaks=1:dim(datas1)[2])
p <- p + scale_y_continuous(expand=c(0,0),labels=row.names(datas1),breaks=1:dim(datas1)[1])
p <- p + theme(axis.ticks=element_blank(),axis.text.x=element_text(size=base_size*1.2,angle=270,hjust=0,colour="gray0"),axis.text.y=element_text(size=base_size*1.2,hjust=1,colour="gray0"))
#legend
p <- p+xlab("") + ylab("") + labs(fill="rate")
p <- p + theme(legend.title=element_text(size=18))
p <- p + theme(legend.key.size = unit(sample_num/4,'cm'),legend.key.width=unit(1,'cm'))
p <- p + theme(legend.text=element_text(colour="black",angle=0,size=15,hjust=3,vjust=3))
ggsave(filename="%s", plot=p,width=sample_num*1.5,height=sample_num*1.5)
'''%(opts.outdir,os.path.join(opts.outdir,'concordance.pairwise.xls'),sample_num,os.path.join(opts.outdir,'concordance.pairwise.heatmap.pdf'))
open(os.path.join(opts.outdir,'concordance.plot.R'),'w').write(rscript)
assert not os.system('%s %s' % (os.path.join(opts.rpath,'Rscript'),os.path.join(opts.outdir,'concordance.plot.R')))

##################################
##### calculate contamination
##################################
if type_index != 'None':
	sample_pair=[]
	for item in patient_sample.keys():
		if len(patient_sample[item])<=1:
			continue
		for each in patient_sample[item]:
			if samplelist[each]['type']=='N':
				for each2 in patient_sample[item]:
					if samplelist[each2]['type']=='T':
						sample_pair.append([each,each2])
	#print sample_pair
	contamination_dir=os.path.join(opts.outdir,'contamination')
	if not os.path.exists(contamination_dir):
		assert not os.system('mkdir '+contamination_dir)
	open(os.path.join(opts.outdir,'contamination.detail.txt'),'w').write('-'*40+'\n')
	open(os.path.join(opts.outdir,'contamination.xls'),'w').write('Sample\tValue\tType\n')
	for each in sample_pair:
		contamination_script= "export export CONPAIR_DIR=%s\n" % (conpair_dir)
		contamination_script += "export PYTHONPATH=${PYTHONPATH}:%s\n" % (os.path.join(conpair_dir,'modules'))
		contamination_script += "python %s -D %s -M %s -N %s -T %s -O %s" % (os.path.join(opts.conpair_dir,'scripts','estimate_tumor_normal_contamination.py'),opts.conpair_dir,genome_files[opts.genome]['markers_file'],os.path.join(pileup_dir,each[0]+'.pileup'),os.path.join(pileup_dir,each[1]+'.pileup'),os.path.join(contamination_dir,'contamination.'+each[0]+'-N.'+each[1]+'-T.txt'))
		open(os.path.join(contamination_dir,'calc-contamination.%s-N.%s-T.sh'%(each[0],each[1])),'w').write(contamination_script+'\n')
		assert not os.system(contamination_script)
		open(os.path.join(opts.outdir,'contamination.detail.txt'),'a').write('Normal:\t%s\nTumor:\t%s\n'%(each[0],each[1]))
		os.system('cat %s >> %s'%(os.path.join(contamination_dir,'contamination.'+each[0]+'-N.'+each[1]+'-T.txt'),os.path.join(opts.outdir,'contamination.detail.txt')))
		open(os.path.join(opts.outdir,'contamination.detail.txt'),'a').write('-'*40+'\n')
		line_num_con=0
		lines_list=open(os.path.join(contamination_dir,'contamination.'+each[0]+'-N.'+each[1]+'-T.txt'),'r').readlines()
		#print lines_list
		value1=float(lines_list[0].strip().split()[-1].strip().replace('%',''))/100
		open(os.path.join(opts.outdir,'contamination.xls'),'a').write('%s\t%f\tNormal\n'%('%s-N.%s-T'%(each[0],each[1]),value1))
		value2=float(lines_list[1].strip().split()[-1].strip().replace('%',''))/100
		open(os.path.join(opts.outdir,'contamination.xls'),'a').write('%s\t%f\tTumor\n'%('%s-N.%s-T'%(each[0],each[1]),value2))
	rscript2 = '''
library(ggplot2)
Args <- commandArgs()
filename=Args[6]
mydatax <- read.table(filename,header=TRUE,sep="\\t")
mydatax$Value <- -log10(mydatax$Value)

base_size <- 20
p <- ggplot(mydatax,aes(x=Sample,y=Value,fill=Type))+geom_bar(position=position_dodge(1.0),width=0.9,stat="identity")
p <- p+scale_x_discrete(limits=mydatax$Sample)
p <- p+scale_y_continuous(expand=c(0,0),limits=c(0,4))
p <- p+theme_gray(base_size = base_size*3)
p <- p+theme(axis.text.x=element_text(size=base_size*2,angle=0,colour="black"))
p <- p+theme(axis.text.y=element_text(size=27))
p <- p + labs(x = "", y= expression(Contamination~~group("(",italic(-log[10]),")")))
p <- p + theme(axis.title.y= element_text(size=33))
p <- p + theme_bw()
#p <- p + theme(panel.grid.major.y=element_line(size=1,linetype =3,color="black"))
#legend
p <- p + theme(legend.title=element_text(size=15))
p <- p + theme(legend.key.size = unit(1.2,"cm"),legend.key.width=unit(1.2,'cm'))
p <- p + theme(legend.text=element_text(colour="black",angle=0,size=15,hjust=3,vjust=3))

th_ch=as.numeric(Args[8])
th=as.numeric(Args[8])
th=-log10(th)
# adding rect line text
#p <- p + geom_hline(aes(yintercept=seq(0,4,by=1)),colour="#990000",linetype="dashed")
p <- p + geom_hline(yintercept=seq(0, 4, by=1),colour="black",linetype="dashed")
p <- p + geom_hline(yintercept=th,colour='red')
p <- p + geom_rect(ymin=0,ymax=th,xmin=0,xmax=nrow(mydatax)+1,fill="#FFED97",alpha=0.1)
p <- p + geom_text(label=th_ch,y=th+0.25,x=nrow(mydatax)+0.22,hjust=0.5, vjust=-0.5,size=4 ,color="red")
p <- p + geom_segment(x=nrow(mydatax)+0.2,y=th+0.25,xend=nrow(mydatax),yend=th,arrow = arrow(length = unit(0.1,"cm")),color="red")
ggsave(filename=Args[7], plot=p,width=nrow(mydatax)*8/6,height=nrow(mydatax))
'''
	open(os.path.join(opts.outdir,'contamination.plot.R'),'w').write(rscript2)
	print '%s %s %s %s %s' % (os.path.join(opts.rpath,'Rscript'),os.path.join(opts.outdir,'contamination.plot.R'),os.path.join(opts.outdir,'contamination.xls'),os.path.join(opts.outdir,'contamination.normal-tumor.bar.pdf'),opts.threshold)
	assert not os.system('%s %s %s %s %s' % (os.path.join(opts.rpath,'Rscript'),os.path.join(opts.outdir,'contamination.plot.R'),os.path.join(opts.outdir,'contamination.xls'),os.path.join(opts.outdir,'contamination.normal-tumor.barplot.pdf'),opts.threshold))
	#assert not os.system('rm -f %s'%(os.path.join(opts.outdir,'contamination.xls')))
