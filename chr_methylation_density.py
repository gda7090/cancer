#!/usr/bin/python
import os.path, os
import argparse
import HTSeq

def parse_option():
	parser = argparse.ArgumentParser(description = 'Pipeline for methylation statistics and plot')
	parser.add_argument('--CX_report', required = True, help = 'The genome wide cytosine report, by bismark')
	parser.add_argument('--methy', required = True, help = 'Binormal tested methylated C result, by mC_identification')
	parser.add_argument('--outdir', required = True, help = 'Output directory.')
	parser.add_argument('--fa', required = True, help= "genome file from UCSC")
	parser.add_argument('--label', required = True, help = "Labels for this sample")
	parser.add_argument('--cutoff', default = 5, type = int, help = 'The minimal depth for C called methylation')
	parser.add_argument('--window', default = 10000, type = int, help = 'The window size(10000bp) for chromosome density calculation.')
	parser.add_argument('--chrlist', required = True, help = "chr names list that are chosen to be analyzed, one column in this file, and the list should not contain too many chromosomes")
	argv = parser.parse_args()
	return argv

def gp2gi(gp, fs, fe, window): ##	fs = feature start position
	k = (gp.pos-fs)/window
	gi = HTSeq.GenomicInterval(gp.chrom, k*window+fs, (k+1)*window+fs, '.')
	return gi
def step_gis(start, end, window, chr):
	gis = [HTSeq.GenomicInterval(chr, i, i + window) for i in range(start, end, window)]
	return gis

def create_dir(dir):
	if not os.path.exists(dir):
		os.mkdir(dir)

def dens(c, mc):
	if c == 0:
		return float(0)
	else:
		return float(mc)/c*100


if __name__ == '__main__':
	argv = parse_option()
	outdir = argv.outdir
	label = argv.label
	cutoff = argv.cutoff

	create_dir(outdir)

	chr_lens = {}
	dict_chrs={}
	for line in open(argv.chrlist):
		if line.strip() == '':
			continue
		else:
			each=line.strip()
			dict_chrs[each]='chromosome'

	for seq_record in HTSeq.FastaReader(argv.fa):
		chr_lens[seq_record.name] = len(seq_record.seq)
	chrlist = []

	types = ['CG', 'CHG', 'CHH']

	Ccount_bin = {'C':{}, 'CG':{}, 'CHG':{}, 'CHH':{}}			# total C
	mCcount_bin = {'C':{}, 'CG':{}, 'CHG':{}, 'CHH':{}}			# mC

	for line in open(argv.CX_report):
		if line.strip() == '':
			continue
		each = line.strip().split('\t')
		if dict_chrs.has_key(each[0]):
			chr = each[0]
		else:
			continue 
		if chr not in chrlist:
			chrlist.append(chr)
		strand = each[2]
		each_iv = HTSeq.GenomicPosition(chr, int(each[1]), strand)
		type = each[5]
		depth = int(each[3]) + int(each[4])

		if not Ccount_bin[type].has_key(chr):
			Ccount_bin[type][chr] = {}

		##	calculate C count in each bin
		gi = gp2gi(each_iv, fs=0, fe=chr_lens[chr], window=argv.window)
		if not Ccount_bin[type][chr].has_key(gi):
			Ccount_bin[type][chr][gi] = {'+':0,'-':0}
		if not depth < argv.cutoff:
			Ccount_bin[type][chr][gi][strand] += 1

	for line in open(argv.methy):
		if line.strip() == '':
			continue
		each = line.strip().split('\t')
		depth = int(each[3]) + int(each[4])
		chr = each[0]
		type = each[5]
		strand = each[2]
		each_iv = HTSeq.GenomicPosition(chr, int(each[1]), strand)

		if not mCcount_bin[type].has_key(chr):
			mCcount_bin[type][chr] = {}

		##	calculate mC count in each bin
		gi = gp2gi(each_iv, fs=0, fe=chr_lens[chr], window=argv.window)
		if not mCcount_bin[type][chr].has_key(gi):
			mCcount_bin[type][chr][gi] = {'+':0, '-':0}
		if not depth < argv.cutoff:
			mCcount_bin[type][chr][gi][strand] += 1

	###	CHR density
	den_dir = outdir + '/chr_density'
	create_dir(den_dir)
	for chr in chrlist:
		if chr_lens[chr] < argv.window:
			continue
		chr_steps = step_gis(start = 0, end = chr_lens[chr], window = argv.window, chr = chr)
		CG_C, CHG_C, CHH_C, CG_mC, CHG_mC, CHH_mC = {},{},{},{},{},{}

		if Ccount_bin['CG'].has_key(chr):
			CG_C = Ccount_bin['CG'][chr]
		if Ccount_bin['CHG'].has_key(chr):
			CHG_C = Ccount_bin['CHG'][chr]
		if Ccount_bin['CHH'].has_key(chr):
			CHH_C = Ccount_bin['CHH'][chr]
		if mCcount_bin['CG'].has_key(chr):
			CG_mC = mCcount_bin['CG'][chr]
		if mCcount_bin['CHG'].has_key(chr):
			CHG_mC = mCcount_bin['CHG'][chr]
		if mCcount_bin['CHH'].has_key(chr):
			CHH_mC = mCcount_bin['CHH'][chr]

		chr_stat = open(den_dir + '/' + chr + '_mCX_density.stat', 'w')
		chr_stat.write('#bin\tplus\tminus\tCG\tCHG\tCHH\n')
		print "chr_stat have been finished\n"
		for i,eachgi in enumerate(chr_steps):
			CG_C_p, CG_C_m, CHG_C_p, CHG_C_m, CHH_C_p, CHH_C_m = 0, 0, 0, 0, 0, 0
			CG_mC_p, CG_mC_m, CHG_mC_p, CHG_mC_m, CHH_mC_p, CHH_mC_m = 0, 0, 0, 0, 0, 0
			if CG_C.has_key(eachgi):
				CG_C_p = CG_C[eachgi]['+']
				CG_C_m = CG_C[eachgi]['-']
			if CHG_C.has_key(eachgi):
				CHG_C_p = CHG_C[eachgi]['+']
				CHG_C_m = CHG_C[eachgi]['-']
			if CHH_C.has_key(eachgi):
				CHH_C_p = CHH_C[eachgi]['+']
				CHH_C_m = CHH_C[eachgi]['-']
			if CG_mC.has_key(eachgi):
				CG_mC_p = CG_mC[eachgi]['+']
				CG_mC_m = CG_mC[eachgi]['-']
			if CHG_mC.has_key(eachgi):
				CHG_mC_p = CHG_mC[eachgi]['+']
				CHG_mC_m = CHG_mC[eachgi]['-']
			if CHH_mC.has_key(eachgi):
				CHH_mC_p = CHH_mC[eachgi]['+']
				CHH_mC_m = CHH_mC[eachgi]['-']
			tmp = '\t'.join([str(each) for each in CG_C_p, CG_C_m, CHG_C_p, CHG_C_m, CHH_C_p, CHH_C_m, CG_mC_p, CG_mC_m, CHG_mC_p, CHG_mC_m, CHH_mC_p, CHH_mC_m])
			C_p, C_m, mC_p, mC_m = CG_C_p + CHG_C_p + CHH_C_p, CG_C_m + CHG_C_m + CHH_C_m, CG_mC_p + CHG_mC_p + CHH_mC_p, CG_mC_m + CHG_mC_m + CHH_mC_m
			CG_b, CHG_b, CHH_b, mCG_b, mCHG_b, mCHH_b = CG_C_p + CG_C_m, CHG_C_p + CHG_C_m, CHH_C_p + CHH_C_m, CG_mC_p + CG_mC_m, CHG_mC_p + CHG_mC_m, CHH_mC_p + CHH_mC_m
			out = '{itr}\t{Cp}\t{Cm}\t{CGb}\t{CHGb}\t{CHHb}\t{temp}\n'.format(
            itr = i, Cp = dens(C_p, mC_p), Cm = 0 - dens(C_m, mC_m),
            CGb = dens(CG_b, mCG_b), CHGb = dens(CHG_b, mCHG_b), CHHb = dens(CHH_b, mCHH_b),
            temp = tmp
            )
			chr_stat.write(out)
		chr_stat.close()

		code='''
name='%s'
''' %(label) + \
                '''
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }
 if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

library('ggplot2')
library('reshape')
setwd("%s")
previous_theme  <- theme_set(theme_bw())
chr='%s'
df=read.table(paste('%s','/',chr,'_mCX_density.stat',sep=''),header=F)
df1<-melt(df[,1:3],id="V1")
p1<-ggplot(df1)+geom_point(aes(x=V1,y=value,colour=variable),size=0.9)+theme(legend.position="",panel.border = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(),axis.text=element_text(vjust = 0.8),axis.title=element_text(face="bold", vjust = 0.2), plot.title=element_text(face="bold",vjust=1))+xlab("")+ylab("")+ggtitle(paste(name,' (chr:',chr,')',sep=''))
p2<-ggplot(df)+geom_smooth(aes(V1,V4),color="#5EC277",se=F)+xlab("")+ylab("")+theme(legend.position="",panel.border = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(),axis.text=element_text(vjust = 0.8),axis.title=element_text(face="bold", vjust = 0.2), plot.title=element_text(face="bold",vjust=1))+ggtitle("mCG/CG")
p3<-ggplot(df)+geom_smooth(aes(V1,V5),color="#5753A5",se=F)+xlab("")+ylab("")+theme(legend.position="",panel.border = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(),axis.text=element_text(vjust = 0.8),axis.title=element_text(face="bold", vjust = 0.2), plot.title=element_text(face="bold",vjust=1))+ggtitle("mCHG/CHG")
p4<-ggplot(df)+geom_smooth(aes(V1,V6),color="#FFA200",se=F)+xlab("")+ylab("")+theme(legend.position="",panel.border = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(),axis.text=element_text(vjust = 0.8),axis.title=element_text(face="bold", vjust = 0.2), plot.title=element_text(face="bold",vjust=1))+ggtitle("mCHH/CHH")
png(paste(name,'.',chr,"_mCX_density.png",sep = ''),type='cairo-png')
multiplot(p1,p2,p3,p4)
dev.off()
pdf(paste(name,'.',chr,"_mCX_density.pdf",sep=''))
multiplot(p1,p2,p3,p4)
dev.off()
''' % (den_dir, chr, den_dir)
		open(den_dir+'/'+chr+'_mCX_density.r','w').write(code)
		os.system('/PUBLIC/software/RNA/R-3.1.2/R-3.1.2/bin/Rscript --vanilla '+den_dir+'/'+chr+'_mCX_density.r')
	#	print 'chr density statistics have been finished\n'
	code='''
options(stringsAsFactors=FALSE)
library('ggplot2')
library('reshape2')
previous_theme <- theme_set(theme_bw())
sample_name = "%s"
list<-read.table("%s")
chrs<-list$V1
setwd("%s")
files = paste(chrs,pattern = "_mCX_density.stat",sep="")

for ( file_name in files)
{
  each_density <- read.delim(file_name, header=FALSE)[-1, 1:3]
  each_density <- melt(each_density,id=1)
  each_density$V1 <- as.numeric(as.character(each_density$V1))
  each_density$variable <- as.character(each_density$variable)
  each_density$value <- as.numeric(each_density$value)
  each_density$ChrLabel <- strsplit(file_name,"_")[[1]][1]
 if(!exists("all_density")){
  all_density <- each_density
 } else {
  all_density <- rbind(all_density, each_density)
 }
}
all_density$ChrLabel <- factor(all_density$ChrLabel, levels=chrs)


all_density = all_density[2:nrow(all_density),]

figure_pdf <- ggplot(data = all_density, aes(x = V1,y = value,color=variable))+ geom_hline(yintercept=0)+labs(title = sample_name)+ geom_point(size = 0.5) +facet_grid(ChrLabel ~ ., scale = "free")
figure_pdf <- figure_pdf + theme(legend.position =0,panel.background = element_blank(),axis.line = element_line()) + labs(x="Chromosome bins",y="Methylation density (%%)") + coord_cartesian(xlim=c(-0.02*max(all_density$V1),1.02*max(all_density$V1)))
forDisplayData <- all_density
figure_png <- ggplot(data = forDisplayData, aes(x = V1,y = value,color=variable))+ geom_hline(yintercept=0)+labs(title = sample_name)+ geom_point(size = 0.8) +facet_grid(ChrLabel ~ ., scale = "free")
figure_png <- figure_png + theme(legend.position =0,panel.background = element_blank(),axis.line = element_line()) + labs(x="Chromosome",y="Methylation density (%%)") + coord_cartesian(xlim=c(-0.02*max(forDisplayData$V1),1.02*max(forDisplayData$V1)))

figure_width = 12
figure_height = floor(length(chrs))
ggsave(paste(sample_name,"all_chr_mC_density_distribution.pdf",sep='.'),width=figure_width,height=figure_height,plot=figure_pdf)
ggsave(paste(sample_name,"all_chr_mC_density_distribution.png",sep='.'),type='cairo-png',width=8,height=figure_height,plot=figure_png)

''' % (label, argv.chrlist, den_dir)
#	print "all_chr_density began\n"
	open(den_dir+'/all_chr_mCX_density.r','w').write(code)
	os.system('/PUBLIC/software/RNA/R-3.1.2/R-3.1.2/bin/Rscript --vanilla '+den_dir+'/all_chr_mCX_density.r')

