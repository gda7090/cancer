import os
	
class SV:
	def __init__(self, pacientID, sampleID, genome_info, TR, projDir, svDir, softwares):
		self.pacientID = pacientID
		self.sampleID = sampleID
		self.mapDir = os.path.join(projDir,'Mapping',sampleID)
		self.statDir = os.path.join(projDir,'Alnstat',sampleID)
		self.svDIr = svDir
		self.breakdancerDir = softwares['breakdancer']
		self.crestDir = softwares['crest']
		self.cnvnatorDir = softwares['cnvnator']
		self.speedseqDir = softwares['speedseq']
		self.ref_chr = genome_info['bychr']
		self.humandbDIR = genome_info['annovardb']
		self.advDir = softwares['advbin']
		self.TR = TR
		self.genome_info = genome_info
		self.refData = genome_info['fasta']
		self.refData2bit = genome_info['fa2bit']
		self.softwares = softwares
		self.annovar = softwares['annovar']
		self.snvdb = 'GeneName,refGene,genomicSuperDups,gff3,avsnp142,cosmic70,clinvar_20150330,gwasCatalog,1000g2014oct_Chinese,1000g2014oct_eas,1000g2014oct_all,esp6500siv2_all,exac03_ALL_EAS,ljb26_sift,ljb26_pp2hvar,ljb26_pp2hdiv,lib26_mt,gerp++gt2,caddgt10/caddindel'
		self.svdb = 'GeneName,refGene,Gencode,cpgIslandExt,cytoBand,dgvMerged'
		self.chroms = [chr for chrset in genome_info['chr_subs'] for chr in chrset]	
		
	def bam2cfg(self):
		#order = 'order finalbam_%s before bam2cfg_%s' % (self.sampleID,self.sampleID)
		order = 'order bam2cfg_%s after finalbam_%s' % (self.sampleID,self.sampleID)
		bam = os.path.join(self.mapDir,self.sampleID+'.final.bam')
		breakdanceroutdir = os.path.join(self.svDIr,'breakdancer')
		return  'cd '+breakdanceroutdir+'\n'+'  \\\n\t'.join(['perl '+os.path.join(self.breakdancerDir,'bam2cfg.pl'),
			'-g -h',
			bam,
			'> '+os.path.join(breakdanceroutdir,self.sampleID+'.breakdancer.cfg')]),order
		
		
	def breakdancerSV(self):
		#order = 'order bam2cfg_%s before breakdancerSV_%s' % (self.sampleID,self.sampleID)
		order = 'order breakdancerSV_%s after bam2cfg_%s' % (self.sampleID,self.sampleID)
		bam = os.path.join(self.mapDir,self.sampleID+'.final.bam')
		breakdanceroutdir = os.path.join(self.svDIr,'breakdancer')
		step0 = 'cd %s' % breakdanceroutdir
		step1 = '  \\\n\t'.join([os.path.join(self.breakdancerDir,'breakdancer-max'),
			'-h',
			'-d '+os.path.join(breakdanceroutdir,self.sampleID+'.breakdancer.SV-supporting'),
			'-g '+os.path.join(breakdanceroutdir,self.sampleID+'.breakdancer.bed'),
			os.path.join(breakdanceroutdir,self.sampleID+'.breakdancer.cfg'),
			'> '+os.path.join(breakdanceroutdir,self.sampleID+'.breakdancer.txt')])
		step2 = '  \\\n\t'.join(['var_sv_breakdancer.filter.pl',
			'-g M -n 6',
			'-a '+os.path.join(breakdanceroutdir,self.sampleID+'.breakdancer.txt'),
			'> '+os.path.join(breakdanceroutdir,self.sampleID+'.breakdancer.flt.txt')])
		step3 = '  \\\n\t'.join(['var_sv_breakdancer.toGff.pl',
			os.path.join(breakdanceroutdir,self.sampleID+'.breakdancer.flt.txt'),
			'> '+os.path.join(breakdanceroutdir,self.sampleID+'.breakdancer.sv.gff')])
		return ' && \\\n'.join([step0,step1,step2,step3]),order

	def breakdancerSvAnno(self):
		breakdanceroutdir = os.path.join(self.svDIr,'breakdancer')
		#order = 'order breakdancerSV_%s before breakdancerSvAnno_%s' % (self.sampleID,self.sampleID)
		order = 'order breakdancerSvAnno_%s after breakdancerSV_%s' % (self.sampleID,self.sampleID)
		sv_db = self.svdb
		if self.genome_info['annovarbuild'] == 'mm10':
			sv_db = 'GeneName,refGene,genomicSuperDups,gff3,snp142chr'
		f = os.path.join(breakdanceroutdir,self.sampleID+'.breakdancer.sv.gff')
		step1 = '  \\\n\t'.join([self.annovar,
			'-z %s' % sv_db,
			'-r %s' % self.refData,
			'-u %s' % self.genome_info['annovarbuild'],
			'-d '+self.humandbDIR,
			f,
			self.sampleID])
		step2 = '  \\\n\t'.join(['sv_cnv.stat.cancer.pl',
			'-s '+self.sampleID,
			f,
			'> '+os.path.join(breakdanceroutdir,self.sampleID+'.breakdancer.sv.stat.xls')])
		return ' && \\\n'.join([step1, step2]), order
	
	def crestSV(self):
		order = ['order crestSV_%s after catSoftClip_%s' % (self.sampleID,self.sampleID),
			'order crestSV_%s after finalbam_%s' % (self.sampleID,self.sampleID)]
		cov = os.path.join(self.mapDir,self.sampleID+'.cover')
		crestoutDir = os.path.join(self.svDIr,'crest')
		bam = os.path.join(self.mapDir,self.sampleID+'.final.bam')
		step1 = ' \\\n\t'.join(['perl '+os.path.join(self.crestDir,'crest_sv_calling.pl'),
			'-cov ' + cov,
			'-outDir ' + crestoutDir,
			'-tumorBam ' + bam,
			'-sampleID ' + self.sampleID+'.crest.sv',
			'-regionList %s' % self.genome_info['chrbed'],
			'-ref ' + self.refData,
			'-bit ' + self.refData2bit])
		return step1+'\n',order

	def crest_SV_sub(self,chroms,i):
		order = ['order crestSV_sub%s_%s after catSoftClip_%s' % (i,self.sampleID,self.sampleID),
			'order crestSV_sub%s_%s after finalbam_%s' % (i,self.sampleID,self.sampleID)]
		cov = os.path.join(self.mapDir,self.sampleID+'.cover')
		crestoutDir = os.path.join(self.svDIr,'crest')
		bam = os.path.join(self.mapDir,self.sampleID+'.final.bam')
		cmds = []
		bychr_suf = '.region.bed'
		if self.genome_info['build'] != 'b37':
			bychr_suf = '.bed'
		for each_chr in chroms:
			cmd_tmp = ' \\\n\t'.join(['perl '+os.path.join(self.crestDir,'crest_sv_calling.pl'),
			'-cov ' + cov,
			'-outDir ' + crestoutDir,
			'-tumorBam ' + bam,
			'-sampleID ' + self.sampleID+'.crest.sv'+'.'+each_chr,
			'-regionList %s' % os.path.join(self.genome_info['bychr'],each_chr + bychr_suf),
			'-ref ' + self.refData,
			'-bit ' + self.refData2bit])
			cmds.append(cmd_tmp)
		return ' && \\\n'.join(cmds)+'\n',order

	def crestSVann(self):
		crestoutDir = os.path.join(self.svDIr,'crest')
		#order = 'order crestSV_%s before crestSvAnno_%s' % (self.sampleID,self.sampleID)
		order = 'order crestSVann_%s after crestSV_%s' % (self.sampleID,self.sampleID)
		sv_db = self.svdb
		if self.genome_info['annovarbuild'] == 'mm10':
			sv_db = 'GeneName,refGene,genomicSuperDups,gff3,snp142chr'
		step1 = '  \\\n\t'.join(['sh '+self.annovar,
			'-t %s' % SVType,
			os.path.join(crestoutDir,self.sampleID+'.crest.sv.gff'),
			self.sampleID])
		step2 = ' \\\n\t'.join(['sv_cnv.stat.cancer.pl',
			'-s '+self.sampleID,
			os.path.join(crestoutDir,self.sampleID+'.crest.sv.hg19_multianno.xls'),
			'> '+os.path.join(crestoutDir,self.sampleID+'.crest.sv.stat.xls')])
		return ' && \\\n'.join([step1,step2]),order

	def crest_catannovar(self,subs):
		crestoutDir = os.path.join(self.svDIr,'crest')
		ann = os.path.join(crestoutDir,self.sampleID+'.crest.sv.hg19_multianno.xls')
		stat = os.path.join(crestoutDir,self.sampleID+'.crest.sv.stat.xls')
		order = ['order crestSVcatann_%s after crestSV_sub%s_%s' % \
			(self.sampleID,each,self.sampleID) for each in subs]
		sv_db = self.svdb
		if self.genome_info['annovarbuild'] == 'mm10':
			sv_db = 'GeneName,refGene,genomicSuperDups,gff3,snp142chr'
		chr_gff = [os.path.join(crestoutDir,self.sampleID+'.crest.sv.%s.gff' % chr) for chr in self.chroms]
		chr_txt = [os.path.join(crestoutDir,self.sampleID+'.crest.sv.%s.predSV.txt' % chr) for chr in self.chroms]
		#step1 = 'cat '+'\\\n\t'.join(chr_gff)+'\\\n\t>'+ os.path.join(crestoutDir,(self.sampleID+'.crest.sv.gff'))
		step1 = 'cat '+'\\\n\t'.join(chr_txt)+'\\\n\t>'+ os.path.join(crestoutDir,(self.sampleID+'.crest.sv.predSV.txt'))
		step2 = 'perl share/software/CREST/var_sv_CREST.toGff.pl ' + os.path.join(crestoutDir,(self.sampleID+'.crest.sv.predSV.txt')) +' > '+os.path.join(crestoutDir,(self.sampleID+'.crest.sv.gff'))
		step3 = '  \\\n\t'.join([self.annovar,
		'-r %s' % self.refData,
		'-u %s' % self.genome_info['annovarbuild'],
		'-z %s' % sv_db,
		'-d %s' % self.humandbDIR,
		os.path.join(crestoutDir,self.sampleID+'.crest.sv.gff'),
		self.sampleID])
		step4 = 'sv_cnv.stat.cancer.pl -s %s %s > %s' %(self.sampleID,ann,stat)
		step5 = 'rm -rf %s \\\n  %s' % ('\\\n  '.join(chr_gff),'\\\n  '.join(chr_txt))
		return ' && \\\n'.join([step1,step2,step3,step4,step5])+'\n',order

	def lumpySV(self):
		order = ['order lumpySV_%s after sambamba_mergesplit_%s' % (self.sampleID,self.sampleID),
			'order lumpySV_%s after finalbam_%s' % (self.sampleID,self.sampleID)]
		bam = os.path.join(self.mapDir,self.sampleID+'.final.bam')
		sbam = os.path.join(self.mapDir,self.sampleID+'.split.bam')
		dbam = os.path.join(self.mapDir,self.sampleID+'.discord.bam')
		step1 = ' \\\n\t'.join([os.path.join(self.speedseqDir,'bin','lumpyexpress'),
			'-B %s' % bam,
			'-S %s' % sbam,
			'-D %s' % dbam,
			'-o %s' % os.path.join(self.svDIr,'lumpy',self.sampleID+'.lumpy.sv.tmp.vcf'),
			'-x %s' % self.genome_info['lumpy-exclude'],
			'-T %s' % os.path.join(self.svDIr,'lumpy'),
			'-K %s' % os.path.join(self.speedseqDir,'bin','speedseq.config'),
			'-P -v -k'])
		step2 = 'grep -v IMPRECISE %s > %s' % (os.path.join(self.svDIr,'lumpy',self.sampleID+'.lumpy.sv.tmp.vcf'),os.path.join(self.svDIr,'lumpy',self.sampleID+'.lumpy.sv.vcf'))
		step3 = ' \\\n\t'.join(['python %s/lumpy_vcf2gff.py' % self.advDir,
			os.path.join(self.svDIr,'lumpy',self.sampleID+'.lumpy.sv.vcf'),
			os.path.join(self.svDIr,'lumpy',self.sampleID+'.lumpy.sv.gff')])
		
		return ' && \\\n'.join([step1,step2,step3]),order

	def lumpySVann(self):
		order = 'order lumpySVann_%s after lumpySV_%s' % (self.sampleID,self.sampleID)
		sv_db = self.svdb
		if self.genome_info['annovarbuild'] == 'mm10':
			sv_db = 'GeneName,refGene,genomicSuperDups,gff3,snp142chr'
		step1 = ' \\\n\t'.join([self.annovar,
			'-r %s' % self.refData,
			'-u %s' % self.genome_info['annovarbuild'],
			'-z %s' % sv_db,
			'-d %s' % self.humandbDIR,
			os.path.join(self.svDIr,'lumpy',self.sampleID+'.lumpy.sv.gff'),
			self.sampleID])
		step2 = ' \\\n\t'.join(['sv_cnv.stat.cancer.pl',
			'-s '+self.sampleID,
			os.path.join(self.svDIr,'lumpy',self.sampleID+'.lumpy.sv.hg19_multianno.xls'),
			'> '+os.path.join(self.svDIr,'lumpy',self.sampleID+'.lumpy.sv.stat.xls')])
		return ' && \\\n'.join([step1,step2]),order


	def cnvnatorCNV(self):
		#order = 'order finalbam_%s before cnvnatorCNV_%s' % (self.sampleID,self.sampleID)
		order = 'order cnvnatorCNV_%s after finalbam_%s' % (self.sampleID,self.sampleID)
		cnvnatoroutDir = os.path.join(self.svDIr,'cnvnator')
		s = os.path.join(cnvnatoroutDir,self.sampleID+'.cnvnator')
		step0 = 'cnvnator -root %s.root -tree %s -unique' % (s,os.path.join(self.mapDir,self.sampleID+'.final.bam'))
		step1 = 'cnvnator -root %s.root -his 100 -d %s' % (s, self.ref_chr)
		step2 = 'cnvnator -root %s.root -stat 100' % s
		step3 = 'cnvnator -root %s.root -partition 100' % s
		step31 = 'cnvnator -root %s.root -call 100 >%s.raw' % (s, s)
		step4 = 'cat %s.raw | grep -v "\\(WARN\\)\\|\\(==\\)" > %s.txt' % (s, s)
		step5 = 'awkopt=\'{feat=$4; if ($4=="deletion") feat="Deletion"; else if \
		($4=="duplication") feat="Duplication"; if ($5>=\'0\') print $1"\\tCNVnator\\t"feat"\\t"$2"\\t"$3"\\t.\\t.\\t.\\tSize="$5";RD="$6";q0="$7};\''
		step6 = ' \\\n'.join(['perl -F"\\t" -ane \'chomp @F;$F[1]=~s/^chr//;$F[1]=~s/[:-]/\\t/g;print join("\\t",$F[1],$F[0],$F[2],$F[3],$F[8])."\\n";\' '+s+'.txt',
			'| awk "$awkopt" | sort -k1,1 -k4n -k5n -u',
			'| awk -v OFS="\\t" \'{print $1,$4-1,$5,$0}\'',
			'| intersectBed -a stdin -b '+self.genome_info['Nblock']+' -f 0.5 -r -v',
			'| perl -F"\\t" -ane \'chomp @F;shift @F;shift @F;shift @F;printf "%s\\n",join("\\t",@F);\'',
			'| awk \'$1~/^[0-9XY]+$/\'',
			'| perl -F"\\t" -ane \'/Size=(.+?);RD=(.+?);q0=(.+)/;if(/^#/|| ($1>=1000 && ($2<0.6 || $2>1.4 ) && $3<0.5 && $F[0]!~/[Y]/) ){print "$_" }\'',
			'| perl -ane \'BEGIN{$CNVID=0} $CNVID++; chomp; @F=split /\\t/; print "$_;CNVID=$CNVID;CNVType=$F[2]\\n";\'',
			'>'+s+'.cnv.gff'])
		step7 = 'rm -rf %s.root' % s
		return ' && \\\n'.join([step0,step1,step2,step3,step31,step4,step5,step6,step7]), order

	def cnvnatorCNVann(self):
		cnvnatoroutDir = os.path.join(self.svDIr,'cnvnator')
		s = os.path.join(cnvnatoroutDir,self.sampleID+'.cnvnator')
		order = 'order cnvnatorCNVann_%s after finalbam_%s' % (self.sampleID,self.sampleID)
		sv_db = self.svdb
		if self.genome_info['annovarbuild'] == 'mm10':
			sv_db = 'GeneName,refGene,genomicSuperDups,gff3,snp142chr'
		step1 = ' \\\n\t'.join([self.annovar,
			'-z %s' % sv_db,
			'-r %s' % self.refData,
			'-u %s' % self.genome_info['annovarbuild'],
			'-t CNVType',
			'-d '+self.humandbDIR,
			s+'.cnv.gff',
			self.sampleID])
		step2 = '\\\n\t'.join(['sv_cnv.stat.cancer.pl',
			'-s '+self.sampleID,
			os.path.join(cnvnatoroutDir,self.sampleID+'.cnvnator.cnv.hg19_multianno.xls'),
			'> '+os.path.join(cnvnatoroutDir,self.sampleID+'.cnvnator.cnv.stat.xls')])
		return ' && \\\n'.join([step1,step2]),order
	
	def freecCNV(self,seq,sex='XY',step=1000):
		#order = ['order freec_cnv_%s after finalbam_%s' % (self.sampleID,self.sampleID),
		#	 'order freec_cnv_%s after depthStat_%s' % (self.sampleID,self.sampleID)]
		chr_flag = ''
		if not self.genome_info['build'].startswith('b'):
			chr_flag = ' --chr'
		freecoutDir = os.path.join(self.svDIr,'freec')
		config = '''
[general]
BedGraphOutput=FALSE
breakPointType = 4
chrFiles = %s
chrLenFile = %s
maxThreads= 6
outputDir = %s
ploidy = 2
gemMappabilityFile = %s
uniqueMatch=TRUE
window = %d
step = %d
sex=%s
telocentromeric = 50000
forceGCcontentNormalization = 1
[sample]
mateFile = %s
inputFormat = pileup
mateOrientation = FR
[BAF]
#to calculate BAF values, you need to provide mateFile in SAMtools pileup format
SNPfile = %s
minimalCoveragePerPosition = 5
shiftInQuality = 33
''' % (self.genome_info['bychr'],self.genome_info['chrlen'],freecoutDir,self.genome_info['mappability'],step*2,step,sex,os.path.join(self.mapDir,self.sampleID+'.final.mpileup.gz'),self.genome_info['dbsnp.txt'])
		open(os.path.join(freecoutDir,'freec.config'),'w').write(config)
		gff = os.path.join(freecoutDir,self.sampleID+'.freec.cnv.gff')
		order = 'order freecCNV_%s after mpileup_%s' % (self.sampleID,self.sampleID)

		step0 = 'sex=`if [ -f "%s" ];then awk \'{if(/F/){print "XX"}else{print "XY"}}\' %s;fi` && \\\n' % \
			(os.path.join(self.statDir,self.sampleID+'.gender'), os.path.join(self.statDir,self.sampleID+'.gender'))
		step0 += 'if [ $sex == "XX" ];then sed -i \'s/sex=XY/sex=XX/g\' %s ;sed -i \'s/sex = XY/sex = XX/g\' %s ;fi' % \
			(os.path.join(freecoutDir,'freec.config'), os.path.join(freecoutDir,'freec.config'))
		step1 = 'freec -conf %s' % os.path.join(freecoutDir,'freec.config')
		step2 = 'python %s/Freec.format.py -i %s -s 1000 -o %s -x $sex -b %s' \
			% (self.advDir,os.path.join(freecoutDir,self.sampleID+'.final.mpileup.gz_CNVs'), \
			os.path.join(freecoutDir,self.sampleID+'.CNV.txt'), self.genome_info['chrbed']) + chr_flag
		step3 = 'grep -v normal %s | awk -F "\\t" -v OFS="\\t" \'{{print $1,"FREEC",$7,$2,$3,".",".",".","CopyNumber="$5";Size="$4";CNVType="$7";Genotype="$8";GTConfidence="$9}}\' > %s' \
			% (os.path.join(freecoutDir,self.sampleID+'.CNV.txt'),gff)
		step4 = 'rm -rf %s %s' %(os.path.join(freecoutDir,'GC_profile.cnp'), \
			os.path.join(freecoutDir,self.sampleID+'.final.mpileup.gz_sample.cpn'))
		step5 = 'gzip -f %s %s' % (os.path.join(freecoutDir,self.sampleID+'.final.mpileup.gz_BAF.txt'), \
			os.path.join(freecoutDir,self.sampleID+'.final.mpileup.gz_ratio.txt'))
		return ' && \\\n'.join([step0,step1,step2,step3,step4,step5])+'\n',order

	def freecCNVann(self):
		order = 'order freecCNVann_%s after freecCNV_%s' % (self.sampleID,self.sampleID)
		freecoutDir = os.path.join(self.svDIr,'freec')
		ann = os.path.join(freecoutDir,self.sampleID+'.freec.cnv.hg19_multianno.xls')
		stat = os.path.join(freecoutDir,self.sampleID+'.freec.cnv.stat.xls')
		sv_db = self.svdb
		if self.genome_info['annovarbuild'] == 'mm10':
			sv_db = 'GeneName,refGene,genomicSuperDups,gff3,snp142chr'
		step1 = ' \\\n\t'.join(['sh '+self.annovar,
			'-t CNVType',
			os.path.join(freecoutDir,self.sampleID+'.freec.cnv.gff'),
			self.sampleID])
		step2 = 'sv_cnv.stat.cancer.pl -s %s %s > %s' % (self.sampleID,ann,stat)
		return ' && \\\n'.join([step1,step2]),order
