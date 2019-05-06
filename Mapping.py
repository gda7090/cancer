import os,sys
import string

class Mapping:

	def __init__(self,pacientID,sampleID,rawDir,qcDir,alignDir,th,lib,genome_info,softwares):
		self.pacientID = pacientID
		self.sampleID = sampleID
		self.lib = lib
		self.qcDir = qcDir
		self.rawDir = rawDir
		self.alignDir = alignDir
		self.mergeBam = os.path.join(alignDir, self.sampleID+'.merge.bam')
		self.mergedBam = os.path.join(alignDir, self.sampleID+'.merged.bam')
		self.nodupBam = os.path.join(alignDir, self.sampleID+'.nodup.bam')
		self.splitBam = os.path.join(alignDir, self.sampleID+'.split.bam')
		self.splitedBam = os.path.join(alignDir, self.sampleID+'.splited.bam')
		self.discordBam = os.path.join(alignDir, self.sampleID+'.discord.bam')
		self.discordedBam = os.path.join(alignDir, self.sampleID+'.discorded.bam')
		self.realignBam = os.path.join(alignDir,self.sampleID+'.realn.bam')
		self.recalBam = os.path.join(alignDir, self.sampleID+'.realn.recal.bam')
		self.finalBam = os.path.join(alignDir, self.sampleID+'.final.bam')
		self.mpilegz = os.path.join(alignDir, self.sampleID+'.final.mpileup.gz')
		self.threads = th
		self.genome_info = genome_info
		self.chroms = [chr for chrset in genome_info['chr_subs'] for chr in chrset]
		self.softwares = softwares
		self.refData = genome_info['fasta']
		self.picardDIR = softwares['picard']
		self.gatkDIR = softwares['gatk']
		self.crestDIR = softwares['crest']
		self.speedseq = softwares['speedseq']
		self.ffpe = softwares['FFPE']

	def construct_fq(self,aLib,qc=True):
		fq1_gz = os.path.join(self.qcDir,self.sampleID+'/'+'_'.join([self.sampleID,aLib,'1.clean.fq.gz']))
		fq2_gz = os.path.join(self.qcDir,self.sampleID+'/'+'_'.join([self.sampleID,aLib,'2.clean.fq.gz']))
		fq1_nz = os.path.join(self.qcDir,self.sampleID+'/'+'_'.join([self.sampleID,aLib,'1.clean.fq']))
		fq2_nz = os.path.join(self.qcDir,self.sampleID+'/'+'_'.join([self.sampleID,aLib,'2.clean.fq']))
		fq1 = fq1_gz
		fq2 = fq2_gz
		sufix = 'fq.gz'
		if not os.path.exists(fq1_gz):
			fq1 = fq1_nz
			fq2 = fq2_nz
			sufix = 'fq'
		if os.path.exists(fq1_nz):
			fq1 = fq1_nz
			fq2 = fq2_nz
			sufix = 'fq'
		if qc:
			fq1 = fq1_nz
			fq2 = fq2_nz
			sufix = 'fq'
		if not qc and ((os.path.exists(fq1_gz) and os.path.exists(fq1_nz)) or (os.path.exists(fq1_gz) and os.path.exists(fq1_nz))):
			print '\033[1;35;40m'
			print "\t!!!Both fq file and fq.gz file exist for lane %s, only one type file is allowed!"%aLib
			print "\t!!!Makesure which one to use and remove the other manually! Or ask help for doublechecker."
			print '\033[0m'
			sys.exit(1)
		return fq1,fq2,sufix

	## aLib  =>  raw libID
	## bLib  =>  split libID

	def construct_rawBam(self,aLib,format='bam'):
		return os.path.join(self.alignDir, self.sampleID+'_'+aLib+'.'+format)

	def construct_sortBam(self,aLib):
		return os.path.join(self.alignDir, self.sampleID+'_'+aLib+'.sort.bam')

	def construct_splitBam(self,aLib,format='bam'):
		return os.path.join(self.alignDir, self.sampleID+'_'+aLib+'.split.'+format)

	def construct_discordBam(self,aLib,format='bam'):
		return os.path.join(self.alignDir, self.sampleID+'_'+aLib+'.discord.'+format)

	def construct_splitFq(self,bLib,sufix='fq'):
		fq1 = os.path.join(self.alignDir,'_'.join([self.sampleID,bLib,'1.clean.'+sufix]))
		fq2 = os.path.join(self.alignDir,'_'.join([self.sampleID,bLib,'2.clean.'+sufix]))
		return fq1,fq2

	def split_fastq(self,aLib,qc=True,num=1):
		order = 'order split_fastq_%s_%s after qc_%s_%s' % (self.sampleID,aLib,self.sampleID,aLib)
		fq1,fq2,sufix = self.construct_fq(aLib,qc=qc)
		if num == 1:
			splitFq1,splitFq2 = self.construct_splitFq(aLib+'-1',sufix=sufix)
			return 'cd %s && \\\nln -sf %s %s && \\\nln -sf %s %s\n' % (self.alignDir, fq1, splitFq1, fq2, splitFq2),order
		else:
			qcstat = os.path.join(self.qcDir,self.sampleID+'/'+'_'.join([self.sampleID,aLib+'.stat']))
			return 'python %s -fq1 %s -fq2 %s -o %s -p %s -n %d -qcstat %s\n' % \
			(os.path.join(self.softwares['summary'],'fastqbig2small.py'),fq1,fq2,self.alignDir,'_'.join([self.sampleID,aLib]),num,qcstat),order

	def bwa_mem(self,aLib,qc=True):
		rawBam = self.construct_rawBam(aLib)
		fq1,fq2,sufix = self.construct_fq(aLib,qc=qc)
#		order = 'order qc_%s_%s before bwa_mem_%s_%s' % (self.sampleID,aLib,self.sampleID,aLib)
		order = 'order bwa_mem_%s_%s after qc_%s_%s' % (self.sampleID,aLib,self.sampleID,aLib)
		return '  \\\n\t'.join(['bwa mem',
			'-t '+self.threads,
			'-M -R \'@RG\\tID:'+self.sampleID+'_'+aLib+'\\tSM:'+self.sampleID+'\\tLB:'+self.sampleID+'\\tPU:'+aLib+'\\tPL:illumina\\tCN:novogene\'',
			self.refData,
			fq1,
			fq2,
			'|samtools view -b -S -t',
			self.refData+'.fai',
			'- >',
			rawBam])+'\n',order

	def bwa_mem_samblaster(self,aLib,bLib,sufix='fq',qc=True,FFPE=False):
		rawSam = self.construct_rawBam(bLib,'sam')
		splitSam = self.construct_splitBam(bLib,'sam')
		discordSam = self.construct_discordBam(bLib,'sam')
		fq1,fq2 = self.construct_splitFq(bLib,sufix=sufix)
		softSam = self.construct_rawBam(bLib,'soft.sam')
		allSam = self.construct_rawBam(bLib,'all.sam')
		order = 'order bwa_mem_samblaster_%s_%s after split_fastq_%s_%s' % (self.sampleID,bLib,self.sampleID,aLib)
		step1 = '  \\\n\t'.join(['bwa mem',
			'-t '+self.threads,
			'-M -R \'@RG\\tID:'+self.sampleID+'_'+aLib+'\\tSM:'+self.sampleID+'\\tLB:'+self.sampleID+'\\tPU:'+aLib+'\\tPL:illumina\\tCN:novogene\'',
			self.refData,
			fq1,
			fq2,
			'|%s' % os.path.join(self.speedseq,'bin','samblaster'),
			'--splitterFile %s' % splitSam,
			'--discordantFile %s' % discordSam,
			'> %s' % rawSam])
		step2 = 'rm -rf %s %s' % (fq1, fq2)
		if FFPE:
			step3 = 'python %s %s %s %s' %(os.path.join(self.ffpe,'filter.py'),rawSam,softSam,allSam)
			step4 = 'mv -f %s %s' %(allSam,rawSam)
			return ' && \\\n'.join([step1,step2,step3,step4]),order
		return ' && \\\n'.join([step1,step2]),order

	def sambamba_sort_bam(self,bLib):
		rawSam = self.construct_rawBam(bLib,'sam')
		sortBam = self.construct_sortBam(bLib)
		order = 'order sambamba_sort_bam_%s_%s after bwa_mem_samblaster_%s_%s' % (self.sampleID,bLib,self.sampleID,bLib)
		step1 = '  \\\n\t'.join(['%s view -S -f bam -l 0' % os.path.join(self.speedseq,'bin','sambamba'),
			'%s' % rawSam,
			'|%s sort' % os.path.join(self.speedseq,'bin','sambamba'),
			'-t 5 -m 6G',
			'--tmpdir %s' % self.alignDir,
			'-o %s /dev/stdin' % (sortBam)])
		step2 = '%s index %s' % (os.path.join(self.speedseq,'bin','sambamba'),sortBam)
		step3 = 'rm -rf %s' % rawSam
		return ' && \\\n'.join([step1,step2,step3]),order
		
	def sambamba_sort_split(self,bLib):
		splitSam = self.construct_splitBam(bLib,'sam')
		splitBam = self.construct_splitBam(bLib,'bam')
		discordSam = self.construct_discordBam(bLib,'sam')
		discordBam = self.construct_discordBam(bLib,'bam')
		order = 'order sambamba_sort_split_%s_%s after bwa_mem_samblaster_%s_%s' % (self.sampleID,bLib,self.sampleID,bLib)
		step1 = '  \\\n\t'.join(['%s view' % os.path.join(self.speedseq,'bin','sambamba'),
				'-S -f bam -l 0 %s' % splitSam,
				'|%s sort' % os.path.join(self.speedseq,'bin','sambamba'),
				'-t 4 -m 1G',
				'--tmpdir=%s' % self.alignDir,
				'-o %s /dev/stdin' % splitBam])
		step2 = '  \\\n\t'.join(['%s view' % os.path.join(self.speedseq,'bin','sambamba'),
				'-S -f bam -l 0 %s' % discordSam,
				'|%s sort' % os.path.join(self.speedseq,'bin','sambamba'),
				'-t 4 -m 1G',
				'--tmpdir=%s' % self.alignDir,
				'-o %s /dev/stdin' % discordBam])
		step3 = '%s index %s' % (os.path.join(self.speedseq,'bin','sambamba'),splitBam)
		step4 = '%s index %s' % (os.path.join(self.speedseq,'bin','sambamba'),discordBam)
		step5 = 'rm -rf %s %s' % (splitSam,discordSam)
		return ' && \\\n'.join([step1,step2,step3,step4,step5])+'\n',order


	def gzip_fq(self,aLib,libs,qc=True):
		fq1,fq2,sufix = self.construct_fq(aLib,qc=qc)
		order = ['order gzip_fq_%s_%s after bwa_mem_samblaster_%s_%s' % (self.sampleID, aLib, self.sampleID, each) for each in libs]
#		order = 'order finalbam_%s before gzip_fq_%s_%s' % (self.sampleID, self.sampleID, aLib)
		step0 = 'cd %s' % os.path.join(self.qcDir,self.sampleID)
		if fq1.endswith('fq'):
			step1 = 'gzip -f %s' % fq1
			step2 = 'gzip -f %s' % fq2
			step3 = 'md5sum %s > %s' % (os.path.basename(fq1)+'.gz', fq1+'.gz.MD5.txt')
			step4 = 'md5sum %s > %s' % (os.path.basename(fq2)+'.gz', fq2+'.gz.MD5.txt')
			return '\n'.join([step0,step1,step2,step3,step4]),order
		else:
			step1 = 'md5sum %s > %s' % (os.path.basename(fq1), fq1+'.MD5.txt')
			step2 = 'md5sum %s > %s' % (os.path.basename(fq2), fq2+'.MD5.txt')
			return ' && \\\n'.join([step0,step1,step2]),order

	def sambamba_mergebam(self,libs,otherBam=''):
		order = []
		bams = []
		rmbams = []
		if otherBam:
			if otherBam == self.mergeBam: 
				os.rename(self.mergeBam,self.mergedBam)
				os.rename(self.mergeBam+'.bai',self.mergedBam+'.bai')
			bams.append(self.mergedBam)
			rmbams.append(self.mergedBam)
		for k in libs:
			bams.append(self.construct_sortBam(k))
			rmbams.append(self.construct_sortBam(k))
			order.append('order sambamba_mergebam_%s after sambamba_sort_bam_%s_%s' % (self.sampleID,self.sampleID,k))
		
		step3 = 'rm -f %s %s\n' % (' '.join(rmbams),' '.join([b+'.bai' for b in rmbams]))
		if len(bams) > 1:
			return '  \\\n\t'.join(['%s merge' % os.path.join(self.speedseq,'bin','sambamba'),
				self.mergeBam]+bams)+' && \\\n' + \
				'%s index %s' % (os.path.join(self.speedseq,'bin','sambamba'),self.mergeBam) + \
				' && \\\n'+step3,order
		if len(bams) == 1:
			return '  \\\n\t'.join(['mv -f',
				bams[0],
				self.mergeBam])+' && \\\n'+'mv -f %s.bai %s.bai\n'%(bams[0],self.mergeBam),order

	def sambamba_mergesplit(self,libs,sBam='',dBam=''):
		order = []
		split_bams = []
		discord_bams = []
		rmbams = []
		if sBam:
			if sBam == self.splitBam: 
				os.rename(self.splitBam,self.splitedBam)
				os.rename(self.splitBam+'.bai',self.splitedBam+'.bai')
			split_bams.append(self.splitedBam)
			rmbams.append(self.splitedBam)
		if dBam:
			if dBam == self.discordBam: 
				os.rename(self.discordBam,self.discordedBam)
				os.rename(self.discordBam+'.bai',self.discordedBam+'.bai')
			discord_bams.append(self.discordedBam)
			rmbams.append(self.discordedBam)
		for k in libs:
			split_bams.append(self.construct_splitBam(k))
			discord_bams.append(self.construct_discordBam(k))
			rmbams.append(self.construct_splitBam(k))
			rmbams.append(self.construct_discordBam(k))
			order.append('order sambamba_mergesplit_%s after sambamba_sort_split_%s_%s' % (self.sampleID,self.sampleID,k))
		
		step3 = 'rm -f %s %s\n' % (' '.join(rmbams),' '.join(['%s.bai'%x for x in rmbams]))
		if len(split_bams) > 1:
			step1 = '  \\\n\t'.join(['%s merge' % os.path.join(self.speedseq,'bin','sambamba'),
				self.splitBam]+split_bams)+' && \\\n' + \
				'%s index %s' % (os.path.join(self.speedseq,'bin','sambamba'),self.splitBam)
			step2 = '  \\\n\t'.join(['%s merge' % os.path.join(self.speedseq,'bin','sambamba'),
				self.discordBam]+discord_bams)+' && \\\n' + \
				'%s index %s' % (os.path.join(self.speedseq,'bin','sambamba'),self.discordBam)
			return ' && \\\n'.join([step1,step2,step3])+'\n',order
		if len(split_bams) == 1:
			return ' \\\n\t'.join(['mv -f',split_bams[0], self.splitBam]) + ' && \\\n' + \
				' \\\n\t'.join(['mv -f',discord_bams[0], self.discordBam]) + ' && \\\n' + \
				' \\\n\t'.join(['mv -f %s.bai %s.bai && \\\n' % (split_bams[0], self.splitBam)]) + \
				' \\\n\t'.join(['mv -f %s.bai %s.bai\n' % (discord_bams[0], self.discordBam)]),order

	def sambamba_markdup(self,bamcount=1,strag='WGS',FFPE=False):
		order = 'order sambamba_markdup_%s after sambamba_mergebam_%s' % (self.sampleID,self.sampleID)
		if FFPE:
			METRICS_FILE = self.mergeBam[:-4]+'.metrics'
			step1 = '  \\\n\t'.join(['java -Xmx10g -jar '+os.path.join(self.picardDIR,'MarkDuplicates.jar'),
				'TMP_DIR=\''+self.alignDir+'\'',
				'INPUT=\''+self.mergeBam+'\'',
				'OUTPUT=\''+self.nodupBam+'\'',
				'METRICS_FILE=\''+METRICS_FILE+'\'',
				'VALIDATION_STRINGENCY=SILENT',
				'ASSUME_SORTED=true',
				'CREATE_INDEX=true',
				'REMOVE_DUPLICATES=true',
				'MAX_RECORDS_IN_RAM=4000000'])
			step2 = 'mv -f %s.bai %s.bai' % (self.nodupBam[:-4],self.nodupBam)
			return ' && \\\n'.join([step1, step2])+'\n',order	
		else:				
			if bamcount == 1:
				step1 = 'cp %s %s' % (self.mergeBam, self.nodupBam)
				step2 = 'cp %s.bai %s.bai' % (self.mergeBam, self.nodupBam)
				return ' && \\\n'.join([step1, step2])+'\n',order
			else:
				if 'WGS' in strag:
					step1 = '  \\\n\t'.join(['python %s' % os.path.join(self.softwares['summary'],'mark_dup.py'),
						'-b %s' % self.mergeBam,
						'-o %s' % self.alignDir,
						'-s %s' % self.sampleID])
					step2 = '%s index %s' % (os.path.join(self.speedseq,'bin','sambamba'),self.nodupBam)
					return ' && \\\n'.join([step1, step2])+'\n',order
				else:
					step1 = ' \\\n\t'.join(['%s markdup' % os.path.join(self.speedseq,'bin','sambamba'),
						'--tmpdir %s' % self.alignDir,self.mergeBam, self.nodupBam])
					step2 = '%s index %s' % (os.path.join(self.speedseq,'bin','sambamba'),self.nodupBam)
					return ' && \\\n'.join([step1, step2])+'\n',order
	def gatk_realign(self):
		order = 'order gatk_realign_%s after sambamba_markdup_%s' % (self.sampleID,self.sampleID)
		step1 = '  \\\n\t'.join([self.softwares['java7']+' -Xmx10g -jar '+os.path.join(self.gatkDIR,'GenomeAnalysisTK.jar'),
			'-T RealignerTargetCreator',
			'-nt 6',
			'-R '+self.refData,
			'-I '+self.nodupBam,
			'-known %s' % self.genome_info['1000indel'],
			'-o '+self.nodupBam+'.intervals'])
		step2 = '  \\\n\t'.join([self.softwares['java7']+' -Xmx10g -jar '+os.path.join(self.gatkDIR,'GenomeAnalysisTK.jar'),
			'-T IndelRealigner',
			'-R '+self.refData,
			'-I '+self.nodupBam,
			'-targetIntervals '+self.nodupBam+'.intervals',
			'-known %s' % self.genome_info['1000indel'],
			'-o '+self.realignBam])
		step3 = 'mv -f %s.bai %s.bai ' % (self.realignBam[:-4],self.realignBam)
		return ' && \\\n'.join([step1,step2,step3])+'\n',order
		
	
	def gatk_bqsr(self):
		order = 'order gatk_bqsr_%s after gatk_realign_%s' % (self.sampleID,self.sampleID)
		step1 = '  \\\n\t'.join([self.softwares['java7']+' -Xmx10g -jar '+os.path.join(self.gatkDIR,'GenomeAnalysisTK.jar'),
			'-T BaseRecalibrator',
			'-nct 6',
			'-R '+self.refData,
			'-I '+self.realignBam,
			'-rf BadCigar',
			'-knownSites %s' % self.genome_info['dbsnp'],
			'-knownSites %s' % self.genome_info['1000indel'],
			'-o '+os.path.join(self.alignDir,'recal_data.table')])
			
		step2 = '  \\\n\t'.join([self.softwares['java7']+' -Xmx10g -jar '+os.path.join(self.gatkDIR,'GenomeAnalysisTK.jar'),
			'-T PrintReads',
			'-R '+self.refData,
			'-I '+self.realignBam,
			'-BQSR '+os.path.join(self.alignDir,'recal_data.table'),
			'-o '+self.recalBam])
		step3 = 'mv -f %s.bai %s.bai' % (self.recalBam[:-4],self.recalBam)
		return ' && \\\n'.join([step1,step2,step3])+'\n',order	
	
	def finalbam(self,reGATK=False):
		if reGATK:
			step1 = 'ln -sf %s %s' % (self.recalBam,self.finalBam)
			step2 = 'ln -sf %s.bai %s.bai' % (self.recalBam,self.finalBam)
			step3 = 'rm -f %s %s.bai %s.intervals %s' % (self.realignBam,self.realignBam,self.mergeBam,os.path.join(self.alignDir,'recal_data.table'))
			order = 'order finalbam_%s after gatk_bqsr_%s' % (self.sampleID,self.sampleID)
			return ' && \\\n'.join([step1,step2,step3])+'\n',order
		else:
			step1 = 'ln -sf %s %s' % (self.nodupBam,self.finalBam)
			step2 = 'ln -sf %s.bai %s.bai' % (self.nodupBam,self.finalBam)
			order = 'order finalbam_%s after sambamba_markdup_%s' % (self.sampleID,self.sampleID)
			return ' && \\\n'.join([step1,step2])+'\n',order
		
	def mpileup(self):
		order = 'order mpileup_%s after finalbam_%s' % (self.sampleID,self.sampleID)
		return ' \\\n\t'.join(['python %s' % os.path.join(self.softwares['advbin'],'mpileup.py'),
			'-bam %s' % self.finalBam,
			'-o %s' % self.alignDir,
			'-ref %s' % self.refData,
			'-sample %s' % self.sampleID]),order
#		return ' \\\n\t'.join([self.softwares['samtools.0.1.18']+' mpileup -q 1 -f '+self.refData,
#			self.finalBam,'|gzip -f - >'+self.mpilegz])+'\n',order

	def extractSoftClip_sub(self,chrs,num):
		cmds = []
		order = 'order extractSoftClip_sub%s_%s after sambamba_mergesplit_%s' % (num,self.sampleID,self.sampleID)
		for chr in chrs:
			cover = os.path.join(self.alignDir,self.sampleID+'.'+chr+'.cover')
			sclip = os.path.join(self.alignDir,self.sampleID+'.'+chr+'.sclip.txt')
			step1 = 'if [ -e %s ]; then rm -rf %s; fi' % (cover,cover)
			step2 = 'if [ -e %s ]; then rm -rf %s; fi' % (sclip,sclip)
			step3 = ' \\\n\t'.join(['perl '+os.path.join(self.crestDIR,'CREST','extractSClip.pl'),
				'-o '+self.alignDir,
				'-i '+self.splitBam,
				'-ref_genome '+self.refData,
				'-r %s -p %s' % (chr,self.sampleID)])
			cmds += [step1,step2,step3]
		return ' && \\\n'.join(cmds),order

	def catSoftClip(self,subs,wtffpe):
			order = ['order catSoftClip_%s after extractSoftClip_sub%s_%s' % \
				(self.sampleID,each,self.sampleID) for each in subs]
			final_cover = os.path.join(self.alignDir,self.sampleID+'.cover')
			final_sclip = os.path.join(self.alignDir,self.sampleID+'.sclip.txt')
			coverlist = [os.path.join(self.alignDir,self.sampleID+'.'+chr+'.cover') for chr in self.chroms]
			scliplist = [os.path.join(self.alignDir,self.sampleID+'.'+chr+'.sclip.txt') for chr in self.chroms]
			step5 = 'echo "No filter cover!!!"'
			if wtffpe == 'Y':
				step5 = 'count=$(less %s|wc -l)' % final_cover +'\n'+ 'if [ $count -gt 10000000 ];then python %s -i %s -o %s && mv %s %s && mv %s %s;else echo "No filter cover";fi' % ('Cancer/script/Cover_Filter.py',final_cover,os.path.join(self.alignDir,self.sampleID+'.flt.cover'),final_cover,os.path.join(self.alignDir,self.sampleID+'.raw.cover'),os.path.join(self.alignDir,self.sampleID+'.flt.cover'),final_cover)
			step1 = 'cat '+' \\\n\t'.join(coverlist) + ' \\\n\t>'+final_cover
			step2 = 'cat '+' \\\n\t'.join(scliplist) + ' \\\n\t>'+final_sclip
			step3 = 'rm -rf '+' \\\n\t'.join(coverlist)
			step4 = 'rm -rf '+' \\\n\t'.join(scliplist)
			return ' && \\\n'.join([step1,step3,step4,step5])+'\n',order
