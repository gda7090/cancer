import os,sys

class QC:

	def __init__(self,sampleID,rawdir,qcDir,fastqs,softwares,databases):
		self.sampleID = sampleID
		self.rawdir = rawdir
		self.qcDir = qcDir
		self.fastqs = fastqs
		self.softwares = softwares
		self.databases = databases

	def rm_rRNA(self,alib,alane,rRNA_rate=10):
		rawfq1 = self.fastqs[alib][alane][0]
		rawfq2 = self.fastqs[alib][alane][1]
		adapt1 = rawfq1[:-5]+'adapter.list.gz'
		adapt2 = rawfq2[:-5]+'adapter.list.gz'
		prefix = '%s_%s_%s' % (self.sampleID,alib,alane)
		cmd = '''cd %s
gzip -dc %s |head -n 4000000 > %s_1.test.fq && \\
gzip -dc %s |head -n 4000000 > %s_2.test.fq && \\
Alignment/bowtie2-2.0.6/bowtie2 \\
  -x /PUBLIC/database/RNA/Med/Database/rRNA/hsa/hsa \\
  -1 %s_1.test.fq -2 %s_2.test.fq --end-to-end --sensitive -p 8 --phred33 --no-mixed -X 600 \\
  -S %s.test.rRNA.sam 2>%s.test.rRNA.stat.txt && \\
rm %s_1.test.fq %s_2.test.fq %s.test.rRNA.sam && \\
maprate=$(tail -n 1 %s.test.rRNA.stat.txt | awk '{print $1}' | awk -F '%%' '{print $1}') && \\
maprate=$(echo $maprate/1|bc) && \\
if [[ $maprate -ge %d ]]; then
  Alignment/bowtie2-2.0.6/bowtie2 \\
    -x /PUBLIC/database/RNA/Med/Database/rRNA/hsa/hsa \\
    -1 %s -2 %s --end-to-end --sensitive -p 8 --phred33 --no-mixed -X 600 \\
    --un-conc-gz %s.unmap.gz \\
    -S %s.rRNA.sam 2>%s.rRNA.stat.txt && \\
	rm %s.rRNA.sam && \\
    mv %s.unmap.1.gz %s_1.fq.gz && \\
    mv %s.unmap.2.gz %s_2.fq.gz
else
  ln -s %s %s_1.fq.gz
  ln -s %s %s_2.fq.gz
fi
'''%(os.path.join(self.rawdir,self.sampleID),rawfq1,prefix,rawfq2,prefix,prefix,prefix,prefix,prefix,prefix,prefix,prefix,prefix,rRNA_rate,rawfq1,rawfq2,prefix,prefix,prefix,prefix,prefix,prefix,prefix,prefix,rawfq1,prefix,rawfq2,prefix)
		if os.path.exists(adapt1) and os.path.exists(adapt2):
			cmd += 'ln -sf %s %s_1.adapter.list.gz\n' % (adapt1,prefix)
			cmd += 'ln -sf %s %s_2.adapter.list.gz\n' % (adapt2,prefix)
		return cmd

	def ln(self,alib,alane,rawopts='',fqcheck=False,rmadapter=False,cutadapter=False,rmdup=False,duplevel=15):
		rawfq1 = self.fastqs[alib][alane][0]
		rawfq2 = self.fastqs[alib][alane][1]
		adapt1 = rawfq1[:-5]+'adapter.list.gz'
		adapt2 = rawfq2[:-5]+'adapter.list.gz'
		cmd = []
		cmd.append('cd %s'%os.path.join(self.rawdir,self.sampleID))

		cutadapt_options = '-n 1 -m 100 -O 10 -e 0.2'
		
		if (not rawopts) and (not cutadapter) and (not rmadapter):######no rm,cut,rawopts####
			cmd.append('ln -sf %s %s_%s_%s_1.fq.gz'%(rawfq1,self.sampleID,alib,alane))
			cmd.append('ln -sf %s %s_%s_%s_2.fq.gz'%(rawfq2,self.sampleID,alib,alane))
#			if os.path.exists(adapt1) and os.path.exists(adapt2) and (not fqcheck):
#				cmd.append('ln -sf %s %s_%s_%s_1.adapter.list.gz'%(adapt1,self.sampleID,alib,alane))
#				cmd.append('ln -sf %s %s_%s_%s_2.adapter.list.gz'%(adapt2,self.sampleID,alib,alane))	
		elif (rawopts) and (not cutadapter) and (not rmadapter):######only rawopts ######
			cmd.append('raw2clean_QC -i %s,%s -o ./ %s' % (rawfq1,rawfq2,rawopts))
			cmd.append('mv %s_1.clean.fq %s_%s_1.fq'%(alib.split("_")[0]+"_"+alib.split("_")[2],self.sampleID,alib))
			#cmd.append('mv %s_%s_1.clean.fq %s_%s_1.fq'%(self.sampleID,alib,self.sampleID,alib))
			cmd.append('gzip -f %s_%s_1.fq'%(self.sampleID,alib))
			cmd.append('mv %s_2.clean.fq %s_%s_2.fq'%(alib.split("_")[0]+"_"+alib.split("_")[2],self.sampleID,alib))
			#cmd.append('mv %s_%s_2.clean.fq %s_%s_2.fq'%(self.sampleID,alib,self.sampleID,alib))
			cmd.append('gzip -f %s_%s_2.fq'%(self.sampleID,alib))

		elif (cutadapter) and (not rmadapter) and (not rawopts):######only cutadapter ######
			step0 = 'rm -f %s_%s_1.fq.gz %s_%s_2.fq.gz %s_%s_1.adapter.list.gz %s_%s_2.adapter.list.gz\n' \
				% (self.sampleID,alib,self.sampleID,alib,self.sampleID,alib,self.sampleID,alib)
			step1 = '\\\n\t'.join(['cutadapt '+cutadapt_options,
				'-b AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
				'-B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT',
				'-o %s_%s_1.fq.gz -p %s_%s_2.fq.gz' %(self.sampleID,alib,self.sampleID,alib),
				rawfq1,rawfq2])
			cmd.append(step0 + step1)
			rawfq1 = os.path.join(self.rawdir,self.sampleID,'%s_%s_1.fq.gz'%(self.sampleID,alib))
			rawfq2 = os.path.join(self.rawdir,self.sampleID,'%s_%s_2.fq.gz'%(self.sampleID,alib))
		elif (rmadapter) and (not cutadapter) and (not rawopts):######only rmadapter ######
			cutadapt_options += ' --discard-trimmed'
			step0 = 'rm -f %s_%s_1.fq.gz %s_%s_2.fq.gz %s_%s_1.adapter.list.gz %s_%s_2.adapter.list.gz\n' \
				% (self.sampleID,alib,self.sampleID,alib,self.sampleID,alib,self.sampleID,alib)
			step1 = '\\\n\t'.join(['cutadapt '+cutadapt_options,
				'-b AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
				'-B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT',
				'-o %s_%s_1.fq.gz -p %s_%s_2.fq.gz' %(self.sampleID,alib,self.sampleID,alib),
				rawfq1,rawfq2])
			cmd.append(step0 + step1)
			rawfq1 = os.path.join(self.rawdir,self.sampleID,'%s_%s_1.fq.gz'%(self.sampleID,alib))
			rawfq2 = os.path.join(self.rawdir,self.sampleID,'%s_%s_2.fq.gz'%(self.sampleID,alib))
		elif rmadapter and rawopts and (not cutadapter):###### rmadapter + rawopts ######
			cutadapt_options += ' --discard-trimmed'
			step0 = 'rm -f %s_%s_1.fq.gz %s_%s_2.fq.gz %s_%s_1.adapter.list.gz %s_%s_2.adapter.list.gz\n' \
				% (self.sampleID,alib,self.sampleID,alib,self.sampleID,alib,self.sampleID,alib)
			step1 = '\\\n\t'.join(['cutadapt '+cutadapt_options,
				'-b AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
				'-B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT',
				'-o %s_%s_1.fq.gz -p %s_%s_2.fq.gz' %(self.sampleID,alib,self.sampleID,alib),
				rawfq1,rawfq2])
			cmd.append(step0 + step1)
			rawfq1 = os.path.join(self.rawdir,self.sampleID,'%s_%s_1.fq.gz'%(self.sampleID,alib))
			rawfq2 = os.path.join(self.rawdir,self.sampleID,'%s_%s_2.fq.gz'%(self.sampleID,alib))
			
			cmd.append('raw2clean_QC -i %s,%s -o ./ %s' % (rawfq1,rawfq2,rawopts))
			cmd.append('mv %s_%s_1.clean.fq %s_%s_1.fq'%(self.sampleID,alib,self.sampleID,alib))
			cmd.append('gzip -f %s_%s_1.fq'%(self.sampleID,alib))
			cmd.append('mv %s_%s_2.clean.fq %s_%s_2.fq'%(self.sampleID,alib,self.sampleID,alib))
			cmd.append('gzip -f %s_%s_2.fq'%(self.sampleID,alib))
		elif cutadapter and rawopts and (not rmadapter):###### cutadapter + rawopts ######
			step0 = 'rm -f %s_%s_1.fq.gz %s_%s_2.fq.gz %s_%s_1.adapter.list.gz %s_%s_2.adapter.list.gz\n' \
				% (self.sampleID,alib,self.sampleID,alib,self.sampleID,alib,self.sampleID,alib)
			step1 = '\\\n\t'.join(['cutadapt '+cutadapt_options,
				'-b AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC',
				'-B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT',
				'-o %s_%s_1.fq.gz -p %s_%s_2.fq.gz' %(self.sampleID,alib,self.sampleID,alib),
				rawfq1,rawfq2])
			cmd.append(step0 + step1)
			rawfq1 = os.path.join(self.rawdir,self.sampleID,'%s_%s_1.fq.gz'%(self.sampleID,alib))
			rawfq2 = os.path.join(self.rawdir,self.sampleID,'%s_%s_2.fq.gz'%(self.sampleID,alib))
			
			cmd.append('raw2clean_QC -i %s,%s -o ./ %s' % (rawfq1,rawfq2,rawopts))
			cmd.append('mv %s_%s_1.clean.fq %s_%s_1.fq'%(self.sampleID,alib,self.sampleID,alib))
			cmd.append('gzip -f %s_%s_1.fq'%(self.sampleID,alib))
			cmd.append('mv %s_%s_2.clean.fq %s_%s_2.fq'%(self.sampleID,alib,self.sampleID,alib))
			cmd.append('gzip -f %s_%s_2.fq'%(self.sampleID,alib))
		if rmdup:
			step2 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['QC_report'],'duplication_rm','removedup_from_fastq.py'),
				'-l %s' % rawfq1,
				'-r %s' % rawfq2,
				'-o %s' % os.path.join(self.rawdir,self.sampleID),
				'-p %s_%s' % (self.sampleID,alib),
				'-d %d' % duplevel])
			cmd.append(step2)

		if fqcheck:
			step1 = '\\\n\t'.join(['%s' % os.path.join(self.softwares['QC_report'],'fqcheck_adapter'),
				'-a %s' % os.path.join(self.softwares['QC_report'],'p7_adapter.fa'),
				'-r %s_%s_1.fq.gz' % (self.sampleID,alib),
				'-l %s_%s_1.adapter.list.gz' % (self.sampleID,alib),
				'-s %s_%s_1.adapter.stat' % (self.sampleID,alib),
				'-c %s_%s_1.adapter.fqcheck' % (self.sampleID,alib)])
			step2 = '\\\n\t'.join(['%s' % os.path.join(self.softwares['QC_report'],'fqcheck_adapter'),
				'-a %s' % os.path.join(self.softwares['QC_report'],'p5_adapter.fa'),
				'-r %s_%s_2.fq.gz' % (self.sampleID,alib),
				'-l %s_%s_2.adapter.list.gz' % (self.sampleID,alib),
				'-s %s_%s_2.adapter.stat' % (self.sampleID,alib),
				'-c %s_%s_2.adapter.fqcheck' % (self.sampleID,alib)])
			cmd.append(step1)
			cmd.append(step2)
		return ' && \\\n'.join(cmd)+'\n'

	def md5_raw(self,alib,alane):
		order = 'order md5_raw_%s_%s_%s after ln_%s_%s_%s' % (self.sampleID,alib,alane,self.sampleID,alib,alane)
		step1 = 'cd %s' % (os.path.join(self.rawdir,self.sampleID))
		step2 = 'md5sum %s_%s_%s_1.fq.gz > %s_%s_%s.fq.gz.MD5.txt' % (self.sampleID,alib,alane,self.sampleID,alib,alane)
		step3 = 'md5sum %s_%s_%s_2.fq.gz >> %s_%s_%s.fq.gz.MD5.txt' % (self.sampleID,alib,alane,self.sampleID,alib,alane)
		return ' && \\\n'.join([step1,step2,step3])+'\n',order

	def qc(self,alib,alane,opts='',fqcheck=False,rmadapter=False,singlecell=False,cutadapter=False):
		order = 'order qc_%s_%s_%s after ln_%s_%s_%s' % (self.sampleID,alib,alane,self.sampleID,alib,alane)
		rawfq1 = self.fastqs[alib][alane][0]
		rawfq2 = self.fastqs[alib][alane][1]
		adapt1 = rawfq1[:-5]+'adapter.list.gz'
		adapt2 = rawfq2[:-5]+'adapter.list.gz'
		fq1 = os.path.join(self.rawdir,self.sampleID,'%s_%s_%s_1.fq.gz'%(self.sampleID,alib,alane))
		fq2 = os.path.join(self.rawdir,self.sampleID,'%s_%s_%s_2.fq.gz'%(self.sampleID,alib,alane))

		cmds = []
		step1 = 'cd %s' % self.qcDir
		cmds.append(step1)
		## single cell
		if singlecell:
			step2 = '\\\n\t'.join(['cutadapt -n 2 -m 60 -O 10 ',
				'-b GTGAGTGATGGTTGAGGTAGTGTGGAG -B GTGAGTGATGGTTGAGGTAGTGTGGAG',
				'-o %s' % os.path.join(self.qcDir,'%s_%s_%s_1.fq.gz'%(self.sampleID,alib,alane)),
				'-p %s' % os.path.join(self.qcDir,'%s_%s_%s_2.fq.gz'%(self.sampleID,alib,alane)),
				'%s %s' % (fq1,fq2)])
			cmds.append(step2)
			fq1 = os.path.join(self.qcDir,'%s_%s_%s_1.fq.gz'%(self.sampleID,alib,alane))
			fq2 = os.path.join(self.qcDir,'%s_%s_%s_2.fq.gz'%(self.sampleID,alib,alane))
			
			##fqcheck##
			step3 = '\\\n\t'.join(['%s' % os.path.join(self.softwares['QC_report'],'fqcheck_adapter'),
				'-a %s' % os.path.join(self.softwares['QC_report'],'p7_adapter.fa'),
				'-r %s_%s_%s_1.fq.gz' % (self.sampleID,alib,alane),
				'-l %s_%s_%s_1.adapter.list.gz' % (self.sampleID,alib,alane),
				'-s %s_%s_%s_1.adapter.stat' % (self.sampleID,alib,alane),
				'-c %s_%s_%s_1.adapter.fqcheck' % (self.sampleID,alib,alane)])
			step4 = '\\\n\t'.join(['%s' % os.path.join(self.softwares['QC_report'],'fqcheck_adapter'),
				'-a %s' % os.path.join(self.softwares['QC_report'],'p5_adapter.fa'),
				'-r %s_%s_%s_2.fq.gz' % (self.sampleID,alib,alane),
				'-l %s_%s_%s_2.adapter.list.gz' % (self.sampleID,alib,alane),
				'-s %s_%s_%s_2.adapter.stat' % (self.sampleID,alib,alane),
				'-c %s_%s_%s_2.adapter.fqcheck' % (self.sampleID,alib,alane)])
			cmds.append(step3)
			cmds.append(step4)
		## qc and adapter
		step5 = 'raw2clean_QC -i %s,%s -o ./ %s ' % (fq1, fq2, opts)
		if fqcheck:
			step5 += '-a %s,%s' % (fq1[:-5]+'adapter.list.gz',fq2[:-5]+'adapter.list.gz')
		else:
			if os.path.exists(adapt1) and os.path.exists(adapt2) and (not rmadapter) and (not cutadapter):
				step5 += '-a %s,%s' % (fq1[:-5]+'adapter.list.gz',fq2[:-5]+'adapter.list.gz')
		cmds.append(step5)
		step6 = 'echo qc_%s_%s_%s done!' % (self.sampleID,alib,alane)
		cmds.append(step6)
		return ' && \\\n'.join(cmds)+'\n',order

	def qc_summary(self,lib2lane):
		order = ['order qc_summary_%s after qc_%s_%s_%s' % (self.sampleID,self.sampleID,lib,lane) \
			for lib in lib2lane for lane in lib2lane[lib]]
		statlist = ['%s_%s_%s.stat'%(self.sampleID,lib,lane) for lib in lib2lane for lane in lib2lane[lib]]
		step1 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'summary_qcstat'),
			','.join(statlist),
			'> %s' % os.path.join(self.qcDir,self.sampleID+'.qcstat')])
		return step1,order

	def md5_clean(self,alib,alane):
		order = 'order md5_clean_%s_%s_%s after qc_%s_%s_%s' % (self.sampleID,alib,alane,self.sampleID,alib,alane)
		step1 = 'cd %s' % self.qcDir
		step2 = 'md5sum %s_%s_%s_1.clean.fq.gz > %s_%s_%s_1.clean.fq.gz.MD5.txt' % (self.sampleID,alib,alane,self.sampleID,alib,alane)
		step3 = 'md5sum %s_%s_%s_2.clean.fq.gz > %s_%s_%s_2.clean.fq.gz.MD5.txt' % (self.sampleID,alib,alane,self.sampleID,alib,alane)
		return ' && \\\n'.join([step1,step2,step3])+'\n',order
		

	def pollution(self,alib,alane):
		order = 'order pollution_%s after ln_%s_%s_%s' % (self.sampleID,self.sampleID,alib,alane)
		step0 = 'cd %s\nif [ ! -f "%s.pollution.xls" ];then' % (self.qcDir, self.sampleID)
		step1 = '\\\n\t'.join(['%s' % os.path.join(self.softwares['QC_report'],'pollution','generate_contamination_QC_v2.pl'),
			'%s/%s_%s_%s_1.fq.gz' %(os.path.join(self.rawdir,self.sampleID),self.sampleID,alib,alane),'1e-8 ',
			'>%s/pollution.sh' % self.qcDir])
		step2 = 'sh %s/pollution.sh' % self.qcDir
		step3 = 'mv Species_Kingdoms_distribution_result.xls %s.pollution.xls\nfi\n' % (self.sampleID)
		return step0+'\n'+' && \\\n'.join([step1,step2,step3])+'\n',order
