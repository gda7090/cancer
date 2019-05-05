import os,sys
import string

class QC:

	def __init__(self,sampleID,rawDir,qcDir,fastqs,softwares):
		self.sampleID = sampleID
		self.rawdir = rawDir
		self.qcdir = qcDir
		self.fastqs = fastqs
		self.softwares = softwares

	def ln(self,alib,rawopts='',fqcheck=False,rmadapter=False,cutadapter=False,rmdup=False,duplevel=15):
		rawfq1 = self.fastqs[alib][0]
		rawfq2 = self.fastqs[alib][1]
		adapt1 = rawfq1[:-5]+'adapter.list.gz'
		adapt2 = rawfq2[:-5]+'adapter.list.gz'
		cmd = []
		cmd.append('cd %s'%os.path.join(self.rawdir,self.sampleID))

		cutadapt_options = '-n 1 -m 100 -O 10 -e 0.2'
		
		if (not rawopts) and (not cutadapter) and (not rmadapter):######no rm,cut,rawopts####
			cmd.append('ln -sf %s %s_%s_1.fq.gz'%(rawfq1,self.sampleID,alib))
			cmd.append('ln -sf %s %s_%s_2.fq.gz'%(rawfq2,self.sampleID,alib))
			if os.path.exists(adapt1) and os.path.exists(adapt2) and (not fqcheck):
				cmd.append('ln -sf %s %s_%s_1.adapter.list.gz'%(adapt1,self.sampleID,alib))
				cmd.append('ln -sf %s %s_%s_2.adapter.list.gz'%(adapt2,self.sampleID,alib))	
		elif (rawopts) and (not cutadapter) and (not rmadapter):######only rawopts ######
			cmd.append('raw2clean_QC -i %s,%s -o ./ %s' % (rawfq1,rawfq2,rawopts))
			cmd.append('mv %s_1.clean.fq %s_%s_1.fq'%(alib.split("_")[0]+"_"+alib.split("_")[2],self.sampleID,alib))
			#cmd.append('mv %s_%s_1.clean.fq %s_%s_1.fq'%(self.sampleID,alib,self.sampleID,alib))
			cmd.append('gzip -f %s_%s_1.fq'%(self.sampleID,alib))
			cmd.append('mv %s_2.clean.fq %s_%s_2.fq'%(alib.split("_")[0]+"_"+alib.split("_")[2],self.sampleID,alib))
			#cmd.append('mv %s_%s_2.clean.fq %s_%s_2.fq'%(self.sampleID,alib,self.sampleID,alib))
			cmd.append('gzip -f %s_%s_2.fq'%(self.sampleID,alib))

		elif (cutadapter) and (not rmadapter) and (not rawopts):######only cutadapter ######
			step0 = 'rm %s_%s_1.fq.gz %s_%s_2.fq.gz %s_%s_1.adapter.list.gz %s_%s_2.adapter.list.gz\n' \
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
			step0 = 'rm %s_%s_1.fq.gz %s_%s_2.fq.gz %s_%s_1.adapter.list.gz %s_%s_2.adapter.list.gz\n' \
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
			step0 = 'rm %s_%s_1.fq.gz %s_%s_2.fq.gz %s_%s_1.adapter.list.gz %s_%s_2.adapter.list.gz\n' \
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
			step0 = 'rm %s_%s_1.fq.gz %s_%s_2.fq.gz %s_%s_1.adapter.list.gz %s_%s_2.adapter.list.gz\n' \
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
			step2 = '\\\n\t'.join(['python /PROJ/HUMAN/share/Cancer/QC_report/duplication_rm/removedup_from_fastq.py',
				'-l %s' % rawfq1,
				'-r %s' % rawfq2,
				'-o %s' % os.path.join(self.rawdir,self.sampleID),
				'-p %s_%s' % (self.sampleID,alib),
				'-d %d' % duplevel])
			cmd.append(step2)

		if fqcheck:
			step1 = '\\\n\t'.join(['/PROJ/HUMAN/share/Cancer_new/QC_report/fqcheck_adapter',
				'-a /PROJ/HUMAN/share/Cancer_new/QC_report/p7_adapter.fa',
				'-r %s_%s_1.fq.gz' % (self.sampleID,alib),
				'-l %s_%s_1.adapter.list.gz' % (self.sampleID,alib),
				'-s %s_%s_1.adapter.stat' % (self.sampleID,alib),
				'-c %s_%s_1.adapter.fqcheck' % (self.sampleID,alib)])
			step2 = '\\\n\t'.join(['/PROJ/HUMAN/share/Cancer_new/QC_report/fqcheck_adapter',
				'-a /PROJ/HUMAN/share/Cancer_new/QC_report/p5_adapter.fa',
				'-r %s_%s_2.fq.gz' % (self.sampleID,alib),
				'-l %s_%s_2.adapter.list.gz' % (self.sampleID,alib),
				'-s %s_%s_2.adapter.stat' % (self.sampleID,alib),
				'-c %s_%s_2.adapter.fqcheck' % (self.sampleID,alib)])
			cmd.append(step1)
			cmd.append(step2)
		return ' && \\\n'.join(cmd)+'\n'


	def md5(self,alib):
		order = 'order md5_%s_%s after ln_%s_%s' % (self.sampleID,alib,self.sampleID,alib)
		step1 = 'cd %s' % (os.path.join(self.rawdir,self.sampleID))
		step2 = 'md5sum %s_%s_1.fq.gz > %s_%s.fq.gz.MD5.txt' % (self.sampleID,alib,self.sampleID,alib)
		step3 = 'md5sum %s_%s_2.fq.gz >> %s_%s.fq.gz.MD5.txt' % (self.sampleID,alib,self.sampleID,alib)
		return ' && \\\n'.join([step1,step2,step3])+'\n',order

	
	def qc(self,alib,opts='',fqcheck=False,rmadapter=False,singlecell=False,cutadapter=False):
		order = 'order qc_%s_%s after ln_%s_%s' % (self.sampleID,alib,self.sampleID,alib)
		rawfq1 = self.fastqs[alib][0]
		rawfq2 = self.fastqs[alib][1]
		adapt1 = rawfq1[:-5]+'adapter.list.gz'
		adapt2 = rawfq2[:-5]+'adapter.list.gz'
		fq1 = os.path.join(self.rawdir,self.sampleID,'%s_%s_1.fq.gz'%(self.sampleID,alib))
		fq2 = os.path.join(self.rawdir,self.sampleID,'%s_%s_2.fq.gz'%(self.sampleID,alib))

		cmds = []
		step1 = 'cd %s' % (os.path.join(self.qcdir,self.sampleID))
		cmds.append(step1)
		## single cell
		if singlecell:
			step2 = '\\\n\t'.join(['cutadapt -n 2 -m 60 -O 10 ',
				'-b GTGAGTGATGGTTGAGGTAGTGTGGAG -B GTGAGTGATGGTTGAGGTAGTGTGGAG',
				'-o %s' % os.path.join(self.qcdir,self.sampleID,'%s_%s_1.fq.gz'%(self.sampleID,alib)),
				'-p %s' % os.path.join(self.qcdir,self.sampleID,'%s_%s_2.fq.gz'%(self.sampleID,alib)),
				'%s %s' % (fq1,fq2)])
			cmds.append(step2)
			fq1 = os.path.join(self.qcdir,self.sampleID,'%s_%s_1.fq.gz'%(self.sampleID,alib))
			fq2 = os.path.join(self.qcdir,self.sampleID,'%s_%s_2.fq.gz'%(self.sampleID,alib))
			
			##fqcheck##
			step3 = '\\\n\t'.join(['/PROJ/HUMAN/share/Cancer_new/QC_report/fqcheck_adapter',
                                '-a /PROJ/HUMAN/share/Cancer_new/QC_report/p7_adapter.fa',
                                '-r %s_%s_1.fq.gz' % (self.sampleID,alib),
                                '-l %s_%s_1.adapter.list.gz' % (self.sampleID,alib),
                                '-s %s_%s_1.adapter.stat' % (self.sampleID,alib),
                                '-c %s_%s_1.adapter.fqcheck' % (self.sampleID,alib)])
                        step4 = '\\\n\t'.join(['/PROJ/HUMAN/share/Cancer_new/QC_report/fqcheck_adapter',
                                '-a /PROJ/HUMAN/share/Cancer_new/QC_report/p5_adapter.fa',
                                '-r %s_%s_2.fq.gz' % (self.sampleID,alib),
                                '-l %s_%s_2.adapter.list.gz' % (self.sampleID,alib),
                                '-s %s_%s_2.adapter.stat' % (self.sampleID,alib),
                                '-c %s_%s_2.adapter.fqcheck' % (self.sampleID,alib)])
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
		step6 = 'echo qc_%s_%s done!' % (self.sampleID,alib)
		cmds.append(step6)
		return ' && \\\n'.join(cmds)+'\n',order

	def gzip_fq(self,alib):
		order = 'order gzip_fq_%s_%s after qc_%s_%s' % (self.sampleID,alib,self.sampleID,alib)
		step1 = 'cd %s' % (os.path.join(self.qcdir,self.sampleID))
		step2 = 'gzip -f %s_%s_1.clean.fq' % (self.sampleID,alib)
		step3 = 'gzip -f %s_%s_2.clean.fq' % (self.sampleID,alib)
		step4 = 'md5sum %s_%s_1.clean.fq.gz > %s_%s_1.clean.fq.gz.MD5.txt' % (self.sampleID,alib,self.sampleID,alib)
		step5 = 'md5sum %s_%s_2.clean.fq.gz > %s_%s_2.clean.fq.gz.MD5.txt' % (self.sampleID,alib,self.sampleID,alib)
		return ' && \\\n'.join([step1,step2,step3,step4,step5])+'\n',order
		

	def pollution(self,alib):
		order = 'order pollution_%s after ln_%s_%s' % (self.sampleID,self.sampleID,alib)
		step0 = 'cd %s\nif [ ! -f "%s.pollution.xls" ];then' % (os.path.join(self.qcdir,self.sampleID), self.sampleID)
		step1 = '\\\n\t'.join(['/PROJ/HUMAN/share/Cancer/script/generate_contamination_QC_v2.pl',
			'%s/%s_%s_1.fq.gz' %(os.path.join(self.rawdir,self.sampleID),self.sampleID,alib),'1e-8 ',
			'>%s/pollution.sh' % os.path.join(self.qcdir,self.sampleID)])
		step2 = 'sh %s/pollution.sh' % os.path.join(self.qcdir,self.sampleID)
		step3 = 'mv Species_Kingdoms_distribution_result.xls %s.pollution.xls\nfi\n' % (self.sampleID)
		return step0+'\n'+' && \\\n'.join([step1,step2,step3])+'\n',order
