import os

class MutationCalling:
		
	def __init__(self, pacientID, sampleID, alignDir, mutationDir, softwares, genome_info ,TR='',flank=True):
		self.pacientID = pacientID
		self.sampleID = sampleID
		self.alignDir = alignDir
		self.mutationDir = mutationDir
		self.gatkDIR = softwares['gatk']
		self.humandbDIR = genome_info['annovardb']
		self.refData = genome_info['fasta']
		self.TR = TR
		self.flank = flank
		self.advDir = softwares['advbin']
		self.genome_info = genome_info
		self.chroms = [chr for chrset in genome_info['chr_subs'] for chr in chrset]
		self.finalBam = os.path.join(alignDir,sampleID+'.final.bam')
		self.softwares = softwares
		self.annovar = softwares['annovar']
		self.snvdb = 'GeneName,refGene,genomicSuperDups,gff3,avsnp142,cosmic70,clinvar_20150330,gwasCatalog,1000g2014oct_Chinese,1000g2014oct_eas,1000g2014oct_all,esp6500siv2_all,exac03_ALL_EAS,ljb26_sift,ljb26_pp2hvar,ljb26_pp2hdiv,ljb26_mt,gerp++gt2'
		self.svdb = 'GeneName,Gencode,cpgIslandExt,cytoBand,dgvMerged'
	
	
	def samtoolsMpileup_sub(self, chrs, num):
		"""call variations by chromosomes.
		the argument [chrs] is 1,2,3,...,X,Y,M"""
		#order = 'order finalbam_%s before samtoolsMpileup_chr%s_%s' % (self.sampleID,chr,self.sampleID)
		order = 'order samtoolsMpileup_sub%s_%s after finalbam_%s' % (num,self.sampleID,self.sampleID)
		cmds = ['cd '+self.mutationDir]
		for chr in chrs:
			rawVcf = os.path.join(self.mutationDir,self.sampleID+'.samtools.%s.var.raw.vcf' % chr)
			if self.flank:
				cmd = '  \\\n\t'.join(['samtools mpileup',
					'-r '+chr,
					'-q 1 -C 50 -t DP,SP,DV -m 2 -F 0.002 -ugf',
					self.refData,
					self.finalBam,
					'|bcftools call ',
					'-vmO v -o '+rawVcf])
			else:
				cmd = '  \\\n\t'.join(['samtools mpileup',
					'-l '+self.TR,
					'-r '+chr,
					'-q 1 -C 50 -t DP,SP,DV -m 2 -F 0.002 -ugf',
					self.refData,
					self.finalBam,
					'|bcftools call ',
					'-vmO v -o '+rawVcf])
			cmds.append(cmd)
		return ' && \\\n'.join(cmds),order

	def cat_sub_all(self,subs):
		"""cat the variations in each chromosome together """
		rawVcf = os.path.join(self.mutationDir,self.sampleID+'.samtools.var.raw.vcf')
		vcflist = [os.path.join(self.mutationDir,'%s.samtools.%s.var.raw.vcf' \
			% (self.sampleID,chr)) for chr in self.chroms]

		order = ['order cat_sub_all_%s after samtoolsMpileup_sub%s_%s' % \
			(self.sampleID,sub,self.sampleID) for sub in subs]
		step1 = 'vcf-concat -p ' + '\\\n\t'.join(vcflist)+'\\\n\t>'+rawVcf
		step2 = 'rm -rf '+'\\\n\t'.join(vcflist)

		return ' && \\\n\t'.join([step1,step2]),order


	def filtersamtoolsCalling(self):
		#order = 'order cat_chr_all_%s before filtersamtoolsCalling_%s' % (self.sampleID,self.sampleID)
		order = 'order filtersamtoolsCalling_%s after cat_sub_all_%s' % (self.sampleID,self.sampleID)
		#order += 'order samtoolsMpileup_%s before filtersamtoolsCalling_%s' % (self.sampleID,self.sampleID)
		rawVcf = os.path.join(self.mutationDir,self.sampleID+'.samtools.var.raw.vcf')
		befltVcf = os.path.join(self.mutationDir,self.sampleID+'.samtools.var.beflt.vcf')
		fltVcf = os.path.join(self.mutationDir,self.sampleID+'.samtools.var.flt.vcf')
		snpVcf = os.path.join(self.mutationDir,self.sampleID+'.samtools.snp.vcf')
		indelVcf = os.path.join(self.mutationDir,self.sampleID+'.samtools.indel.vcf')
		step1 = '  \\\n\t'.join(['bcftools filter',
			'-s FLTER -i "%QUAL>30 && DP>4 && MQ>40" ',
			rawVcf,
			'>',
			befltVcf])
		step2 = 'awk -F "\\t" \'{if(/^#/){print}else{if($7=="PASS"){print}}}\' %s > %s ' % (befltVcf,fltVcf)
		step3 = "awk '/^#/ || !/INDEL;/' %s > %s" %(fltVcf,snpVcf)
		step4 = "awk '/^#/ || /INDEL;/' %s > %s" %(fltVcf,indelVcf)
		step5 = 'rm -f %s %s' % (fltVcf,befltVcf)
		return ' && \\\n'.join([step1,step2,step3,step4,step5]),order
	
	def gatk_calling(self):
		#order = 'order finalbam_%s before gatk_calling_%s' % (self.sampleID,self.sampleID)
		order = 'order gatk_calling_%s after finalbam_%s' % (self.sampleID,self.sampleID)
		step1 = ' \\\n'.join([self.softwares['java7']+' -Xmx6g -jar '+os.path.join(self.gatkDIR,'GenomeAnalysisTK.jar'),
			'-T HaplotypeCaller',
			'-R %s' % self.refData,
			'-I %s' % self.finalBam,
			'-L %s' % self.TR,
			'-nct 8',
			'--genotyping_mode DISCOVERY',
			'-stand_emit_conf 10',
			'-stand_call_conf 30',
			'-o %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.var.vcf')])
		if self.flank:
			step1 = ' \\\n'.join([self.softwares['java7']+' -Xmx6g -jar '+os.path.join(self.gatkDIR,'GenomeAnalysisTK.jar'),
				'-T HaplotypeCaller',
				'-R %s' % self.refData,
				'-I %s' % self.finalBam,
				'-nct 8',
				'--genotyping_mode DISCOVERY',
				'-stand_emit_conf 10',
				'-stand_call_conf 30',
				'-o %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.var.vcf')])
		
		return step1, order
	
	def gatk_variantion_filter(self):
		#order = 'order gatk_calling_%s before gatk_variantion_filter_%s' % (self.sampleID,self.sampleID)
		order = 'order gatk_variantion_filter_%s after gatk_calling_%s' % (self.sampleID,self.sampleID)
		#step1 = ' \\\n'.join(['java -Xmx6g -jar '+os.path.join(self.gatkDIR,'GenomeAnalysisTK.jar'),
			  #'-T VariantRecalibrator',
			  #'-R %s' % self.refData,
			  #'-input %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.var.vcf'),
			  #'-resource:hapmap,known=false,training=true,truth=true,prior=15.0 /PROJ/HUMAN/share/Cancer_pipeline_test/new_pip/hapmap_3.3.hg19.vcf',
			  ##'-resource:omni,known=false,training=true,truth=true,prior=12.0 omni.vcf',
			  ##'-resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G.vcf',
			  #'-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /PROJ/HUMAN/share/Cancer_pipeline_test/new_pip/00-All.vcf',
			  #'-an DP',
			  #'-an QD',
			  #'-an FS',
			  #'-an MQRankSum',
			  #'-an ReadPosRankSum',
			  #'-mode SNP',
			  #'-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0',
			  #'-recalFile %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.var.vcf.recalibrate_SNP.recal'),
			  #'-tranchesFile %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.var.vcf.recalibrate_SNP.tranches'),
			  #'-rscriptFile %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.var.vcf.recalibrate_SNP_plots.R')])
		#step2 = ' \\\n'.join(['java -Xmx6g -jar '+os.path.join(self.gatkDIR,'GenomeAnalysisTK.jar'),
			  #'-T ApplyRecalibration',
			  #'-R %s' % self.refData,
			  #'-input %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.var.vcf'),
			  #'-mode SNP',
			  #'--ts_filter_level 99.0',
			  #'-recalFile %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.var.vcf.recalibrate_SNP.recal'),
			  #'-tranchesFile %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.var.vcf.recalibrate_SNP.tranches'),
			  #'-o %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.var.recalibrated_snps_raw_indels.vcf')])
		#step3 = ' \\\n'.join(['java -Xmx6g -jar '+os.path.join(self.gatkDIR,'GenomeAnalysisTK.jar'),
			  #'-T VariantRecalibrator',
			  #'-R %s' % self.refData,
			  #'-input %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.var.recalibrated_snps_raw_indels.vcf'),
			  #'-resource:mills,known=true,training=true,truth=true,prior=12.0 /PROJ/HUMAN/share/Cancer_pipeline_test/new_pip/Mills_and_1000G_gold_standard.indels.b37.vcf',
			  #'-an DP',
			  #'-an FS',
			  #'-an MQRankSum',
			  #'-an ReadPosRankSum',
			  #'-mode INDEL',
			  #'-tranche 100.0',
			  #'--maxGaussians 4',
			  #'-recalFile %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.var.vcf.recalibrate_INDEL.recal'),
			  #'-tranchesFile %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.var.vcf.recalibrate_INDEL.tranches'),
			  #'-rscriptFile %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.var.vcf.recalibrate_INDEL_plots.R')])
		#step4 = ' \\\n'.join(['java -Xmx6g -jar '+os.path.join(self.gatkDIR,'GenomeAnalysisTK.jar'),
			  #'-T ApplyRecalibration',
			  #'-R %s' % self.refData,
			  #'-input %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.var.recalibrated_snps_raw_indels.vcf'),
			  #'-mode INDEL',
			  #'--ts_filter_level 99.0',
			  #'-recalFile %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.var.vcf.recalibrate_INDEL.recal'),
			  #'-tranchesFile %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.var.vcf.recalibrate_INDEL.tranches'),
			  #'-o %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.var.recalibrated_variants.vcf')])
		step5 = ' \\\n'.join([self.softwares['java7']+' -Xmx2g -jar '+os.path.join(self.gatkDIR,'GenomeAnalysisTK.jar'),
			'-T SelectVariants',
			'-R %s' % self.refData,
			'-V %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.var.vcf'),
			#'-L 20',
			'-selectType SNP',
			'-o %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.var.raw_snp.vcf')])
		step6 = ' \\\n'.join([self.softwares['java7']+' -Xmx2g -jar '+os.path.join(self.gatkDIR,'GenomeAnalysisTK.jar'),
			'-T VariantFiltration',
			'-R %s' % self.refData,
			'-V %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.var.raw_snp.vcf'),
			'--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0"',
			'--filterName "my_snp_filter"',
			'-o %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.snp.befil.vcf'),
			'2>/dev/null'])
		step7 = ' \\\n'.join([self.softwares['java7']+' -Xmx2g -jar '+os.path.join(self.gatkDIR,'GenomeAnalysisTK.jar'),
			'-T SelectVariants',
			'-R %s' % self.refData,
			'-V %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.var.vcf'),
			#'-L 20',
			'-selectType INDEL',
			'-o %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.var.raw_indel.vcf')])
		step8 = ' \\\n'.join([self.softwares['java7']+' -Xmx2g -jar '+os.path.join(self.gatkDIR,'GenomeAnalysisTK.jar'),
			'-T VariantFiltration',
			'-R %s' % self.refData,
			'-V %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.var.raw_indel.vcf'),
			'--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"',
			'--filterName "my_indel_filter"',
			'-o %s' % os.path.join(self.mutationDir,self.sampleID+'.GATK.indel.befil.vcf'),
			'2>/dev/null'])
		step9 = 'awk -F "\\t" \'{if(/^#/){print}else{if($7=="PASS" && $1!="hs37d5"){print}}}\' %s > %s ' % (os.path.join(self.mutationDir,self.sampleID+'.GATK.indel.befil.vcf'),os.path.join(self.mutationDir,self.sampleID+'.GATK.indel.vcf'))
		step10 = 'awk -F "\\t" \'{if(/^#/){print}else{if($7=="PASS" && $1!="hs37d5"){print}}}\' %s > %s ' % (os.path.join(self.mutationDir,self.sampleID+'.GATK.snp.befil.vcf'),os.path.join(self.mutationDir,self.sampleID+'.GATK.snp.vcf'))
		step11 = 'rm -f %s %s %s %s' % (os.path.join(self.mutationDir,self.sampleID+'.GATK.var.raw_snp.vcf'),os.path.join(self.mutationDir,self.sampleID+'.GATK.var.raw_indel.vcf'),os.path.join(self.mutationDir,self.sampleID+'.GATK.snp.befil.vcf'),os.path.join(self.mutationDir,self.sampleID+'.GATK.indel.befil.vcf'))
		step12 = 'rm -f %s.idx %s.idx %s.idx %s.idx' % (os.path.join(self.mutationDir,self.sampleID+'.GATK.var.raw_snp.vcf'),os.path.join(self.mutationDir,self.sampleID+'.GATK.var.raw_indel.vcf'),os.path.join(self.mutationDir,self.sampleID+'.GATK.snp.befil.vcf'),os.path.join(self.mutationDir,self.sampleID+'.GATK.indel.befil.vcf'))
		return ' && \\\n'.join([step5,step6,step7,step8,step9,step10,step11,step12]), order

	def annotatVcf(self,m):
		if m == 'samtools':
			#order = 'order filtersamtoolsCalling_%s before annotatVcf_%s' % (self.sampleID,self.sampleID)
			order = 'order annotatVcf_%s after filtersamtoolsCalling_%s' % (self.sampleID,self.sampleID)
		else:
			#order = 'order gatk_variantion_filter_%s before annotatVcf_%s' % (self.sampleID,self.sampleID)
			order = 'order annotatVcf_%s after gatk_variantion_filter_%s' % (self.sampleID,self.sampleID)
		snp = '  \\\n\t'.join([self.annovar,
			'-b SNP',
			'-r %s' % self.refData,
			'-u %s' % self.genome_info['annovarbuild'],
			'-d '+self.humandbDIR,
			'-z %s,caddgt10' % self.snvdb,
			os.path.join(self.mutationDir,self.sampleID+'.'+m+'.snp.vcf'),
			self.sampleID])
		indel = '  \\\n\t'.join([self.annovar,
			'-b INDEL',
			'-r %s' % self.refData,
			'-u %s' % self.genome_info['annovarbuild'],
			'-z %s,caddindel' % self.snvdb,
			'-d '+self.humandbDIR,
			os.path.join(self.mutationDir,self.sampleID+'.'+m+'.indel.vcf'),
			self.sampleID])
		gz_raw = 'bgzip -f %s && \\\ntabix -p vcf %s.gz' %(os.path.join(self.mutationDir,self.sampleID+'.'+m+'.var.raw.vcf'),os.path.join(self.mutationDir,self.sampleID+'.'+m+'.var.raw.vcf'))
		gz_snp = 'bgzip -f %s && \\\ntabix -p vcf %s.gz' %(os.path.join(self.mutationDir,self.sampleID+'.'+m+'.snp.vcf'),os.path.join(self.mutationDir,self.sampleID+'.'+m+'.snp.vcf'))
		gz_indel = 'bgzip -f %s && \\\ntabix -p vcf %s.gz' %(os.path.join(self.mutationDir,self.sampleID+'.'+m+'.indel.vcf'),os.path.join(self.mutationDir,self.sampleID+'.'+m+'.indel.vcf'))
		gz_xls1 = 'gzip -f %s %s' % (os.path.join(self.mutationDir,self.sampleID+'.'+m+'.snp.annovar.hg19_multianno.xls'),\
			os.path.join(self.mutationDir,self.sampleID+'.'+m+'.indel.annovar.hg19_multianno.xls'))
		gz_xls2 = 'gzip -f %s %s' % (os.path.join(self.mutationDir,self.sampleID+'.'+m+'.snp.pathway.xls'),\
			os.path.join(self.mutationDir,self.sampleID+'.'+m+'.indel.pathway.xls'))
		return ' && \\\n'.join([snp,indel,gz_raw,gz_snp,gz_indel,gz_xls1,gz_xls2])+'\n',order

	def merge_vcf(self,m,sample_list,projdir,sample_order):
		order = []
		comd = []
		snp_list = []
		indel_list = []
		mutdir = os.path.split(self.mutationDir)[0]
		merge_dir = os.path.join(mutdir,'merge_vcf')
		if not os.path.exists(merge_dir):
			os.mkdir(merge_dir)
		counts = os.path.join(merge_dir,'counts')
		if not os.path.exists(counts):
			comd.append('mkdir -p %s' % counts)
		pics = os.path.join(merge_dir,'pics')
		if not os.path.exists(pics):
			comd.append('mkdir -p %s' % pics)
		### samples in qc_list
		qc_list = []
		for line in open(os.path.join(projdir,'qc_list')):
			array = line.strip().split('\t')
			if array[2] not in qc_list:
				qc_list.append(array[2])
		sample_tmp = []
		for each in sample_order:
			if each in qc_list:
				sample_tmp.append(each)
		qc_list = sample_tmp
		##############################
		snp_list = [each+'\t'+os.path.join(mutdir,each+'.'+m+'/'+each+'.'+m+'.snp.reformated.vcf.gz.stat.txt') for each in qc_list]
		indel_list = [each+'\t'+os.path.join(mutdir,each+'.'+m+'/'+each+'.'+m+'.indel.reformated.vcf.gz.stat.txt') for each in qc_list]
		for k in sample_list:
			#order.append('order annotatVcf_%s before merge_vcf' % k)
			order.append('order merge_vcf after annotatVcf_%s' % k)
		comd.append('cd %s' % merge_dir)
		open(os.path.join(merge_dir,'snplist'),'w').write('\n'.join(snp_list)+'\n')
		open(os.path.join(merge_dir,'indellist'),'w').write('\n'.join(indel_list)+'\n')
		comd.append('python %s/summary_snp.py %s/snplist %s/counts'%(self.advDir,merge_dir,merge_dir))
		comd.append('python %s/summary_indel.py %s/indellist %s/counts'%(self.advDir,merge_dir,merge_dir))
		comd.append('Rscript %s/snp_indel_plot.R counts/  pics/' % self.advDir)
		comd.append('for i in `ls pics/*.png`;do sp=${i%.*};convert $sp.png $sp.JPEG;done')
		snpvcf_list = [each+'\t'+os.path.join(mutdir,each+'.'+m+'/'+each+'.'+m+'.snp.vcf') for each in qc_list]
		open(os.path.join(merge_dir,'snp.vcflist'),'w').write('\n'.join(snpvcf_list)+'\n')
		if len(qc_list) > 100:
			comd.append('python %s/hcluster_snp.py -i %s/snp.vcflist -o %s -r %s:1-20000000\n'%(self.advDir,merge_dir,merge_dir,self.chroms[0]))
		else:
			comd.append('python %s/hcluster_snp.py -i %s/snp.vcflist -o %s -r %s\n'%(self.advDir,merge_dir,merge_dir,self.chroms[0]))
		
		return ' && \\\n'.join(comd),order 
	

