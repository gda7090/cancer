import os,sys
import string

class RNAVariants:

	def __init__(self,sampleID,qcDir,alignDir,mutDir,fusionDir,softwares,databases):
		self.sampleID = sampleID
		self.qcDir = qcDir
		self.alignDir = alignDir
		self.mutDir = mutDir
		self.fusionDir = fusionDir
		self.bam = os.path.join(alignDir, self.sampleID+'.bam')
		self.softwares = softwares
		self.databases = databases

	def construct_fq(self,alib,alane):
		fq1_gz = os.path.join(self.qcDir,'_'.join([self.sampleID,alib,alane,'1.clean.fq.gz']))
		fq2_gz = os.path.join(self.qcDir,'_'.join([self.sampleID,alib,alane,'2.clean.fq.gz']))
		return fq1,fq2

	def process_bam(self):
		order = 'order process_bam_%s after mapping_rnaseq_%s' % (self.sampleID,self.sampleID)
		step1 = '\\\n\t'.join([self.softwares['java'],
			'-jar %s %s' % (os.path.join(self.softwares['picard']),'MarkDuplicates'),
			'VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true REMOVE_DUPLICATES=true MAX_FILE_HANDLES=1000',
			'TMP_DIR=%s' % self.mutDir,
			'I=%s' % self.bam,
			'O=%s' % os.path.join(self.mutDir,self.sampleID+'.nodup.bam')])
		step2 = '\\\n\t'.join([self.softwares['java'],
			'-jar %s' % self.softwares['GATK'],
			'-T SplitNCigarReads -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS',
			'-R %s' % self.databases['fasta'],
			'-I %s' % os.path.join(self.mutDir,self.sampleID+'.nodup.bam'),
			'-o %s' % os.path.join(self.mutDir,self.sampleID+'.splitN.bam')])
		step3 = 'rm %s' % os.path.join(self.mutDir,self.sampleID+'.nodup.bam')
		return ' && \\\n'.join([step1]),order

	def mutation_calling(self,soft='samtools'):
		order = 'order mutation_calling_%s after process_bam_%s' % (self.sampleID,self.sampleID)
		if soft == 'samtools':
			step1 = '\\\n\t'.join(['samtools mpileup',
				'-q 1 -C 50 -t DP,SP,DV -m 2 -F 0.002 -ugf',
				self.databases['fasta'], 
				os.path.join(self.mutDir,self.sampleID+'.splitN.bam'),
				'|bcftools call -vmO v',
				'-o %s' % os.path.join(self.mutDir,self.sampleID+'.raw.vcf')])
		elif soft == 'gatk':
			step1 = '\\\n\t'.join([self.softwares['java'],
				'-Xmx6g -jar %s' % self.softwares['GATK'],
				'-T HaplotypeCaller -stand_call_conf 20 -stand_emit_conf 20 -ploidy 2',
				'-rf BadCigar -dontUseSoftClippedBases --filter_reads_with_N_cigar',
				'-R %s' % self.databases['fasta'],
				'-I %s' % os.path.join(self.mutDir,self.sampleID+'.splitN.bam'),
				'-L %s' % self.databases['chrom_bed'],
				'-o %s' % os.path.join(self.mutDir,self.sampleID+'.raw.vcf')])
				
		return ' && \\\n'.join([step1]),order

	def filter_annotation(self,soft='samtools'):
		order = 'order filter_annotation_%s after mutation_calling_%s' % (self.sampleID,self.sampleID)
		if soft == 'samtools':
			step1 = '\\\n\t'.join(['bcftools filter',
				'-s FLTER -i "%QUAL>20 && DP>4 && MQ>30"',
				os.path.join(self.mutDir,self.sampleID+'.raw.vcf'),
				'> %s' % os.path.join(self.mutDir,self.sampleID+'.filter.vcf')])
			step2 = 'awk -F "\\t" \'{if(/^#/){print}else{if($7=="PASS" && !/INDEL;/){print}}}\' \\\n\t%s \\\n\t> %s ' % \
				(os.path.join(self.mutDir,self.sampleID+'.filter.vcf'),os.path.join(self.mutDir,self.sampleID+'.snp.vcf'))
			step3 = 'awk -F "\\t" \'{if(/^#/){print}else{if($7=="PASS" && /INDEL;/){print}}}\' \\\n\t%s \\\n\t> %s ' % \
				(os.path.join(self.mutDir,self.sampleID+'.filter.vcf'),os.path.join(self.mutDir,self.sampleID+'.indel.vcf'))
		else:
			step1 = '\\\n\t'.join([self.softwares['java'],
				'-Xmx1g -jar %s' % self.softwares['GATK'],
				'-T VariantFiltration -window 35 -cluster 3 ',
				'-filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0"',
				'-R %s' % self.databases['fasta'],
				'-V %s' % os.path.join(self.mutDir,self.sampleID+'.raw.vcf'),
				'-o %s' % os.path.join(self.mutDir,self.sampleID+'.filter.vcf')])
			step2 = '\\\n\t'.join([self.softwares['java'],
				'-Xmx1g -jar %s' % self.softwares['GATK'],
				'-T SelectVariants -selectType SNP',
				'-R %s' % self.databases['fasta'],
				'-V %s' % os.path.join(self.mutDir,self.sampleID+'.filter.vcf'),
				'-o %s' % os.path.join(self.mutDir,self.sampleID+'.snp.vcf')])
			step3 = '\\\n\t'.join([self.softwares['java'],
				'-Xmx1g -jar %s' % self.softwares['GATK'],
				'-T SelectVariants -selectType INDEL',
				'-R %s' % self.databases['fasta'],
				'-V %s' % os.path.join(self.mutDir,self.sampleID+'.filter.vcf'),
				'-o %s' % os.path.join(self.mutDir,self.sampleID+'.indel.vcf')])
		step4 = '\\\n\t'.join([self.softwares['annovar'],
			'-b SNP -u %s' % (self.databases['annovar_version']),
			'-r %s' % self.databases['fasta'],
			'-z %s' % self.self.databases['annovar_db'],
			'-d %s' % self.annovar_dbs,
			'%s %s' % (os.path.join(self.mutDir,self.sampleID+'.snp.vcf'),self.sampleID)])
		step5 = '\\\n\t'.join([self.softwares['annovar'],
			'-b INDEL -u %s' % (self.databases['annovar_version']),
			'-r %s' % self.databases['fasta'],
			'-z %s' % self.self.databases['annovar_db'],
			'-d %s' % self.annovar_dbs,
			'%s %s' % (os.path.join(self.mutDir,self.sampleID+'.indel.vcf'),self.sampleID)])
		step6 = 'rm %s' % os.path.join(self.mutDir,self.sampleID+'.filter.vcf')
		step7 = 'gzip -f %s' % os.path.join(self.mutDir,self.sampleID+'.raw.vcf')
		step8 = 'bgzip -f %s %s' % (os.path.join(self.mutDir,self.sampleID+'.snp.vcf'),os.path.join(self.mutDir,self.sampleID+'.indel.vcf'))
		step9 = 'tabix -p vcf %s.gz && \\\ntabix -p vcf %s.gz' % (os.path.join(self.mutDir,self.sampleID+'.snp.vcf'), \
			os.path.join(self.mutDir,self.sampleID+'.indel.vcf'))
		return ' && \\\n'.join([step1,step2,step3,step4,step5,step6,step7,step8,step9]),order

	def fusion_mapping(self,lib2lane):
		order = ['order fusion_mapping_%s after qc_%s_%s_%s' % (self.sampleID,self.sampleID,lib,lane) \
			for lib in lib2lane for lane in lib2lane[lib]]
		step1 = ''
		return ' && \\\n'.join([step1]),order

	def fusion_calling(self):
		order = 'order fusion_calling_%s after fusion_mapping_%s' % (self.sampleID,self.sampleID)
		step1 = ''
		return ' && \\\n'.join([step1]),order

