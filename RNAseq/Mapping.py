import os,sys
import string

class Mapping:

	def __init__(self,sampleID,qcDir,alignDir,softwares,databases):
		self.sampleID = sampleID
		self.qcDir = qcDir
		self.alignDir = alignDir
		self.bam = os.path.join(alignDir, self.sampleID+'.bam')
		self.chroms = [chr for chrset in databases['chr_subs'] for chr in chrset]
		self.softwares = softwares
		self.databases = databases

	def construct_fq(self,alib,alane):
		fq1_gz = os.path.join(self.qcDir,'_'.join([self.sampleID,alib,alane,'1.clean.fq.gz']))
		fq2_gz = os.path.join(self.qcDir,'_'.join([self.sampleID,alib,alane,'2.clean.fq.gz']))
		return fq1,fq2


	def mapping_rnaseq(self,lib2lane,soft='histat2'):
		order = ['order mapping_rnaseq_%s after qc_%s_%s_%s' % (self.sampleID,self.sampleID,lib,lane) \
			for lib in lib2lane for lane in lib2lane[lib]]
		clean_fqs = [self.construct_fq(lib,lane) for lib in lib2lane for lane in lib2lane[lib]]
		clean_fqs1 = [x[0] for x in clean_fqs]
		clean_fqs2 = [x[1] for x in clean_fqs]
		if soft == 'tophat2':
			cmd = '''cd %s
tophat-2.0.9/tophat2 \\
  -p 8 -N 2 --read-edit-dist 2 -a 8 -m 0 -g 20 -i 20 -I 1000000 -r 50 --mate-std-dev 20 --no-mixed \\
  --rg-id %s --rg-sample %s --rg-library %s --rg-platform illumina --rg-platform-unit PU \\
  -G %s --library-type fr-unstranded \\
  --tmp-dir %s/tophat2 \\
  -o %s/tophat2 \\
  %s \\
  %s \\
  %s \\
  2> %s/tophat2/%s.tophat2.log && \\
mv %s/tophat2/accepted_hits.bam %s && \\
samtools index %s
''' % (self.alignDir,self.sampleID,self.sampleID,lib2lane.keys()[0],self.databases['gtf'],self.alignDir,self.alignDir,self.databases['tophat2_index'],','.join(clean_fqs1),','.join(clean_fqs2),self.alignDir,self.sampleID,self.alignDir,self.bam,self.bam)
		elif soft == 'hisat2':
			cmd = '''cd %s
HISAT2/hisat2-2.0.5/hisat2 \\
  -p 4 --dta -t --phred33 --rg-id %s --rg "SM:%s" --rg "PL:illumina" --rg "LB:%s" \\
  -x %s \\
  -1 %s \\
  -2 %s \
  |samtools-1.3.1/bin/samtools sort -O BAM --threads 4 \\
  -o %s && \\
  samtools index %s
''' % (self.alignDir,self.sampleID,self.sampleID,lib2lane.keys()[0],self.databases['hisat2_index'],','.join(clean_fqs1),','.join(clean_fqs2),self.bam,self.bam)
		elif soft == 'star':
			rg = ['ID:%s_%s_%s SM:%s LB:%s PL:illumina'%(self.sampleID,lib,lane,self.sampleID,lib) \
				for lib in lib2lane for lane in lib2lane[lib]]
			cmd = '''cd %s
rm -rf %s/star
STAR2.5/bin/Linux_x86_64/STAR --readFilesCommand zcat \\
   --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --limitGenomeGenerateRAM 10000000000 \\
  --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterMismatchNmax 2 \\
  --runThreadN 8  --outSAMtype BAM SortedByCoordinate --outFileNamePrefix %s. \\
  --outSAMattrRGline %s \\
  --genomeDir %s \\
  --outTmpDir %s/star \\
  --readFilesIn %s %s && \\
mv %s/start/%s.Aligned.sortedByCoord.out.bam %s && \\
samtools index %s
''' % (self.alignDir,self.alignDir,self.sampleID,' , '.join(rg),self.databases['star_index'],self.alignDir,','.join(clean_fqs1),','.join(clean_fqs2),self.alignDir,self.sampleID,self.bam,self.bam)
		return cmd,order

	def cal_readcount(self,gtf=None):
		order = 'order cal_readcount_%s after mapping_rnaseq_%s' % (self.sampleID,self.sampleID)
		step1 = '\\\n\t'.join([self.softwares['featurecount'],
			'-T 8 -F GTF -t exon -g gene_id -s 0 -Q 10 -C -p',
			'-a %s' % self.databases['gtf'],
			'-o %s' % os.path.join(self.alignDir,self.sampleID+'.readcount'),
			self.bam])
		step2 = '\\\n\t'.join([os.path.join(self.softwares['bin'],'regionstat'),
			'--count %s' % os.path.join(self.alignDir,self.sampleID+'.readcount'),
			'--summary %s.summary' % os.path.join(self.alignDir,self.sampleID+'.readcount'),
			'--outfile %s' % os.path.join(self.alignDir,self.sampleID+'.mapregion.xls')])
		step3 = '\\\n\t'.join([os.path.join(self.softwares['bin'],'mapregion'),
			'--region %s' % os.path.join(self.alignDir,self.sampleID+'.mapregion.xls'),
			'--prefix %s' % os.path.join(self.alignDir,self.sampleID),
			'--sampleid %s' % self.sampleID])
		return ' && \\\n'.join([step1,step2,step3]),order

	def mapping_summary(self):
		order = 'order mapping_summary_%s after mapping_rnaseq_%s' % (self.sampleID,self.sampleID)
		step1 = '\\\n\t'.join([os.path.join(self.softwares['bin'],'bamstat'),
			os.path.join(self.qcDir,self.sampleID+'.qcstat'),
			os.path.join(self.alignDir,self.sampleID+'.mapregion.xls'),
			'>%s' % os.path.join(self.alignDir,self.sampleID+'.mapstat')])
		return ' && \\\n'.join([step1]),order

class Assembly:

	def __init__(self,sampleID,alignDir,assemDir,assemRootDir,softwares,databases):
		self.sampleID = sampleID
		self.alignDir = alignDir
		self.assemDir = assemDir
		self.assemRootDir = assemRootDir
		self.bam = os.path.join(alignDir, self.sampleID+'.bam')
		self.gtf = os.path.join(assemDir, self.sampleID+'.gtf')

	def assembly(self,soft='cufflinks'):
		order = 'order assembly_%s after mapping_rnaseq_%s' % (self.sampleID,self.sampleID)
		if soft == 'cufflinks':
			step1 = '\\\n\t'.join([self.softwares['cufflinks'],
				'-p 4 -u --library-type fr-firststrand',
				'-g %s' % self.databases['gtf'],
				'-o %s' % self.assemDir,
				self.bam])
			step2 = 'mv %s %s'%(os.path.join(self.assemDir,'transcripts.gtf'), self.gtf)
		elif soft == "stringtie":
			step1 = '\\\n\t'.join([self.softwares['stringtie'],
				'-p 4 -l %s' % self.sampleID,
				'-G %s' % self.databases['gtf'],
				'-o %s' % self.assemDir,
				self.bam])
			step2 = 'sed -i \'1,2d\' %s' % os.path.join(self.assemDir,self.sampleID+'.gtf')
		return ' && \\\n'.join([step1,step2]),order

	def merge_assembly(self,samples,soft='cufflinks',mode='gene'):
		order = ['order merge_assembly after assembly_%s' % each  for each in samples]
		open(os.path.join(self.assemRootDir,'merge_gtf','gtf.list'),'w').write( \
			'\n'.join([os.path.join(self.assemRootDir,each,each+'.gtf') for each in samples])+'\n')
		step1 = '\\\n\t'.join([self.softwares['cuffmerge'],
			'-o %s' % os.path.join(self.assemRootDir,'merge_gtf'),
			'-g %s' % self.databases['gtf'],
			'-s %s' % self.databases['fasta'],
			os.path.join(self.assemRootDir,'merge_gtf','gtf.list')])
		step2 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'compare_gtf.py'),
			self.databases['gtf'],
			os.path.join(self.assemRootDir,'merge_gtf','merged.gtf'),
			os.path.join(self.assemRootDir,'merge_gtf')])
		if mode == 'gene':
			step3 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'gff_merge_sort.py'),
				os.path.join(self.assemRootDir,'merge_gtf','novel_genes.gtf'),
				self.databases['gtf'],
				'> %s' % os.path.join(self.assemRootDir,'merge_gtf','known_novel.merge.gtf')])
		else:
			step3 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'gff_merge_sort.py'),
				os.path.join(self.assemRootDir,'merge_gtf','novel_transcripts.gtf'),
				self.databases['gtf'],
				'> %s' % os.path.join(self.assemRootDir,'merge_gtf','known_novel.merge.gtf')])
		step4 = 'rm %s' % ' '.join([os.path.join(self.assemRootDir,each,each+'.gtf') for each in samples])
		return ' && \\\n'.join([step1,step2,step3]),order

	def extract_fasta(self):
		order = 'order extract_fasta after merge_assembly'
		step1 = '\\\n\t'.join(['perl %s' % os.path.join(self.softwares['bin'],'extractFA'),
			os.path.join(self.assemRootDir,'merge_gtf','novel_transcripts.gtf'),
			self.databases['fasta'],
			os.path.join(self.assemRootDir,'lncRNA_filter','novel_transcripts.fasta')])
		return step1,order


	def lncRNA_cpc(self):
		order = 'order lncRNA_cpc after extract_fasta'
		step1 = '\\\n\t'.join(['perl %s' % os.path.join(self.softwares['bin'],'cpc'),
			'-seq %s' %  os.path.join(self.assemRootDir,'lncRNA_filter','novel_transcripts.fasta'),
			'-db eu -strand plus',
			'-outdir %s' % os.path.join(self.assemRootDir,'lncRNA_filter','CPC')])
		return step1,order

	def lncRNA_cnci(self):
		order = 'order lncRNA_cnci after extract_fasta'
		step1 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'cnci'),
			'-f %s' %  os.path.join(self.assemRootDir,'lncRNA_filter','novel_transcripts.fasta'),
			'-m ve -p 1',
			'-o %s' % os.path.join(self.assemRootDir,'lncRNA_filter','CNCI')])
		return step1,order

	def lncRNA_pfam(self):
		order = 'order lncRNA_pfam after extract_fasta'
		step1 = 'PATH=hmmer-3.1b1/bin:$PATH\nexport PERL5LIB=perl5/lib/perl5:perl5/lib/perl5/x86_64-linux:perl5/lib/perl5/x86_64-linux-thread-multi:Perl-5.18.2/lib/perl5/5.18.2/:Perl-5.18.2/lib/perl5/site_perl/5.18.2/x86_64-linux/:software/PfamScan'
		step2 = '\\\n\t'.join(['perl %s' % os.path.join(self.softwares['bin'],'seq2protein'),
			os.path.join(self.assemRootDir,'lncRNA_filter','novel_transcripts.fasta'),
			">%s" % os.path.join(self.assemRootDir,'lncRNA_filter','PFAM','novel_transcripts.protien.fasta')])
		step3 = '\\\n\t'.join(['perl %s' % os.path.join(self.softwares['bin'],'pfam'),
			'-fasta %s' %  os.path.join(self.assemRootDir,'lncRNA_filter','PFAM','novel_transcripts.fasta'),
			'-dir %s' % self.databases['pfam_hmm'],
			'-out %s' % os.path.join(self.assemRootDir,'lncRNA_filter','PFAM','Pfam_scan.out'),
			'-famB -cpu 2'])
		return step1,order

	def lncRNA_identification(self):
		order = 'order lncRNA_identification after merge_assembly'
		step1 = '\\\n\t'.join(['perl %s' % os.path.join(self.softwares['bin'],'lnc_venn'),
			'-gtf %s' % os.path.join(self.assemRootDir,'merge_gtf','novel_transcripts.gtf'),
			'-outdir %s' % os.path.join(self.assemRootDir,'lncRNA_filter'),
			'-CNCI %s' % os.path.join(self.assemRootDir,'lncRNA_filter','CNVI','CNCI.noncoding.id.txt'),
			'-CPC %s' % os.path.join(self.assemRootDir,'lncRNA_filter','CPC','CPC.noncoding.id.txt'),
			'-PFAM %s' % os.path.join(self.assemRootDir,'lncRNA_filter','PFAM','PFAM.noncoding.id.txt')])
		step2 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'get_new_gtf'),
			os.path.join(self.assemRootDir,'lncRNA_filter','noncoding.result.id'),
			os.path.join(self.assemRootDir,'lncRNA_filter','coding.result.id'),
			os.path.join(self.assemRootDir,'merge_gtf','novel_transcripts.gtf'),
			os.path.join(self.assemRootDir,'lncRNA_filter')])
		step3 = '\\\n\t'.join(['python %s' % os.path.join(self.softwares['bin'],'lnc_classify'),
			'--lnc_gtf %s' % os.path.join(self.assemRootDir,'lncRNA_filter','Novel_lncRNA.gtf'),
			'--out_dir %s' % os.path.join(self.assemRootDir,'lncRNA_filter')])
		return ' && \\\n'.join([step1,step2,step3]),order

	def lncRNA_signatures(self):
		order = 'order lncRNA_signatures after lncRNA_identification'
		step1 = ''
		return ' && \\\n'.join([step1]),order

class Quantification:

	def __init__(self,sampleID,qcDir,alignDir,expDir,softwares,databases):
		self.sampleID = sampleID
		self.qcDir = qcDir
		self.alignDir = alignDir
		self.assemDir = assemDir
		self.expDir = expDir
		self.bam = os.path.join(alignDir, self.sampleID+'.bam')
		self.gtf = os.path.join(assemDir, self.sampleID+'.gtf')

	def construct_fq(self,alib,alane):
		fq1_gz = os.path.join(self.qcDir,'_'.join([self.sampleID,alib,alane,'1.clean.fq.gz']))
		fq2_gz = os.path.join(self.qcDir,'_'.join([self.sampleID,alib,alane,'2.clean.fq.gz']))
		return fq1,fq2

	def prepare_rsem_index(self,gtf):
		order = 'order prepare_rsem_index after merge_assembly'
		step1 = ''
		return step1,order


	def rsem_quantification(self,assembly=False):
		order = ['order rsem_quantification_%s after qc_%s_%s_%s' % (self.sampleID,self.sampleID,lib,lane) \
			for lib in lib2lane for lane in lib2lane[lib]]
		rsem_index = self.databases['rsem_index']
		if assembly:
			order += ['order rsem_quantification_%s after prepare_rsem_index' % self.sampleID]
			rsem_index = os.path.join(self.expDir,'RSEM_INDEX')
		clean_fqs = [self.construct_fq(lib,lane,qc) for lib in lib2lane for lane in lib2lane[lib]]
		clean_fqs1 = [x[0] for x in clean_fqs]
		clean_fqs2 = [x[1] for x in clean_fqs]
		step1 = '\\\n\t'.join([self.softwares['rsem_exp'],
			'-p 8 --paired-end --star --estimate-rspd --forward-prob 0.5 --time --star-gzipped-read-file',
			'--star --star-path %s' % os.path.dirname(self.softwares['star']),
			'--temporary-folder %s' % self.expDir,
			','.join(clean_fqs1),
			','.join(clean_fqs2),
			rsem_index,
			os.path.join(self.expDir,self.sampleID)])
		step2 = 'rm %s' % os.path.join(self.expDir,self.sampleID+'transcripts.bam')
		return ' && \\\n'.join([step1,step2]),order

	def htseq_quantification(self,gtf=None):
		order = 'order htseq_quantification_%s after mapping_rnaseq_%s' % (self.sampleID,self.sampleID)
		gtf_quant = self.databases['gtf']
		if gtf:
			gtf_quant = gtf
		step1 = '\\\n\t'.join([self.softwares['featurecount'],
			'-T 8 -F GTF -t exon -g gene_id -s 0 -Q 10 -C -p',
			'-a %s' % gtf_quant,
			'-o %s' % os.path.join(self.expDir,self.sampleID+'.readcount'),
			self.bam])
		return step1,order

	def cuffquant_quantification(self,gtf=None):
		order = 'order cuffquant_quantification after mapping_rnaseq_%s' %(self.sampleID)
		gtf_quant = self.databases['gtf']
		if gtf:
			gtf_quant = gtf
		step1 = '\\\n\t'.join([self.softwares['cuffquant'],
			'--library-type fr-firststrand -p 2 --max-bundle-frags 1165754',
			'-o %s' % self.expDir,
			gtf_quant, self.bam])
		return step1,order

	def stringtie_quantification(self,gtf=None):
		order = 'order stringtie_quantification after mapping_rnaseq_%s' %(self.sampleID)
		gtf_quant = self.databases['gtf']
		if gtf:
			gtf_quant = gtf
		step1 = ''
		return step1,order

	def quantification(self,gtf=None,assembly=False,soft='htseq'):
		gtf_quant = self.databases['gtf']
		if gtf:
			gtf_quant = gtf
		if soft == 'htseq':
			script,order = self.htseq_quantification(gtf)
		elif soft == 'rsem':
			script,order = self.rsem_quantification(assembly)
		elif soft == 'cuffquant':
			script,order = self.cuffquant_quantification(gtf)
		elif soft == 'stringtie':
			script,order = self.stringtie_quantification(gtf)
		return script,order

	def quantification_summary(self,samples,soft='rsem'):
		order = ['order quantification_summary after %s_quantification_%s'%(soft,each) \
			for each in samples]
		step1 = ''
		return step1,order

		
