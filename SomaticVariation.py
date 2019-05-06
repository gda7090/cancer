import os
	
class Somatic:
	def __init__(self, pacientID, sampleID, Nname, projDir, genome_info, somaticDir, softwares, TR,seqstrag,tFreq=0.1, nFreq=0.05, flank=True):
		self.pacientID = pacientID
		self.sampleID = sampleID
		self.Nname = Nname
		mapDir = os.path.join(projDir,'Mapping')
		self.mapDir = mapDir
		self.statDir = os.path.join(projDir,'Alnstat')
		self.Tbam = os.path.join(os.path.join(mapDir,sampleID),sampleID+'.final.bam')
		self.Nbam = os.path.join(os.path.join(mapDir,Nname),Nname+'.final.bam')
		self.Tsbam = os.path.join(os.path.join(mapDir,sampleID),sampleID+'.split.bam')
		self.Nsbam = os.path.join(os.path.join(mapDir,Nname),Nname+'.split.bam')
		self.Tdbam = os.path.join(os.path.join(mapDir,sampleID),sampleID+'.discord.bam')
		self.Ndbam = os.path.join(os.path.join(mapDir,Nname),Nname+'.discord.bam')
		self.Nmpilegz = os.path.join(os.path.join(mapDir,Nname),Nname+'.final.mpileup.gz')
		self.Tmpilegz = os.path.join(os.path.join(mapDir,sampleID),sampleID+'.final.mpileup.gz')
		self.TR = TR
		self.tFreq = str(tFreq)
		self.nFreq = str(nFreq)
		self.refData = genome_info['fasta']
		self.refData2bit = genome_info['fa2bit']
		self.genome_info = genome_info
		self.chroms = [chr for chrset in genome_info['chr_subs'] for chr in chrset]
		self.somaticDir = somaticDir
		self.somaticsnvDir = os.path.join(somaticDir,'varScan')
		self.somaticsnvDirsam = os.path.join(somaticDir,'samtools')
		self.somaticsvDirc = os.path.join(somaticDir,'crest')
		self.somaticsvDirl = os.path.join(somaticDir,'lumpy')
		self.somaticsvDirb = os.path.join(somaticDir,'breakdancer')
		self.somaticcnvDirf = os.path.join(somaticDir,'freec')
		self.somaticcnvDirv = os.path.join(somaticDir,'varScanCNV')
		self.somaticcnvDire = os.path.join(somaticDir,'ExomeCNV')
		self.somaticsnvDirm = os.path.join(somaticDir,'muTect')
		self.somaticsnvDirs = os.path.join(somaticDir,'Strelka')
		self.varscanDir = softwares['varscan']
		self.breakdancerDir = softwares['breakdancer']
		self.gatkDir = softwares['gatk']
		self.muTectDir = softwares['mutect']
		self.strelkaDir = softwares['strelka']
		self.exomecnvDir = softwares['exomecnv']
		self.freecDir = softwares['freec']
		self.speedseqDir = softwares['speedseq']
		self.humandbDir = genome_info['annovardb']
		self.advDir = softwares['advbin']
		self.crestDir = softwares['crest']
		self.seqstrag = seqstrag
		self.platform = {'WGS':'HiseqX'}.get(seqstrag,'Hiseq')
		self.flank = flank
		self.softwares = softwares
		self.annovar = softwares['annovar']
		self.somatic_snvdb = 'GeneName,refGene,Gencode,wgRna,targetScanS,tfbsConsSites,genomicSuperDups,gff3,avsnp142,cosmic70,clinvar_20150330,gwasCatalog,1000g2014oct_eas,1000g2014oct_all,esp6500siv2_all,exac03_ALL_EAS,ljb26_sift,ljb26_pp2hvar,ljb26_pp2hdiv,ljb26_mt,gerp++gt2'
		self.somatic_svdb = 'GeneName,refGene,Gencode,cpgIslandExt,cytoBand,dgvMerged'
		self.ffpe = softwares['FFPE']	
	def somaticsamtoolsMpileup(self):
		#order = ['order finalbam_%s before somaticsamtoolsMpileup_%s' % (self.Nname,self.sampleID),
		#		 'order finalbam_%s before somaticsamtoolsMpileup_%s' % (self.sampleID,self.sampleID)]
		order = ['order somaticsamtoolsMpileup_%s after finalbam_%s' % (self.sampleID,self.Nname),
			 'order somaticsamtoolsMpileup_%s after finalbam_%s' % (self.sampleID,self.sampleID)]
		if self.TR and self.flank:
			return '  \\\n\t'.join(['samtools mpileup -q 1 -t DP,DV,SP',
				'-l '+self.TR,
				'-uvf '+self.refData,
				self.Nbam,
				self.Tbam,
				'> '+os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.mpileup')]), order
		else:
			return '  \\\n\t'.join(['samtools mpileup -q 1 -t DP,DV,SP',
				'-uvf '+self.refData,
				self.Nbam,
				self.Tbam,
				'> '+os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.mpileup')]), order
	def somaticsamtoolscalling(self):
		#order = 'order somaticsamtoolsMpileup_%s before somaticsamtoolscalling_%s' % (self.sampleID,self.sampleID)
		order = 'order somaticsamtoolscalling_%s after somaticsamtoolsMpileup_%s' % (self.sampleID,self.sampleID)
		step1 = 'bcftools call -vc %s > %s' % (os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.mpileup'),
			os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.mpileup.vcf'))
		step2 = 'python Cancer_pipeline/v1.5/filter.after_adjust.py %s %s' % (os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.mpileup.vcf'),
		#step2 = 'filter_somatic_snv.pl -normalDP 8 -tumorDP 8 %s > %s' % (os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.mpileup.vcf'),
			os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.vcf'))
		step3 = 'awk \'/^#/ || /INDEL;/\' %s > %s' % (os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.vcf'),
			os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.indel.vcf'))
		step4 = 'awk \'/^#/ || !/INDEL;/\' %s > %s' % (os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.vcf'),
			os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.snv.beflt.vcf'))
		step5 = 'awk \'!/^#/ {print $1,$2,$2}\' %s > %s' % (os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.snv.beflt.vcf'),
			os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.snv.beflt.vcf.list'))
		step6 = 'bam-readcount %s -q 1 -b 13 -f %s -l %s > %s' % (self.Tbam,
			self.refData,
			os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.snv.beflt.vcf.list'),
			os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.snv.beflt.vcf.readcount'))
		step7 = 'fpfilter_vcf.pl %s %s --output-basename %s' % (os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.snv.beflt.vcf'),
			os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.snv.beflt.vcf.readcount'),
			os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.snv.beflt.vcf.fpfilter'))
		step8 = 'awk \'/^#/\' %s > %s' % (os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.snv.beflt.vcf'), os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.snv.vcf'))
		step9 = 'cat %s >> %s' % (os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.snv.beflt.vcf.fpfilter.pass'), os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.snv.vcf'))
		step10 = 'rm -f %s ' % '\\\n\t'.join([os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.mpileup'),os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.vcf'),os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.snv.beflt.vcf'),os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.snv.beflt.vcf.list'),os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.snv.beflt.vcf.readcount'),os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.snv.beflt.vcf.fpfilter')])
		return ' && \\\n'.join([step1,step2,step3,step4,step5,step6,step7,step8,step9,step10]), order

	def somaticsamtoolsannovar(self):
		#order = 'order somaticsamtoolscalling_%s before somaticsamtoolsannovar_%s' %(self.sampleID,self.sampleID)
		order = 'order somaticsamtoolsannovar_%s after somaticsamtoolscalling_%s' %(self.sampleID,self.sampleID)
		somatic_snv_db1 = self.somatic_snvdb + ',caddgt10'
		somatic_snv_db2 = self.somatic_snvdb + ',caddindel'
		if self.genome_info['annovarbuild'] == 'mm10':
			somatic_snv_db1 = somatic_snv_db2 ='GeneName,refGene,genomicSuperDups,gff3,snp142chr'
		step1 = '  \\\n\t'.join([self.annovar,
			'-a somatic',
			'-r %s' % self.refData,
			'-u %s' % self.genome_info['annovarbuild'],
			'-b SNP',
			'-z %s' % somatic_snv_db1,
			'-d %s' % self.humandbDir,
			os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.snv.vcf'), 
			'%s,%s' % (self.Nname,self.sampleID)])
		step2 = '  \\\n\t'.join([self.annovar,
			'-a somatic',
			'-r %s' % self.refData,
			'-u %s' % self.genome_info['annovarbuild'],
			'-b INDEL',
			'-z %s' % somatic_snv_db2,
			'-d %s' % self.humandbDir,
			os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.indel.vcf'),
			'%s,%s' % (self.Nname,self.sampleID)])
		step3 = ' \\\n\t'.join(['python %s' % os.path.join(self.advDir,'vcf2maf.py'),
			'-v %s' % os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.snv.reformated.vcf.gz'),
			'-o %s' % os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.snv.maf'),
			'-t %s -n %s -s %s -p %s -m samtools -f avsnp142 -x caddgt10' % (self.sampleID,self.Nname,self.seqstrag,self.platform)])
		step4 = ' \\\n\t'.join(['python %s' % os.path.join(self.advDir,'vcf2maf.py'),
			'-v %s' % os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.indel.reformated.vcf.gz'),
			'-o %s' % os.path.join(self.somaticsnvDirsam,self.sampleID+'.samtools.somatic.indel.maf'),
			'-t %s -n %s -s %s -p %s -m samtools -f avsnp142 -x caddindel' % (self.sampleID,self.Nname,self.seqstrag,self.platform)])
		return ' && \\\n'.join([step1,step2,step3,step4])+'\n', order
	
	def varscan(self):
		#order = 'order finalbam_%s before varscan_%s' % (self.sampleID,self.sampleID)
		order = ['order varscan_%s after finalbam_%s' % (self.sampleID,self.sampleID),
			'order varscan_%s after finalbam_%s' % (self.sampleID,self.Nname)]
		step1 = '  \\\n\t'.join([self.softwares['somatic_varScan'],
			'-r %s' % self.TR,
			self.sampleID,
			self.Nbam,
			self.Tbam,
			self.somaticsnvDir])
		if self.flank:
			step1 = '  \\\n\t'.join([self.softwares['somatic_varScan'],
				self.sampleID,
				self.Nbam,
				self.Tbam,
				self.somaticsnvDir])
		step2 = 'mv %s %s' % (os.path.join(self.somaticsnvDir,self.sampleID+'.varScan.somatic.snp.vcf'), \
				os.path.join(self.somaticsnvDir,self.sampleID+'.varScan.somatic.snv.vcf'))
		return ' && \\\n'.join([step1,step2])+'\n',order
	
	def somaticvarscanannovar(self):
		#order = 'order varscan_%s before somaticvarscanannovar_%s' % (self.sampleID,self.sampleID)
		order = 'order somaticvarscanannovar_%s after varscan_%s' % (self.sampleID,self.sampleID)
		somatic_snv_db1 = self.somatic_snvdb + ',caddgt10'
		somatic_snv_db2 = self.somatic_snvdb + ',caddindel'
		if self.genome_info['annovarbuild'] == 'mm10':
			somatic_snv_db1 = somatic_snv_db2 = 'GeneName,refGene,genomicSuperDups,gff3,snp142chr' 
		step1 = '  \\\n\t'.join([self.annovar,
			'-a somatic',
			'-r %s' % self.refData,
			'-u %s' % self.genome_info['annovarbuild'],
			'-b SNP',
			'-z %s' % somatic_snv_db1,
			'-d %s' % self.humandbDir,
			os.path.join(self.somaticsnvDir,self.sampleID+'.varScan.somatic.snv.vcf'), 
			'%s,%s' % (self.Nname,self.sampleID)])
		step2 = '  \\\n\t'.join([self.annovar,
			'-a somatic',
			'-r %s' % self.refData,
			'-u %s' % self.genome_info['annovarbuild'],
			'-b INDEL',
			'-z %s' % somatic_snv_db2,
			'-d %s' % self.humandbDir,
			os.path.join(self.somaticsnvDir,self.sampleID+'.varScan.somatic.indel.vcf'),
			'%s,%s' % (self.Nname,self.sampleID)])
		step3 = ' \\\n\t'.join(['python %s' % os.path.join(self.advDir,'vcf2maf.py'),
			'-v %s' % os.path.join(self.somaticsnvDir,self.sampleID+'.varScan.somatic.snv.reformated.vcf.gz'),
			'-o %s' % os.path.join(self.somaticsnvDir,self.sampleID+'.varScan.somatic.snv.maf'),
			'-t %s -n %s -s %s -p %s -m varScan -f avsnp142 -x caddgt10' % (self.sampleID,self.Nname,self.seqstrag,self.platform)])
		step4 = ' \\\n\t'.join(['python %s' % os.path.join(self.advDir,'vcf2maf.py'),
			'-v %s' % os.path.join(self.somaticsnvDir,self.sampleID+'.varScan.somatic.indel.reformated.vcf.gz'),
			'-o %s' % os.path.join(self.somaticsnvDir,self.sampleID+'.varScan.somatic.indel.maf'),
			'-t %s -n %s -s %s -p %s -m varScan -f avsnp142 -x caddindel' % (self.sampleID,self.Nname,self.seqstrag,self.platform)])
		return ' && \\\n'.join([step1,step2,step3,step4])+'\n',order
	
	def exomeCNV_Ndepth(self):
		order1 = 'order exomeCNV_Ndepth_%s after finalbam_%s' % (self.sampleID,self.Nname)
		step1 = ' \\\n\t'.join([self.softwares['java7']+' -Djava.io.tmpdir=%s -Xmx5g -jar %s' % \
			(self.somaticcnvDire,os.path.join(self.gatkDir,'GenomeAnalysisTK.jar')),
			'-T DepthOfCoverage -rf BadCigar',
			'-mmq 1',
			'-mbq 13',
			'-omitLocusTable',
			'-omitBaseOutput',
			'-R %s' % self.refData,
			'-I %s' % self.Nbam,
			'-L %s' % self.TR,
			'-o %s' % os.path.join(self.somaticcnvDire,self.sampleID+'.normal')])
		return step1, order1
	
	def exomeCNV_Tdepth(self):
		order1 = 'order exomeCNV_Tdepth_%s after finalbam_%s' % (self.sampleID,self.sampleID)
		step1 = ' \\\n\t'.join([self.softwares['java7']+' -Djava.io.tmpdir=%s -Xmx5g -jar %s' % \
			(self.somaticcnvDire,os.path.join(self.gatkDir,'GenomeAnalysisTK.jar')),
			'-T DepthOfCoverage -rf BadCigar',
			'-mmq 1',
			'-mbq 13',
			'-omitLocusTable',
			'-omitBaseOutput',
			'-R %s' % self.refData,
			'-I %s' % self.Tbam,
			'-L %s' % self.TR,
			'-o %s' % os.path.join(self.somaticcnvDire,self.sampleID+'.tumor')])
		return step1, order1
		
	
	def exomeCNV_CNV(self):
		order1 = 'order exomeCNV_CNV_%s after exomeCNV_Tdepth_%s' % (self.sampleID,self.sampleID)
		order2 = 'order exomeCNV_CNV_%s after exomeCNV_Ndepth_%s' % (self.sampleID,self.sampleID)
		somatic_svdb = self.somatic_svdb
		if self.genome_info['annovarbuild'] == 'mm10':
			somatic_svdb = 'GeneName,refGene,genomicSuperDups,gff3,snp142chr'
		step0 = 'cd %s' % self.somaticcnvDire
		step1 = '%s %s -n %s -t %s -i %s -l 150 -r %s -o %s' % (self.softwares['r3.0.3'], \
			os.path.join(self.softwares['exomecnv'],'exomeCNV.CNV.call.CBS.R'), \
			os.path.join(self.somaticcnvDire,self.sampleID+'.normal.sample_interval_summary'), \
			os.path.join(self.somaticcnvDire,self.sampleID+'.tumor.sample_interval_summary'), \
			self.sampleID,'0.3',self.somaticcnvDire)
		step2 = '%s %s %s %s' % (self.softwares['r3.0.3'], \
			os.path.join(self.softwares['exomecnv'],'exomeCNV.CNV.call.CBS2gff.R'), \
			os.path.join(self.somaticcnvDire,self.sampleID+'.CNV.CBS.txt'),
			os.path.join(self.somaticcnvDire,self.sampleID+'.ExomeCNV.somatic.cnv.gff'))
		step3 = 'cp %s %s ' % (os.path.join(self.somaticcnvDire,self.sampleID+'.CNV.CBS.txt'), \
			os.path.join(self.somaticcnvDire,self.sampleID+'.somatic.CNV.txt'))
		step4 = '  \\\n\t'.join([self.annovar,
			'-z %s' % somatic_svdb,
			'-r %s' % self.refData,
			'-u %s' % self.genome_info['annovarbuild'],
			'-t copy.number',
			'-d %s' % self.humandbDir,
			os.path.join(self.somaticcnvDire,self.sampleID+'.ExomeCNV.somatic.cnv.gff'),
			self.sampleID])
		step5 = 'sv_cnv.stat.cancer.pl %s/%s.ExomeCNV.somatic.cnv.hg19_multianno.xls >%s/%s.ExomeCNV.somatic.cnv.stat.xls' % \
			(self.somaticcnvDire,self.sampleID,self.somaticcnvDire,self.sampleID)
		return ' && \\\n'.join([step1,step2,step3,step4,step5]), [order1,order2]
		
	
	def somatic_muTect(self):
		#order = ['order finalbam_%s before somatic_muTect_%s' % (self.Nname,self.sampleID),
		#	'order finalbam_%s before somatic_muTect_%s' % (self.sampleID,self.sampleID)]
		order = ['order somatic_muTect_%s after finalbam_%s' % (self.sampleID,self.Nname),
			'order somatic_muTect_%s after finalbam_%s' % (self.sampleID,self.sampleID)]
		cal_param = ['--dbsnp %s' % self.genome_info['dbsnp.mutect']]
		if 'cosmic.mutect' in self.genome_info:
			cal_param.append('--cosmic %s' % self.genome_info['cosmic.mutect'])
		cal_param.append('>/dev/null')
		if self.genome_info['annovarbuild'] == 'mm10':
			cal_param = []
		step1 = '  \\\n\t'.join([self.softwares['java6']+' -Xmx6G -Djava.io.tmpdir=%s -jar %s/muTect.jar' % (self.somaticsnvDirm,self.muTectDir),
			'-T MuTect',
			'-rf BadCigar',
			'-dt NONE',
			'--reference_sequence %s' % self.refData,
			'--intervals %s' % self.TR,
			'--input_file:normal %s' % self.Nbam,
			'--input_file:tumor %s' % self.Tbam,
			'--vcf %s' % os.path.join(self.somaticsnvDirm,self.sampleID+'.muTect.call.vcf')] + \
			cal_param)
#			'--coverage_file %s' % os.path.join(self.somaticsnvDirm,'coverage.wig.txt')])
		if self.flank:
			step1 = '  \\\n\t'.join([self.softwares['java6']+' -Xmx6G -Djava.io.tmpdir=%s -jar %s/muTect.jar' % (self.somaticsnvDirm,self.muTectDir),
				'-T MuTect',
				'-rf BadCigar',
				'-dt NONE',
				'--reference_sequence %s' % self.refData,
				'--input_file:normal %s' % self.Nbam,
				'--input_file:tumor %s' % self.Tbam,
				'--vcf %s' % os.path.join(self.somaticsnvDirm,self.sampleID+'.muTect.call.vcf')] + \
				cal_param)
		step2 = 'grep -v REJECT  %s | awk -F "\\t" \'$1!="hs37d5"\' > %s' % (os.path.join(self.somaticsnvDirm,self.sampleID+'.muTect.call.vcf'),os.path.join(self.somaticsnvDirm,self.sampleID+'.muTect.somatic.snv.vcf'))
		step2_1 = 'rm -rf %s' % (os.path.join(self.somaticsnvDirm,self.sampleID+'.muTect.call.vcf'))
		step3 = 'names=`grep "#CHROM" %s | awk \'{print $(NF-1)","$NF}\'`' % os.path.join(self.somaticsnvDirm,self.sampleID+'.muTect.somatic.snv.vcf')
		step4 = '  \\\n\t'.join([self.annovar,
			'-a somatic',
			'-b SNP',
			'-z %s,caddgt10' % self.somatic_snvdb,
			'-d %s' % self.humandbDir,
			os.path.join(self.somaticsnvDirm,self.sampleID+'.muTect.somatic.snv.vcf'),
			'$names'])
		step5 = ' \\\n\t'.join(['python %s' % os.path.join(self.advDir,'vcf2maf.py'),
			'-v %s' % os.path.join(self.somaticsnvDirm,self.sampleID+'.muTect.somatic.snv.reformated.vcf.gz'),
			'-o %s' % os.path.join(self.somaticsnvDirm,self.sampleID+'.muTect.somatic.snv.maf'),
			'-t %s -n %s -s %s -p %s -m muTect -f avsnp142 -x caddgt10' % (self.sampleID,self.Nname,self.seqstrag,self.platform)])

		return ' && \\\n'.join([step1,step2,step2_1,step3,step4,step5]),order
	
	
	
	def somatic_muTect_sub(self,chroms,i,FFPE=False):
		#order = ['order finalbam_%s before somatic_muTect_%s' % (self.Nname,self.sampleID),
		#	'order finalbam_%s before somatic_muTect_%s' % (self.sampleID,self.sampleID)]
		order = ['order somatic_muTect_sub%s_%s after finalbam_%s' % (i,self.sampleID,self.Nname),
			'order somatic_muTect_sub%s_%s after finalbam_%s' % (i,self.sampleID,self.sampleID)]
		cal_param = ['--dbsnp %s' % self.genome_info['dbsnp.mutect']]
		if 'cosmic.mutect' in self.genome_info:
			cal_param.append('--cosmic %s' % self.genome_info['cosmic.mutect'])
		cal_param.append('>/dev/null')
		cmds = []
		if self.genome_info['annovarbuild'] == 'mm10':
			cal_param = []
		if FFPE:
			for each_chr in chroms:
				cmd_tmp = ' \\\n\t'.join([self.softwares['java6']+' -Xmx6G -Djava.io.tmpdir=%s -jar %s/muTect.jar' % (self.somaticsnvDirm,self.muTectDir),
					'-T MuTect --max_alt_alleles_in_normal_qscore_sum 5  --initial_tumor_lod 8 --tumor_lod 1.3 --fraction_contamination 0.02 --required_maximum_alt_allele_mapping_quality_score 40 --strand_artifact_power_threshold 0.05 --gap_events_threshold 8 --pir_median_threshold 10 --pir_mad_threshold 3.0  -rf BadCigar -dt NONE -L %s' % each_chr,
					'--reference_sequence %s' % self.refData,
					'--intervals %s' % self.TR,
					'--input_file:normal %s' % self.Nbam,
					'--input_file:tumor %s' % self.Tbam,
					'--vcf %s' % os.path.join(self.somaticsnvDirm,self.sampleID+'.%s.muTect.call.vcf' % each_chr)] + \
					cal_param)
				if self.flank:
					cmd_tmp =' \\\n\t'.join([self.softwares['java6']+' -Xmx6G -Djava.io.tmpdir=%s -jar %s/muTect.jar' % (self.somaticsnvDirm,self.muTectDir),
						'-T MuTect --max_alt_alleles_in_normal_qscore_sum 5  --initial_tumor_lod 8 --tumor_lod 1.3 --fraction_contamination 0.02 --required_maximum_alt_allele_mapping_quality_score 40 --strand_artifact_power_threshold 0.05 --gap_events_threshold 11 --pir_median_threshold 10 --pir_mad_threshold 3.0  -rf BadCigar -dt NONE -L %s' % each_chr,
						'--reference_sequence %s' % self.refData,
						'--input_file:normal %s' % self.Nbam,
						'--input_file:tumor %s' % self.Tbam,
						'--vcf %s' % os.path.join(self.somaticsnvDirm,self.sampleID+'.%s.muTect.call.vcf' % each_chr)] + \
						cal_param)
				cmds.append(cmd_tmp)
		else:
			for each_chr in chroms:
				cmd_tmp = ' \\\n\t'.join([self.softwares['java6']+' -Xmx6G -Djava.io.tmpdir=%s -jar %s/muTect.jar' % (self.somaticsnvDirm,self.muTectDir),
					'-T MuTect -rf BadCigar -dt NONE -L %s' % each_chr,
					'--reference_sequence %s' % self.refData,
					'--intervals %s' % self.TR,
					'--input_file:normal %s' % self.Nbam,
					'--input_file:tumor %s' % self.Tbam,
					'--vcf %s' % os.path.join(self.somaticsnvDirm,self.sampleID+'.%s.muTect.call.vcf' % each_chr)] + \
					cal_param)
				if self.flank:
					cmd_tmp =' \\\n\t'.join([self.softwares['java6']+' -Xmx6G -Djava.io.tmpdir=%s -jar %s/muTect.jar' % (self.somaticsnvDirm,self.muTectDir),
						'-T MuTect -rf BadCigar -dt NONE -L %s' % each_chr,
						'--reference_sequence %s' % self.refData,
						'--input_file:normal %s' % self.Nbam,
						'--input_file:tumor %s' % self.Tbam,
						'--vcf %s' % os.path.join(self.somaticsnvDirm,self.sampleID+'.%s.muTect.call.vcf' % each_chr)] + \
						cal_param)
				cmds.append(cmd_tmp)
		return ' && \\\n'.join(cmds)+'\n',order


	def somatic_muTect_catannovar(self,subs,FFPE=False):
		order = ['order somatic_muTect_catannovar_%s after somatic_muTect_sub%s_%s' % \
			(self.sampleID,each,self.sampleID) for each in subs]
		chr_vcfs = [os.path.join(self.somaticsnvDirm,self.sampleID+'.%s.muTect.call.vcf' % chr) for chr in self.chroms]
		raw_vcf = os.path.join(self.somaticsnvDirm,self.sampleID+'.muTect.call.vcf')
		flt_vcf = os.path.join(self.somaticsnvDirm,self.sampleID+'.muTect.somatic.snv.vcf')
		somatic_snv_db1 = self.somatic_snvdb + ',caddgt10'
		if self.genome_info['annovarbuild'] == 'mm10':
			somatic_snv_db1 = 'GeneName,refGene,genomicSuperDups,gff3,snp142chr'
		step1 = 'vcf-concat -p '+'\\\n\t'.join(chr_vcfs)+'\\\n\t>'+raw_vcf
		step2 = 'grep -v REJECT %s > %s' % (raw_vcf,flt_vcf)
		step3 = 'names=`grep "#CHROM" %s | awk \'{print $(NF-1)","$NF}\'`' % flt_vcf
		step4 = '  \\\n\t'.join([self.annovar,
			'-a somatic',
			'-r %s' % self.refData,
			'-u %s' % self.genome_info['annovarbuild'],
			'-b SNP',
			'-z %s' % somatic_snv_db1,
			'-d %s' % self.humandbDir,
			flt_vcf,'$names'])
		step5 = ' \\\n\t'.join(['python %s' % os.path.join(self.advDir,'vcf2maf.py'),
			'-v %s' % os.path.join(self.somaticsnvDirm,self.sampleID+'.muTect.somatic.snv.reformated.vcf.gz'),
			'-o %s' % os.path.join(self.somaticsnvDirm,self.sampleID+'.muTect.somatic.snv.maf'),
			'-t %s -n %s -s %s -p %s -m muTect -f avsnp142 -x caddgt10' % (self.sampleID,self.Nname,self.seqstrag,self.platform)])
		step6 = 'rm -rf %s \\\n  %s \\\n  %s' % ('\\\n  '.join(chr_vcfs),'\\\n  '.join(['%s.idx'%each for each in chr_vcfs]), raw_vcf)
		step_all = [step1,step2,step3,step4,step5,step6]
		if FFPE:
			filter_list = os.path.join(self.somaticsnvDirm,self.sampleID+'.filter.list')
			step7 = 'echo -e "%s\\t%s\\t%s" > %s'%(self.sampleID,os.path.join(self.somaticsnvDirm,self.sampleID+'.muTect.somatic.snv.annovar.hg19_multianno.xls'),self.Tmpilegz,filter_list)
			step_all.append(step7)
			step8 = 'python %s -i %s -t \'C>T,G>A\' -o %s' %(os.path.join(self.ffpe,'CT.Filter.py'),filter_list,self.somaticsnvDirm)
			step_all.append(step8)
		return ' && \\\n'.join(step_all)+'\n',order	
		#return ' && \\\n'.join([step1,step2,step3,step4,step5,step6])+'\n',order
	
	def somatic_Strelka(self):
		#order = ['order finalbam_%s before somatic_Strelka_%s' % (self.Nname,self.sampleID),
		#	'order finalbam_%s before somatic_Strelka_%s' % (self.sampleID,self.sampleID)]		
		order = ['order somatic_Strelka_%s after finalbam_%s' % (self.sampleID,self.Nname),
			'order somatic_Strelka_%s after finalbam_%s' % (self.sampleID,self.sampleID)]		
		if 'WGS' in self.seqstrag:
			config = self.strelkaDir+'/etc/strelka_config_bwa_default_WGS.ini'
		else:
			config = self.strelkaDir+'/etc/strelka_config_bwa_default_TR.ini'
		somatic_snv_db1 = self.somatic_snvdb + ',caddgt10'
		somatic_snv_db2 = self.somatic_snvdb + ',caddindel'
		if self.genome_info['annovarbuild'] == 'mm10':
			somatic_snv_db1 = somatic_snv_db2 = 'GeneName,refGene,genomicSuperDups,gff3,snp142chr'
		step0 = 'if [ -d "%s/build" ];then rm -rf %s/build; fi' % (self.somaticsnvDirs,self.somaticsnvDirs)
		step1 = '  \\\n\t'.join([self.strelkaDir+'/bin/configureStrelkaWorkflow.pl',
			'--tumor=%s' % self.Tbam,
			'--normal=%s' % self.Nbam,
			'--ref=%s' % self.refData,
			'--config=%s' % config,
			'--output-dir=%s' % self.somaticsnvDirs+'/build'])
		step2 = 'make -C %s -j 8' % (self.somaticsnvDirs+'/build')		
		step3 = 'awk -F "\\t" \'$1!="hs37d5"\' %s/results/passed.somatic.indels.vcf > %s/%s.Strelka.somatic.indel.vcf' % (self.somaticsnvDirs+'/build',self.somaticsnvDirs,self.sampleID)
		step4 = 'awk -F "\\t" \'$1!="hs37d5"\' %s/results/passed.somatic.snvs.vcf > %s/%s.Strelka.somatic.snv.vcf' % (self.somaticsnvDirs+'/build',self.somaticsnvDirs,self.sampleID)
		step5 = '  \\\n\t'.join([self.annovar,
			'-a somatic',
			'-r %s' % self.refData,
			'-u %s' % self.genome_info['annovarbuild'],
			'-b SNP',
			'-z %s' % somatic_snv_db1,
			'-d %s' % self.humandbDir,
			'-v vcf4old',
			os.path.join(self.somaticsnvDirs,self.sampleID+'.Strelka.somatic.snv.vcf'),
			'%s,%s' % (self.Nname,self.sampleID)])
		step6 = '  \\\n\t'.join([self.annovar,
			'-a somatic',
			'-r %s' % self.refData,
			'-u %s' % self.genome_info['annovarbuild'],
			'-b INDEL',
			'-z %s' % somatic_snv_db2,
			'-d %s' % self.humandbDir,
			'-v vcf4old',
			os.path.join(self.somaticsnvDirs,self.sampleID+'.Strelka.somatic.indel.vcf'),
			'%s,%s' % (self.Nname,self.sampleID)])
		step7 = ' \\\n\t'.join(['python %s' % os.path.join(self.advDir,'vcf2maf.py'),
			'-v %s' % os.path.join(self.somaticsnvDirs,self.sampleID+'.Strelka.somatic.snv.reformated.vcf.gz'),
			'-o %s' % os.path.join(self.somaticsnvDirs,self.sampleID+'.Strelka.somatic.snv.maf'),
			'-t %s -n %s -s %s -p %s -m Strelka -f avsnp142 -x caddgt10' % (self.sampleID,self.Nname,self.seqstrag,self.platform)])
		step8 = ' \\\n\t'.join(['python %s' % os.path.join(self.advDir,'vcf2maf.py'),
			'-v %s' % os.path.join(self.somaticsnvDirs,self.sampleID+'.Strelka.somatic.indel.reformated.vcf.gz'),
			'-o %s' % os.path.join(self.somaticsnvDirs,self.sampleID+'.Strelka.somatic.indel.maf'),
			'-t %s -n %s -s %s -p %s -m Strelka -f avsnp142 -x caddindel' % (self.sampleID,self.Nname,self.seqstrag,self.platform)])
		step9 = 'rm -rf %s/build' % self.somaticsnvDirs
		
		return ' && \\\n'.join([step0,step1,step2,step3,step4,step5,step6,step7,step8,step9])+'\n',order
		
	def sobam2cfg(self):
		#order = ['order finalbam_%s before sobam2cfg_%s' % (self.sampleID,self.sampleID),
		#	'order finalbam_%s before sobam2cfg_%s' % (self.Nname,self.sampleID)]
		order = ['order sobam2cfg_%s after finalbam_%s' % (self.sampleID,self.sampleID),
			'order sobam2cfg_%s after finalbam_%s' % (self.sampleID,self.Nname)]
		step1 = 'cd %s' % self.somaticsvDirb
		step2 = ' \\\n\t'.join(['perl '+os.path.join(self.breakdancerDir,'bam2cfg.pl'),
			'-g -h',
			self.Tbam,
			self.Nbam,
			'> '+os.path.join(self.somaticsvDirb,self.sampleID+'.breakdancer.cfg')])
		return ' && \\\n'.join([step1,step2])+'\n',order
		
		
	def sobreakdancerSV(self):
		#order = 'order sobam2cfg_%s before sobreakdancerSV_%s' % (self.sampleID,self.sampleID)
		order = 'order sobreakdancerSV_%s after sobam2cfg_%s' % (self.sampleID,self.sampleID)
		
		step0 = 'cd %s' % self.somaticsvDirb
		step1 = '  \\\n\t'.join([os.path.join(self.breakdancerDir,'breakdancer-max'),
			'-h',
			'-d '+os.path.join(self.somaticsvDirb,self.sampleID+'.breakdancer.SV-supporting'),
			'-g '+os.path.join(self.somaticsvDirb,self.sampleID+'.breakdancer.bed'),
			os.path.join(self.somaticsvDirb,self.sampleID+'.breakdancer.cfg'),
			'> '+os.path.join(self.somaticsvDirb,self.sampleID+'.breakdancer.txt')])
		step2 = '  \\\n\t'.join(['var_sv_breakdancer.filter.pl',
				'-g M -n 6',
				'-b %s' % self.Nbam,
				'-a '+os.path.join(self.somaticsvDirb,self.sampleID+'.breakdancer.txt'),
				'> '+os.path.join(self.somaticsvDirb,self.sampleID+'.breakdancer.flt.txt')])
		step3 = '  \\\n\t'.join(['var_sv_breakdancer.toGff.pl',
				os.path.join(self.somaticsvDirb,self.sampleID+'.breakdancer.flt.txt'),
				'> '+os.path.join(self.somaticsvDirb,self.sampleID+'.breakdancer.somatic.sv.gff')])
		return ' && \\\n'.join([step0,step1,step2,step3])+'\n',order

	def sobreakdancerSvAnno(self):
		#order = 'order sobreakdancerSV_%s before sobreakdancerSvAnno_%s' % (self.sampleID,self.sampleID)
		order = 'order sobreakdancerSvAnno_%s after sobreakdancerSV_%s' % (self.sampleID,self.sampleID)
		somatic_sv_db = self.somatic_svdb
		if self.genome_info['annovarbuild'] == 'mm10':
			somatic_sv_db = 'GeneName,refGene,genomicSuperDups,gff3,snp142chr'
		step1 = '  \\\n\t'.join([self.annovar,
			'-r %s' % self.refData,
			'-u %s' % self.genome_info['annovarbuild'],
			'-z %s' % somatic_sv_db,
			'-d %s' % self.humandbDir,
			os.path.join(self.somaticsvDirb,self.sampleID+'.breakdancer.somatic.sv.gff'),
			self.sampleID])
		step2 = '  \\\n\t'.join(['sv_cnv.stat.cancer.pl',
			'-s '+self.sampleID,
			os.path.join(self.somaticsvDirb,self.sampleID+'.breakdancer.somatic.sv.hg19_multianno.xls'),
			'> '+os.path.join(self.somaticsvDirb,self.sampleID+'.breakdancer.somatic.sv.stat.xls')])
		return ' && \\\n'.join([step1, step2])+'\n',order

	def diff_sclip(self,wtffpe):
		order = ['order diff_sclip_%s after catSoftClip_%s' % (self.sampleID,self.Nname),
			'order diff_sclip_%s after finalbam_%s' % (self.sampleID,self.Nname),
			'order diff_sclip_%s after catSoftClip_%s' % (self.sampleID,self.sampleID),
			'order diff_sclip_%s after finalbam_%s' % (self.sampleID,self.sampleID)]
		final_somatic_cover = os.path.join(self.somaticsvDirc,self.sampleID+'.somatic.cover')
		step1 = ' \\\n\t'.join(['perl '+os.path.join(self.crestDir,'CREST','countDiff.pl'),
			'-d '+os.path.join(self.mapDir,self.sampleID,self.sampleID+'.cover'),
			'-g '+os.path.join(self.mapDir,self.Nname,self.Nname+'.cover'),
			'> '+final_somatic_cover])
		step2 = 'echo "No filter"'
		if wtffpe == 'Y':
			step2 = 'count=$(less %s|wc -l)' % final_somatic_cover +'\n'+ 'if [ $count -gt 10000000 ];then python %s -i %s -o %s && mv %s %s && mv %s %s;else echo "No filter cover";fi' % ('Cancer/script/Cover_Filter.py',final_somatic_cover,os.path.join(self.somaticsvDirc,self.sampleID+'.somatic.flt.cover'),final_somatic_cover,os.path.join(self.somaticsvDirc,self.sampleID+'.somatic.raw.cover'),os.path.join(self.somaticsvDirc,self.sampleID+'.somatic.flt.cover'),final_somatic_cover)
		return ' && \\\n'.join([step1,step2]) + '\n',order
	
	def crest_somaticSV(self):
		order = 'order crest_somaticSV_%s after diff_sclip_%s' % (self.sampleID,self.sampleID)
		step1 = '\\\n\t'.join(['perl '+os.path.join(self.crestDir,'crest_sv_calling.pl'),
			'-cov ' + os.path.join(self.somaticsvDirc,self.sampleID+'.somatic.cover'),
			'-outDir ' + self.somaticsvDirc,
			'-tumorBam ' + self.Tbam,
			'-normalBam ' + self.Nbam,
			'-sampleID ' + self.sampleID+'.crest.somatic.sv',
			'-regionList %s' % self.genome_info['chrbed'],
			'-ref ' + self.refData,
			'-bit ' + self.refData2bit])
		return step1+'\n',order

	def crest_somaticSV_sub(self,chroms,i):
		order = 'order somatic_crest_sub%s_%s after diff_sclip_%s' % (i,self.sampleID,self.sampleID)
		cmds = []
		bychr_suf = '.region.bed'
		if self.genome_info['build'] != 'b37':
			bychr_suf = '.bed'
		for each_chr in chroms:
			cmd_tmp = '\\\n\t'.join(['perl '+os.path.join(self.crestDir,'crest_sv_calling.pl'),
			'-cov ' + os.path.join(self.somaticsvDirc,self.sampleID+'.somatic.cover'),
			'-outDir ' + self.somaticsvDirc,
			'-tumorBam ' + self.Tbam,
			'-normalBam ' + self.Nbam,
			'-sampleID ' + self.sampleID+'.crest.somatic.sv'+'.'+each_chr,
			'-regionList %s' % os.path.join(self.genome_info['bychr'],each_chr + bychr_suf),
			'-ref ' + self.refData,
			'-bit ' + self.refData2bit])
			cmds.append(cmd_tmp)
		return ' && \\\n'.join(cmds)+'\n',order

	def crest_somaticSVann(self):
		gff = os.path.join(self.somaticsvDirc,self.sampleID+'.crest.somatic.sv.gff')
		ann = os.path.join(self.somaticsvDirc,self.sampleID+'.crest.somatic.sv.hg19_multianno.xls')
		stat = os.path.join(self.somaticsvDirc,self.sampleID+'.crest.somatic.sv.stat.xls')

		order = 'order crest_somaticSVann_%s after crest_somaticSV_%s' % (self.sampleID,self.sampleID)
		somatic_sv_db = self.somatic_svdb
		if self.genome_info['annovarbuild'] == 'mm10':
			somatic_sv_db = 'GeneName,refGene,genomicSuperDups,gff3,snp142chr'
		step1 = '  \\\n\t'.join([self.annovar,
			'-r %s' % self.refData,
			'-u %s' % self.genome_info['annovarbuild'],
			'-z %s' % somatic_sv_db,
			'-d %s' % self.humandbDir,
			os.path.join(self.somaticsvDirc,self.sampleID+'.crest.somatic.sv.gff'),
			self.sampleID])
		step2 = 'sv_cnv.stat.cancer.pl -s %s %s > %s' %(self.sampleID,ann,stat)
		return ' && \\\n'.join([step1,step2])+'\n',order

	def somatic_crest_catannovar(self,subs):
		ann = os.path.join(self.somaticsvDirc,self.sampleID+'.crest.somatic.sv.hg19_multianno.xls')
		stat = os.path.join(self.somaticsvDirc,self.sampleID+'.crest.somatic.sv.stat.xls')
		order = ['order somatic_crest_catannovar_%s after somatic_crest_sub%s_%s' % \
			(self.sampleID,each,self.sampleID) for each in subs]
		somatic_sv_db = self.somatic_svdb
                if self.genome_info['annovarbuild'] == 'mm10':
                        somatic_sv_db = 'GeneName,refGene,genomicSuperDups,gff3,snp142chr'
		chr_gff = [os.path.join(self.somaticsvDirc,self.sampleID+'.crest.somatic.sv.%s.gff' % chr) for chr in self.chroms]
		chr_txt = [os.path.join(self.somaticsvDirc,self.sampleID+'.crest.somatic.sv.%s.predSV.txt' % chr) for chr in self.chroms]
		#step1 = 'cat '+'\\\n\t'.join(chr_gff)+'\\\n\t>'+ os.path.join(self.somaticsvDirc,(self.sampleID+'.crest.somatic.sv.gff'))
		step1 = 'cat '+'\\\n\t'.join(chr_txt)+'\\\n\t>'+ os.path.join(self.somaticsvDirc,(self.sampleID+'.crest.somatic.sv.predSV.txt'))
		step2 = 'perl share/software/CREST/var_sv_CREST.toGff.pl ' + os.path.join(self.somaticsvDirc,(self.sampleID+'.crest.somatic.sv.predSV.txt')) + ' > '+ os.path.join(self.somaticsvDirc,(self.sampleID+'.crest.somatic.sv.gff'))
		step3 = '  \\\n\t'.join([self.annovar,
		'-r %s' % self.refData,
		'-u %s' % self.genome_info['annovarbuild'],
		'-z %s' % somatic_sv_db,
		'-d %s' % self.humandbDir,
		os.path.join(self.somaticsvDirc,self.sampleID+'.crest.somatic.sv.gff'),
		self.sampleID])
		step4 = 'sv_cnv.stat.cancer.pl -s %s %s > %s' %(self.sampleID,ann,stat)
		step5 = 'rm -rf %s \\\n  %s' % ('\\\n  '.join(chr_gff),'\\\n  '.join(chr_txt))
		return ' && \\\n'.join([step1,step2,step3,step4,step5])+'\n',order

	def lumpy_somaticSV(self):
		order = ['order lumpy_somaticSV_%s after finalbam_%s' % (self.sampleID,self.sampleID),
			'order lumpy_somaticSV_%s after finalbam_%s' % (self.sampleID,self.Nname),
			'order lumpy_somaticSV_%s after sambamba_mergesplit_%s' % (self.sampleID,self.sampleID),
			'order lumpy_somaticSV_%s after sambamba_mergesplit_%s' % (self.sampleID,self.Nname)]
		step1 = ' \\\n\t'.join([os.path.join(self.speedseqDir,'bin','lumpyexpress'),
			'-B %s,%s' % (self.Tbam,self.Nbam),
			'-S %s,%s' % (self.Tsbam,self.Nsbam),
			'-D %s,%s' % (self.Tdbam,self.Ndbam),
			'-o %s' % os.path.join(self.somaticsvDirl,self.sampleID+'.lumpy.somatic.sv.tmp.vcf'),
			'-x %s' % self.genome_info['lumpy-exclude'],
			'-T %s' % self.somaticsvDirl,
			'-K %s' % os.path.join(self.speedseqDir,'bin','speedseq.config'),
			'-P -v -k'])
		step2 = 'python Cancer/summary/lumpy_filter1.py -vcf %s -o %s' % (os.path.join(self.somaticsvDirl,self.sampleID+'.lumpy.somatic.sv.tmp.vcf'),os.path.join(self.somaticsvDirl,self.sampleID+'.lumpy.somatic.sv.vcf'))
		step3 = ' \\\n\t'.join(['python %s/lumpy_vcf2gff.py' % self.advDir,
			os.path.join(self.somaticsvDirl,self.sampleID+'.lumpy.somatic.sv.vcf'),
			os.path.join(self.somaticsvDirl,self.sampleID+'.lumpy.somatic.sv.gff')])
		step4 = 'gzip %s' % os.path.join(self.somaticsvDirl,self.sampleID+'.lumpy.somatic.sv.vcf')
		return ' && \\\n'.join([step1,step2,step3,step4])+'\n',order

	def lumpy_somaticSVann(self):
		order = 'order lumpy_somaticSVann_%s after lumpy_somaticSV_%s' % (self.sampleID,self.sampleID)
		ann = os.path.join(self.somaticsvDirl,self.sampleID+'.lumpy.somatic.sv.hg19_multianno.xls')
		stat = os.path.join(self.somaticsvDirl,self.sampleID+'.lumpy.somatic.sv.stat.xls')
		somatic_sv_db = self.somatic_svdb
		if self.genome_info['annovarbuild'] == 'mm10':
			somatic_sv_db = 'GeneName,refGene,genomicSuperDups,gff3,snp142chr'
		step1 = ' \\\n\t'.join([self.annovar,
			'-r %s' % self.refData,
			'-u %s' % self.genome_info['annovarbuild'],
			'-z %s' % somatic_sv_db,
			'-d %s' % self.humandbDir,
			os.path.join(self.somaticsvDirl,self.sampleID+'.lumpy.somatic.sv.gff'),
			self.sampleID])
		step2 = 'sv_cnv.stat.cancer.pl -s %s %s > %s' %(self.sampleID,ann,stat)
		return ' && \\\n'.join([step1,step2])+'\n',order
		

	def freec_somaticCNV(self,sex='XY',step=1000):
		cfg = '''
[general]
BedGraphOutput = FALSE
breakPointType = 2
chrFiles = %s
chrLenFile = %s
maxThreads = 6
ploidy = 2
gemMappabilityFile = %s
uniqueMatch = TRUE
sex = %s
telocentromeric = 50000
forceGCcontentNormalization = %d
window = %d
step = %d
outputDir = %s
breakPointThreshold = %s
noisyData = %s
readCountThreshold = %d
printNA = %s

[sample]
mateFile = %s
inputFormat = pileup
mateOrientation = FR

[control]
mateFile = %s
inputFormat = pileup
mateOrientation = FR

[BAF]
SNPfile = %s
minimalCoveragePerPosition = 5
shiftInQuality = 33
'''
		chr_flag = ''
		if not self.genome_info['build'].startswith('b'):
			chr_flag = ' --chr'
		fGC = 0; window = step*2; step = step; bp = 1.5; noisy = 'TRUE'; rc = 50; pna = 'FALSE'
		if 'WGS' in self.seqstrag:
			fGC = 1; bp = 0.8; noisy = 'FALSE'; rc = 10; pna = 'TRUE'
		config = cfg % (self.genome_info['bychr'],self.genome_info['chrlen'],self.genome_info['mappability'],sex,fGC,window,step,self.somaticcnvDirf,bp,noisy,rc,pna,self.Tmpilegz,self.Nmpilegz,self.genome_info['dbsnp.txt'])
		if not 'WGS' in self.seqstrag:
			config += '\n[target]\ncaptureRegions = %s\n' % self.TR
		open(os.path.join(self.somaticcnvDirf,'somatic_freec.config'),'w').write(config)
		gff = os.path.join(self.somaticcnvDirf,self.sampleID+'.freec.somatic.cnv.gff')

		open(os.path.join(self.somaticcnvDirf,'plot.config'),'w').write('%s\t%s\t%s\n'%(self.sampleID,os.path.join(self.somaticcnvDirf,self.sampleID+'.final.mpileup.gz_ratio.txt'),os.path.join(self.somaticcnvDirf,self.sampleID+'.final.mpileup.gz_BAF.txt')))
		order = ['order freec_somaticCNV_%s after mpileup_%s' % (self.sampleID,self.Nname),
			 'order freec_somaticCNV_%s after mpileup_%s' % (self.sampleID,self.sampleID)]
		step0 = 'sex=`if [ -f "%s" ];then awk \'{if(/F/){print "XX"}else{print "XY"}}\' %s;fi` && \\\n' % \
			(os.path.join(self.statDir,self.sampleID,self.sampleID+'.gender'), os.path.join(self.statDir,self.sampleID,self.sampleID+'.gender'))
		step0 += 'if [ $sex == "XX" ];then sed -i \'s/sex=XY/sex=XX/g\' %s ;sed -i \'s/sex = XY/sex = XX/g\' %s ;fi' % \
			(os.path.join(self.somaticcnvDirf,'somatic_freec.config'), os.path.join(self.somaticcnvDirf,'somatic_freec.config'))
		step1 = 'freec -conf %s' % os.path.join(self.somaticcnvDirf,'somatic_freec.config')
		if 'WGS' in self.seqstrag:
			step2 = 'python %s/plot_cnv_profile.py -c %s -s 20 -g %s -o %s' % \
			(self.advDir,os.path.join(self.somaticcnvDirf,'plot.config'),self.sampleID,self.somaticcnvDirf)
		else:
			step2 = 'python %s/plot_cnv_profile.py -c %s -s 50 -g %s -o %s' % \
			(self.advDir,os.path.join(self.somaticcnvDirf,'plot.config'),self.sampleID,self.somaticcnvDirf)
		step3 = 'python %s/Freec.somatic.format.py -i %s -s %d -o %s -x $sex -b %s' % \
			(self.advDir,os.path.join(self.somaticcnvDirf,self.sampleID+'.final.mpileup.gz_CNVs'), \
			 step, os.path.join(self.somaticcnvDirf,self.sampleID+'.somatic.CNV.txt'),self.genome_info['chrbed']) + chr_flag
		step4 = 'grep somatic %s | grep -v normal | awk -F "\\t" -v OFS="\\t" \'{{print $1,"FREEC",$7,$2,$3,".",".",".","CopyNumber="$5";Size="$4";CNVType="$7";Genotype="$8";GTConfidence="$9}}\' > %s' \
			% (os.path.join(self.somaticcnvDirf,self.sampleID+'.somatic.CNV.txt'),gff)
		step5 = 'rm -rf %s %s %s' % (os.path.join(self.somaticcnvDirf,'GC_profile.cnp'), \
			os.path.join(self.somaticcnvDirf,self.sampleID+'.final.mpileup.gz_sample.cpn'), \
			os.path.join(self.somaticcnvDirf,self.Nname+'.final.mpileup.gz_control.cpn'))
		step6 = 'gzip -f %s \\\n\t%s \\\n\t%s \\\n\t%s' % ( \
			os.path.join(self.somaticcnvDirf,self.sampleID+'.final.mpileup.gz_BAF.txt'), \
			os.path.join(self.somaticcnvDirf,self.sampleID+'.final.mpileup.gz_ratio.txt'), \
			os.path.join(self.somaticcnvDirf,self.sampleID+'.final.mpileup.gz_normal_BAF.txt'), \
			os.path.join(self.somaticcnvDirf,self.sampleID+'.final.mpileup.gz_normal_ratio.txt'))
		return ' && \\\n'.join([step0,step1,step2,step3,step4,step5,step6])+'\n',order

	def freec_somaticCNVann(self):
		order = 'order freec_somaticCNVann_%s after freec_somaticCNV_%s' % (self.sampleID,self.sampleID)
		ann = os.path.join(self.somaticcnvDirf,self.sampleID+'.freec.somatic.cnv.hg19_multianno.xls')
		stat = os.path.join(self.somaticcnvDirf,self.sampleID+'.freec.somatic.cnv.stat.xls')
		somatic_sv_db = self.somatic_svdb
		if self.genome_info['annovarbuild'] == 'mm10':
			somatic_sv_db = 'GeneName,refGene,genomicSuperDups,gff3,snp142chr'
		step1 = '  \\\n\t'.join([self.annovar,
			'-r %s' % self.refData,
			'-u %s' % self.genome_info['annovarbuild'],
			'-z %s' % somatic_sv_db,
			'-t CNVType',
			'-d %s' % self.humandbDir,
			os.path.join(self.somaticcnvDirf,self.sampleID+'.freec.somatic.cnv.gff'),
			self.sampleID])
		step2 = 'sv_cnv.stat.cancer.pl -s %s %s > %s' % (self.sampleID,ann,stat)
		return ' && \\\n'.join([step1,step2])+'\n',order
		
	def varscan_somaticCNV(self):
		order = ['order varscan_somaticCNV_%s after mpileup_%s' % (self.sampleID,self.Nname),
			 'order varscan_somaticCNV_%s after mpileup_%s' % (self.sampleID,self.sampleID)]
		prefix = os.path.join(self.somaticcnvDirv,self.sampleID)
		step1 = 'gzip -dc %s > %s' % (self.Tmpilegz,os.path.join(self.somaticcnvDirv,self.sampleID+'.Tumor.mpileup'))
		step2 = 'gzip -dc %s > %s' % (self.Nmpilegz,os.path.join(self.somaticcnvDirv,self.sampleID+'.Normal.mpileup'))
		step3 = ' \\\n\t'.join(['java -Xmx5G -jar %s copynumber' % os.path.join(self.softwares['varscan'],'VarScan.jar'),
			os.path.join(self.somaticcnvDirv,self.sampleID+'.Normal.mpileup'),
			os.path.join(self.somaticcnvDirv,self.sampleID+'.Tumor.mpileup'),
			os.path.join(self.somaticcnvDirv,self.sampleID+'.varScan')])
		step4 = ' \\\n\t'.join(['java -Xmx5G -jar %s copyCaller' % os.path.join(self.softwares['varscan'],'VarScan.jar'),
			os.path.join(self.somaticcnvDirv,self.sampleID+'.varScan.copynumber'),
			'--output-file %s' % os.path.join(self.somaticcnvDirv,self.sampleID+'.varScan.copyCaller.out')])
		step5 = ' \\\n\t'.join(['%s' % self.softwares['r3.2.1'],
			os.path.join(self.softwares['varscan'],'cbs.segment.R'),
			os.path.join(self.somaticcnvDirv,self.sampleID+'.varScan.copyCaller.out'),
			os.path.join(self.somaticcnvDirv,self.sampleID+'.CBS.segments.txt'),
			self.sampleID])
		step6 = ' \\\n\t'.join(['perl %s' % os.path.join(self.softwares['varscan'],'mergeSegments.pl'),
			os.path.join(self.somaticcnvDirv,self.sampleID+'.CBS.segments.txt'),
			'--ref-arm-sizes %s' % self.genome_info['armsize'],
			'--output-basename %s.merge.segments'%prefix])
		step7 = 'grep -v neutral %s.merge.segments.events.tsv >%s.somatic.CNV.txt'%(prefix,prefix)
		step8 = 'awk -F"\\t" -v OFS="\\t" \'{if(NR>1 && $9>200 ){print $1,"VARSCAN2",$8,$2,$3,".",".",".","CopyNumber="$5";Size="$9";CNVType="$8";size_class="$10";chrom_arm="$11";arm_fraction="$12";chrom_fraction="$13}}\' %s.somatic.CNV.txt > %s.varScan.somatic.cnv.gff' % \
			(prefix,prefix)
		step9 = 'rm -rf %s %s %s' % (os.path.join(self.somaticcnvDirv,self.sampleID+'.Tumor.mpileup'),os.path.join(self.somaticcnvDirv,self.sampleID+'.Normal.mpileup'),os.path.join(self.somaticcnvDirv,self.sampleID+'.varScan.copynumber'))
		return ' && \\\n'.join([step1,step2,step3,step4,step5,step6,step7,step8,step9])+'\n',order

	def varscan_somaticCNVann(self):
		order = 'order varscan_somaticCNVann_%s after varscan_somaticCNV_%s' % (self.sampleID,self.sampleID)
		somatic_sv_db = self.somatic_svdb
		if self.genome_info['annovarbuild'] == 'mm10':
			somatic_sv_db = 'GeneName,refGene,genomicSuperDups,gff3,snp142chr'
		step1 = ' \\\n\t'.join([self.annovar,
			'-r %s' % self.refData,
			'-u %s' % self.genome_info['annovarbuild'],
			'-z %s' % somatic_sv_db,
			'-t CNVType',
			'-d %s' % self.humandbDir,
			os.path.join(self.somaticcnvDirv,self.sampleID+'.varScan.somatic.cnv.gff'),
			self.sampleID])
		step2 = ' \\\n\t'.join(['sv_cnv.stat.cancer.pl',
			'-s %s' % self.sampleID,
			os.path.join(self.somaticcnvDirv,self.sampleID+'.varScan.somatic.cnv.hg19_multianno.xls'),
			'>%s' % os.path.join(self.somaticcnvDirv,self.sampleID+'.varScan.somatic.cnv.stat.xls')])
		return ' && \\\n'.join([step1,step2])+'\n',order

	def spectrum(self,sample_list,spectrumdir,analydir,snvm,indelm,sample_order):
		order = []
		cmds = []
		list = []
		snv_stat = []
		indel_stat = []
		dirtmp = analydir+'/Somatic'
		#### samples in qc_list
		pinfo = {}
		for line in open(os.path.join(analydir,'qc_list')):
			array = line.strip().split('\t')
			if array[1] not in pinfo:
				pinfo[array[1]] = {'T':{},'N':''}
			if array[6] == 'T':
				pinfo[array[1]]['T'][array[2]] = 1
			else:
				pinfo[array[1]]['N'] = array[2]
			
		sample_tmp = []
		for eachp in pinfo:
			if pinfo[eachp]['N'] == '':
				continue
			for each in pinfo[eachp]['T']:
				if each in sample_order:
					if each not in sample_tmp:
						sample_tmp.append(each)
		qc_list = sample_tmp
		########################################################
		if snvm=='samtools' or indelm =='samtools':
			#order = ['order somaticsamtoolsannovar_%s before spectrum' % each for each in sample_list]
			order = ['order spectrum after somaticsamtoolsannovar_%s' % each for each in sample_list]
			list = ['%s %s/%s/somatools/%s.samtools.somatic.snv.reformated.vcf.gz' % (each,dirtmp,each,each) for each in qc_list]
			snv_stat = ['%s %s/%s/somatools/%s.samtools.somatic.snv.reformated.vcf.gz.stat.txt' % (each,dirtmp,each,each) for each in qc_list]
			indel_stat = ['%s %s/%s/somatools/%s.samtools.somatic.indel.reformated.vcf.gz.stat.txt' % (each,dirtmp,each,each) for each in qc_list]
		elif snvm=='varScan' or indelm =='varScan':
			#order = ['order somaticvarscanannovar_%s before spectrum' % each for each in sample_list]
			order = ['order spectrum after somaticvarscanannovar_%s' % each for each in sample_list]
			list = ['%s %s/%s/varScan/%s.varScan.somatic.snv.reformated.vcf.gz' % (each,dirtmp,each,each) for each in qc_list]
			snv_stat = ['%s %s/%s/varScan/%s.varScan.somatic.snv.reformated.vcf.gz.stat.txt' % (each,dirtmp,each,each) for each in qc_list]
			indel_stat = ['%s %s/%s/varScan/%s.varScan.somatic.indel.reformated.vcf.gz.stat.txt' % (each,dirtmp,each,each) for each in qc_list]
		elif snvm=='muTect' or indelm =='Strelka':
			#order = ['order somatic_muTect_%s before spectrum' % each for each in sample_list]
			#order += ['order somatic_Strelka_%s before spectrum' % each for each in sample_list]
			order = ['order spectrum after somatic_muTect_catannovar_%s' % each for each in sample_list]
			order += ['order spectrum after somatic_Strelka_%s' % each for each in sample_list]
			list = ['%s %s/%s/muTect/%s.muTect.somatic.snv.reformated.vcf.gz' % (each,dirtmp,each,each) for each in qc_list]
			snv_stat = ['%s %s/%s/muTect/%s.muTect.somatic.snv.reformated.vcf.gz.stat.txt' % (each,dirtmp,each,each) for each in qc_list]
			indel_stat = ['%s %s/%s/Strelka/%s.Strelka.somatic.indel.reformated.vcf.gz.stat.txt' % (each,dirtmp,each,each) for each in qc_list]
		else:
			print 'Tools for somatic snp and indel are not properly choosed!'
			sys.exit(1)

		open(spectrumdir+'/snv.list','w').write('\n'.join(snv_stat)+'\n')
		open(spectrumdir+'/indel.list','w').write('\n'.join(indel_stat)+'\n')
		cmds.append('python %s/summary_snp.py %s/snv.list %s' % (self.advDir, spectrumdir, spectrumdir))
		cmds.append('python %s/summary_indel.py %s/indel.list %s' % (self.advDir, spectrumdir, spectrumdir))
		cmds.append('rm -f %s/snp_features.xls %s/indel_features.xls' % (spectrumdir, spectrumdir))
		cmds.append('mv -f %s/snp_function.stat.xls %s/somatic_snv.function.stat.xls' % (spectrumdir, spectrumdir))
		cmds.append('mv -f %s/indel_function.stat.xls %s/somatic_indel.function.stat.xls' % (spectrumdir, spectrumdir))
		return ' && \\\n'.join(cmds)+'\n',order

