import os
class AdvAnalysis:
	def __init__(self,Tsamples,Nsamples,T2pid,T2Nsample,softwares,advdir,mapDir,mutDir,svDir,somaticDir,mutM,svM,cnvM,sosnvM,soindelM,sosvM,socnvM,seqstrag,genome_info):
		self.Tsamples = Tsamples
		self.Nsamples = Nsamples
		self.T2pid = T2pid
		self.T2Nsample = T2Nsample
		self.softwares = softwares
		self.advbin = softwares['advbin']
		self.advdir = advdir
		self.mapDir = mapDir
		self.mutDir = mutDir
		self.svDir = svDir
		self.somaticDir = somaticDir
		self.mutM = mutM
		self.svM = svM
		self.cnvM = cnvM
		self.sosnvM = sosnvM
		self.soindelM = soindelM
		self.socnvM = socnvM
		self.sosvM = sosvM
		self.seqstrag = seqstrag
		self.genome_info = genome_info
		self.refData = genome_info['fasta']
		self.platform = {'WGS':'HiseqX'}.get(seqstrag,'Hiseq2500')
		self.prepareDir = os.path.join(advdir,'00.prepare')
		self.signatureDir = os.path.join(advdir,'Signatures')
		self.driverDir = os.path.join(advdir,'DriverGenes')
		self.musicDir = os.path.join(advdir,'Smg')
		self.gisticDir = os.path.join(advdir,'Gistic')
		self.fusionDir = os.path.join(advdir,'FusionGenes')
		self.evolutionDir = os.path.join(advdir,'Evolution')
		self.absoluteDir = os.path.join(advdir,'Absolute')
		self.cloneDir = os.path.join(advdir,'Clone')
		self.circosDir = os.path.join(advdir,'Circos')
		self.DriverGenesDir = os.path.join(advdir,'DriverGenes')
		self.predisposeDir = os.path.join(advdir,'PredisposeGenes')
		self.mrtDir = os.path.join(advdir,'Mrt')
		self.oncodriveDir = os.path.join(advdir,'OncodriveClust')
		self.targetDir = os.path.join(advdir,'DrugTarget')
		self.resistanceDir = os.path.join(advdir,'DrugResistance')
		self.noncodingDir = os.path.join(advdir,'Noncoding')

	def contruct_sovcf(self,sample,method,type):
		return os.path.join(self.somaticDir,sample,method,'%s.%s.somatic.%s.reformated.vcf.gz'%(sample,method,type))

	def contruct_somaf(self,sample,method,type):
		return os.path.join(self.somaticDir,sample,method,'%s.%s.somatic.%s.maf'%(sample,method,type))

	def contruct_vcf(self,sample,method,type):
		return os.path.join(self.mutDir,sample+'.'+method,'%s.%s.%s.reformated.vcf.gz'%(sample,method,type))

	def contruct_sogff(self,sample,method,type):
		return os.path.join(self.somaticDir,sample,method,'%s.%s.somatic.%s.gff'%(sample,method,type))

	def contruct_socnv(self,sample,method):
		return os.path.join(self.somaticDir,sample,method,'%s.somatic.CNV.txt'%(sample))

	def contruct_gff(self,sample,method,type):
		return os.path.join(self.svDir,sample,method,'%s.%s.%s.gff'%(sample,method,type))

	def contruct_bam(self,sample):
		return os.path.join(self.mapDir,sample,sample+'.final.bam')

	def vcf2maf(self,sample,method,type):
		vcf = self.contruct_sovcf(sample,method,type)
		maf = self.contruct_somaf(sample,method,type)
		step1 = '\\\n\t'.join(['python %s'%(os.path.join(self.advbin,'vcf2maf.py')),
			'-v %s' % vcf,
			'-o %s' % maf,
			'-t %s -n %s' % (sample,self.T2Nsample[sample]),
			'-s %s -p %s' % (self.seqstrag,self.platform),
			'-m %s' % method])
		return step1+'\n'

	def combine_maf(self):
		snv_maf = [self.contruct_somaf(each,self.sosnvM,'snv') for each in self.Tsamples]
		indel_maf = [self.contruct_somaf(each,self.soindelM,'indel') for each in self.Tsamples]
		step1 = [self.vcf2maf(each,self.sosnvM,'snv') for each in self.Tsamples]
		step2 = [self.vcf2maf(each,self.soindelM,'indel') for each in self.Tsamples]
		step3 = 'awk \'NR==1||FNR>1\' %s > %s' % \
			(' '.join(snv_maf+indel_maf), os.path.join(self.prepareDir,'Somatic_mutation.maf'))
		return step3+'\n'
#		return ' && \\\n'.join(step1+step2+[step3])

	def signature_spectrum(self,rank=3):
		order = 'order signature_spectrum after combine_maf'
		step1 = '\\\n\t'.join(['python %s'%(os.path.join(self.advbin,'mutation_signature.py')),
			'-m %s' % os.path.join(self.prepareDir,'Somatic_mutation.maf'),
			'-o %s' % self.signatureDir,
			'-g %s' % self.refData,
			'-r %s -a 4 -b 5 -i 10 -j 12 -s 15' % rank])
#			'-c %s' % os.path.join(self.prepareDir,'coverage_bases.txt')])
		return step1+'\n',order

	def driverGene(self):
		order = 'order driverGene after combine_maf'
		step1 = '\\\n\t'.join(['MutSig/MutSigCV_1.4/run_MutSigCV.sh',
			os.path.join(self.prepareDir,'Somatic_mutation.maf'),
			'MutSig/dbfiles/exome_full192.coverage.txt',
			'MutSig/dbfiles/gene.covariates.txt',
			os.path.join(self.driverDir,'MutSigCV'),
			'MutSig/dbfiles/mutation_type_dictionary_file.txt',
			'driver_mutation/data/chr_files_hg19'])
		step2 = 'cut -f 1,8,9,10,14,15 %s > %s' % \
			(os.path.join(self.driverDir,'MutSigCV.sig_genes.txt'),os.path.join(self.driverDir,'MutSigCV.sig_genes.xls'))
		return ' && \\\n'.join([step1,step2])+'\n',order

	def smg(self):
		order = 'order smg after combine_maf'
		info = ['%s\t%s\t%s\n' % (each, self.contruct_bam(each), self.contruct_bam(self.T2Nsample[each])) for each in self.Tsamples]
		open(os.path.join(self.musicDir,'smg.config'),'w').write(''.join(info))
		open(os.path.join(self.prepareDir,'smg.config'),'w').write(''.join(info))
		step1 = '\\\n\t'.join(['python %s'%(os.path.join(self.advbin,'music_driver.py')),
			'-c %s' % (os.path.join(self.musicDir,'smg.config')),
			'-m %s' % (os.path.join(self.prepareDir,'Somatic_mutation.maf')),
			'-g %s' % self.refData,
			'-o %s' % self.musicDir])
		step2 = '\\\n\t'.join(['python %s'%(os.path.join(self.advbin,'landscape_mutations.py')),
			'-m %s' % (os.path.join(self.prepareDir,'Somatic_mutation.maf')),
			'-o %s' % os.path.join(self.musicDir,'landscape'),
			'-p Smg -i 11',
			'-v %s' % (os.path.join(self.musicDir,'coverage','total_covgs')),
			'-g %s' %(os.path.join(self.musicDir,'smg.tmp'))])
		return ' && \\\n'.join([step1,step2])+'\n',order

	def gistic(self):
		order = 'order gistic after combine_maf'
		info = ['%s\t%s\n'%(each,self.contruct_socnv(each,self.socnvM)) for each in self.Tsamples]
		open(os.path.join(self.gisticDir,'somatic_cnv.list'),'w').write(''.join(info))
		open(os.path.join(self.prepareDir,'somatic_cnv.list'),'w').write(''.join(info))
		if self.socnvM == 'freec':
			markerfile = 'GISTIC/markerfiles/freec.wgs_w250.markersfile.txt'
			if 'WGS' in self.seqstrag:
				markerfile = 'GISTIC/markerfiles/freec.wgs_w1000.markersfile.txt'
			step1 = '\\\n\t'.join(['python %s'%(os.path.join(self.advbin,'gistic_cnv.py')),
				'-i %s' % os.path.join(self.gisticDir,'somatic_cnv.list'),
				'-m %s' % markerfile,
				'-o %s' % self.gisticDir,
				'-c 0 -s 1 -e 2 -p 5 -n 4 -v cnv'])
		elif self.socnvM == 'ExomeCNV':
			step1 = '\\\n\t'.join(['python %s'%(os.path.join(self.advbin,'gistic_cnv.py')),
				'-i %s' % os.path.join(self.gisticDir,'somatic_cnv.list'),
				'-o %s' % self.gisticDir,
				'-c 0 -s 1 -e 2 -p -1 -n 7 -v cnv'])
		return step1+'\n',order

	def fusion_gene(self,type='sv'):
		order = 'order fusion_gene after combine_maf'
		svann = ['%s\t%s\n' % (each,self.contruct_sogff(each,self.sosvM,type).replace('gff','hg19_multianno.xls')) for each in self.Tsamples]
		fusionlist = [os.path.join(self.fusionDir,'tmp',each+'.fusionGene.xls') for each in self.Tsamples]
		open(os.path.join(self.fusionDir,'somatic_svann.list'),'w').write(''.join(svann))
		open(os.path.join(self.prepareDir,'somatic_svann.list'),'w').write(''.join(svann))
		step1_1 = 'cd %s\nmkdir -p %s/tmp\nwhile read id ann\ndo\n' % (self.fusionDir, self.fusionDir)
		step1_2 = 'var_sv_putative.fusion.gene.pl $ann >%s/tmp/$id.fusionGene.txt && \\' % self.fusionDir
		step1_3 = '\npython %s %s/tmp/$id.fusionGene.txt $id >%s/tmp/$id.fusionGene.xls\ndone<somatic_svann.list' % \
			(os.path.join(self.advbin,'fusion_format.py'), self.fusionDir, self.fusionDir)
		step1 = step1_1+step1_2+step1_3
		step2 = 'awk \'NR==1 || FNR>1\' %s >%s' % (' '.join(fusionlist), os.path.join(self.fusionDir,'FusionGenes.SVbased.xls'))
		return ' && \\\n'.join([step1,step2])+'\n',order

	def phylip_evolution(self,pid,Tsamplelist,func='coding'):
		order = 'order phylip_evolution_%s after combine_maf' % pid
		return '\\\n\t'.join(['python %s' % os.path.join(self.advbin,'phylogenetic_tree_phylip.py'),
				'-i %s ' % (os.path.join(self.prepareDir,'Somatic_mutation.maf')),
				'-s %s -p %s -f coding' % (','.join(Tsamplelist), pid),
				'-o %s' % os.path.join(self.evolutionDir,pid)])+'\n',order

	def absolutes(self):
		order = 'order absolutes after combine_maf'
		maf_cnv = ['%s\t%s\t%s\n' %(each,self.contruct_somaf(each,self.sosnvM,'snv'),self.contruct_socnv(each,self.socnvM)) for each in self.Tsamples]
		options = ''
		if self.socnvM == 'freec':
			options = '-c 0 -s 1 -e 2 -p 5 -n 4 -v cnv -r all'
		elif self.socnvM == 'ExomeCNV':
			options = '-c 0 -s 1 -e 2 -p 5 -n 7 -v cnv -r all'
		open(os.path.join(self.absoluteDir,'maf_cnv.list'),'w').write(''.join(maf_cnv))
		return '\\\n\t'.join(['python %s' % (os.path.join(self.advbin,'Absolutes_analysis.py')),
			'-i %s' % os.path.join(self.absoluteDir,'maf_cnv.list'),
			'-o %s' % self.absoluteDir,
			options])+'\n',order

	def expands(self):
		order = 'order expands after combine_maf'
		maf_cnv = ['%s\t%s\t%s\n' %(each,self.contruct_somaf(each,self.sosnvM,'snv'),self.contruct_socnv(each,self.socnvM)) for each in self.Tsamples]
		options = ''
		if self.socnvM == 'freec':
			options = '-c 0 -s 1 -e 2 -n 4 -v cnv -r gene'
		elif self.socnvM == 'ExomeCNV':
			options = '-c 0 -s 1 -e 2 -n 7 -v cnv -r gene'
		return '\\\n\t'.join(['#python %s' % (os.path.join(self.advbin,'Expands_analysis.py')),
			'-i %s' % os.path.join(self.cloneDir,'expands','maf_cnv.list'),
			'-o %s' % os.path.join(self.cloneDir,'expand'),
			options])+'\n',order

	def sciclones(self):
		order = 'order sciclones after combine_maf'
		maf_cnv = ['%s\t%s\t%s\n' %(each,self.contruct_somaf(each,self.sosnvM,'snv'),self.contruct_socnv(each,self.socnvM)) for each in self.Tsamples]
		options = ''
		if self.socnvM == 'freec':
			options = '-c 0 -s 1 -e 2 -n 4 -v cnv -r gene'
		elif self.socnvM == 'ExomeCNV':
			options = '-c 0 -s 1 -e 2 -n 7 -v cnv -r gene'
		return '\\\n\t'.join(['python %s' % (os.path.join(self.advbin,'Sciclones_analysis.py')),
			'-i %s' % os.path.join(self.cloneDir,'sciclone','maf_cnv.list'),
			'-o %s' % os.path.join(self.cloneDir,'sciclone'),
			options])+'\n',order
		
	def pyclones(self):
		order = 'order pyclones after combine_maf'
		maf_cnv = ['%s\t%s\t%s\n' %(each,self.contruct_somaf(each,self.sosnvM,'snv'),self.contruct_socnv(each,self.socnvM)) for each in self.Tsamples]
		options = ''
		if self.socnvM == 'freec':
			options = '-c 0 -s 1 -e 2 -n 4 -v cnv -r gene -g 7'
		elif self.socnvM == 'ExomeCNV':
			options = '-c 0 -s 1 -e 2 -n 7 -v cnv -r gene'
		return '\\\n\t'.join(['python %s' % (os.path.join(self.advbin,'Pyclones_analysis.py')),
			'-i %s' % os.path.join(self.cloneDir,'pyclone','maf_cnv.list'),
			'-o %s' % os.path.join(self.cloneDir,'pyclone'),
			options])+'\n',order


	def predispose_filter(self):
		snp_list = ['%s\t%s\n' % (each, self.contruct_vcf(each,self.mutM,'snp')) for each in self.Nsamples]
		indel_list = ['%s\t%s\n' % (each, self.contruct_vcf(each,self.mutM,'indel')) for each in self.Nsamples]
		open(os.path.join(self.predisposeDir,'snplist'),'w').write(''.join(snp_list))
		open(os.path.join(self.predisposeDir,'indellist'),'w').write(''.join(indel_list))
		script='''cd %s\nmkdir -p predispose
while read Nid Nvcf
do
\tpython %s \\
\t\t-c $Nvcf \\
\t\t-o predispose/$Nid.predispose_genes.snvs.xls -n $Nid \\
\t\t-x SIFT_score,SIFT_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,MutationTaster_score,MutationTaster_pred,caddgt10 \\
\t\t-b AnnotationDB/CGC/CGC.CancerType \\
\t\t-m %s -f avsnp142
done<%s/snplist
while read Nid Nvcf
do
\tpython %s \\
\t\t-c $Nvcf \\
\t\t-o predispose/$Nid.predispose_genes.indels.xls -n $Nid \\
\t\t-x SIFT_score,SIFT_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,MutationTaster_score,MutationTaster_pred,caddindel \\
\t\t-b AnnotationDB/CGC/CGC.CancerType \\
\t\t-m %s -f avsnp142
done<%s/indellist
awk 'NR==1 || FNR>1' predispose/*.predispose_genes.snvs.xls predispose/*.predispose_genes.indels.xls > %s/Predispose_genes.snvs_indels.xls
''' % (self.predisposeDir,os.path.join(self.advbin,'predispose_genes_filter.py'),self.mutM,self.predisposeDir,os.path.join(self.advbin,'predispose_genes_filter.py'),self.mutM,self.predisposeDir,self.predisposeDir)
		return script

	def drivergenes_filter(self):
		sosnv_list = ['%s\t%s\n' % (each, self.contruct_sovcf(each,self.sosnvM,'snv')) for each in self.Tsamples]
		soindel_list = ['%s\t%s\n' % (each, self.contruct_sovcf(each,self.soindelM,'indel')) for each in self.Tsamples]
		open(os.path.join(self.DriverGenesDir,'sosnvlist'),'w').write(''.join(sosnv_list))
		open(os.path.join(self.DriverGenesDir,'soindellist'),'w').write(''.join(soindel_list))
		script='''cd %s\nmkdir -p drivergene
while read Tid Tvcf
do
\tpython %s -v $Tvcf \\
\t\t-o drivergene/$Tid.DriverGenes.snvs.xls -t $Tid \\
\t\t-x SIFT_score,SIFT_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,MutationTaster_score,MutationTaster_pred,caddgt10 \\
\t\t-b AnnotationDB/CGC/CancerGenes.db \\
\t\t-m %s -f avsnp142
done<%s/sosnvlist
while read Tid Tvcf
do
\tpython %s -v $Tvcf \\
\t\t-o drivergene/$Tid.DriverGenes.indels.xls -t $Tid \\
\t\t-x SIFT_score,SIFT_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,MutationTaster_score,MutationTaster_pred,caddindel \\
\t\t-b AnnotationDB/CGC/CancerGenes.db \\
\t\t-m %s -f avsnp142
done<%s/soindellist
awk 'NR==1 || FNR>1' drivergene/*.DriverGenes.snvs.xls drivergene/*.DriverGenes.indels.xls > %s/DriverGenes.snvs_indels.xls
''' % (self.DriverGenesDir,os.path.join(self.advbin,'driver_genes_filter.py'),self.sosnvM,self.DriverGenesDir,os.path.join(self.advbin,'driver_genes_filter.py'),self.soindelM,self.DriverGenesDir,self.DriverGenesDir)
		return script

	def circos(self,sample,sv=False):
		order = 'order circos_%s after combine_maf'%sample
		if sv:
			return '\\\n\t'.join(['python %s' % os.path.join(self.advbin,'cancer_circos.py'),
				'-b %s' % self.contruct_bam(sample),
				'-p %s' % self.contruct_sovcf(sample,self.sosnvM,'snv'),
				'-i %s' % self.contruct_sovcf(sample,self.soindelM,'indel'),
				'-c %s' % self.contruct_socnv(sample,self.socnvM).replace('somatic.CNV.txt','final.mpileup.gz_ratio.txt'),
				'-v %s' % self.contruct_sogff(sample,self.sosvM,'sv').replace('gff','hg19_multianno.xls'),
				'-o %s' % self.circosDir,
				'-s %s' % sample])+'\n',order
		else:
			return '\\\n\t'.join(['python %s' % os.path.join(self.advbin,'cancer_circos.py'),
				'-b %s' % self.contruct_bam(sample),
				'-p %s' % self.contruct_sovcf(sample,self.sosnvM,'snv'),
				'-i %s' % self.contruct_sovcf(sample,self.soindelM,'indel'),
				'-c %s' % self.contruct_socnv(sample,self.socnvM).replace('somatic.CNV.txt','final.mpileup.gz_ratio.txt'),
				'-o %s' % self.circosDir,
				'-s %s' % sample])+'\n',order

	def mrt_music(self):
		order = 'order mrt_music after smg'
		return ' \\\n\t'.join(['python %s' % os.path.join(self.advbin,'music_mrt.py'),
			'-m %s' % os.path.join(self.prepareDir,'Somatic_mutation.maf'),
			'-o %s' % self.mrtDir,
			'-s %s' % os.path.join(self.musicDir,'Significantly_Mutated_Genes.xls')])+'\n',order

	def oncodriveclust(self):
		order = 'order oncodriveclust after smg'
		return ' \\\n\t'.join(['python %s' % os.path.join(self.advbin,'oncodriveclust.py'),
			'-m %s' % os.path.join(self.prepareDir,'Somatic_mutation.maf'),
			'-o %s' % self.oncodriveDir])+'\n',order

	def drug_target(self):
		order = 'order drug_target after combine_maf'
		return ' \\\n\t'.join(['python %s' % os.path.join(self.advbin,'target_drug.py'),
			'-m %s' % os.path.join(self.prepareDir,'Somatic_mutation.maf'),
			'-o %s' % self.targetDir])+'\n',order

	def drug_resistance(self):
		order = 'order drug_resistance after combine_maf'
		return ' \\\n\t'.join(['python %s' % os.path.join(self.advbin,'DrugResistance.py'),
			'-m %s' % os.path.join(self.prepareDir,'Somatic_mutation.maf'),
			'-o %s' % self.resistanceDir,
			'-r %s' % 'Cancer/Drug/database/Resistance_drug_database'])+'\n',order


	def noncoding_filter(self):
		order = 'order noncoding_filter after combine_maf'
		return ' \\\n\t'.join(['python %s' % os.path.join(self.advbin,'noncoding_filter.py'),
			'-i %s' % os.path.join(self.prepareDir,'Somatic_mutation.maf'),
			'-o %s' % self.noncodingDir])+'\n',order

		self.noncodingDir = os.path.join(advdir,'Noncoding')

		
