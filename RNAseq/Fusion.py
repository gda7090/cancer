import os,sys

class Fusion:

	def __init__(self,sampleID,cleanDir,fusionDir,softwares,databases):
		self.sampleID = sampleID
		self.cleanDir = cleanDir
		self.fusionDir = fusionDir
		self.softwares = softwares
		self.databases = databases
		self.fq1_gz = os.path.join(self.cleanDir,'_'.join([self.sampleID,'1.clean.fq.gz']))
		self.fq2_gz = os.path.join(self.cleanDir,'_'.join([self.sampleID,'2.clean.fq.gz']))

	def call_fusion(self,software='starfusion'):
		order = 'order call_fusion_%s after mapping_rnaseq_%s' % (self.sampleID,self.sampleID)
		if software == 'starfusion':
			step1 = '\\\n\t'.join(['STAR-Fusion/STAR-Fusion',
				'--genome_lib_dir %s'%(self.databases),
				'--left_fq %s'%(self.fq1_gz),
				'--right_fq %s'%(self.fq2_gz),
				'--output_dir %s/starfusion/%s'%(self.fusionDir,self.sampleID)])
			step2 = 'cut -f 1 %s/starfusion/%s/star-fusion.fusion_candidates.final.abridged |grep -v "#" > %s/starfusion/%s/%s.fusionlist'%(self.fusionDir,self.sampleID,self.fusionDir,self.sampleID,self.sampleID)
			step3 = '\\\n\t'.join(['FusionInspector_v0.5.0_FULL/FusionInspector',
				'--fusions %s/starfusion/%s/%s.fusionlist'%(self.fusionDir,self.sampleID,self.sampleID),
				'--genome_lib %s'%(self.databases),
				'--left_fq %s'%(self.fq1_gz),
				'--right_fq %s'%(self.fq2_gz),
				'--out_dir %s/starfusion/%s'%(self.fusionDir,self.sampleID),
				'--out_prefix %s --prep_for_IGV'%(self.sampleID)])
			step4 = 'mkdir -p %s/starfusion/%s/final_fusion_genes'%(self.fusionDir,self.sampleID)
			step5 = 'mv  %s/starfusion/%s/star-fusion.predict.intermediates_dir %s/starfusion/%s/star-fusion.filter.intermediates_dir %s/starfusion/%s/final_fusion_genes'%(self.fusionDir,self.sampleID,self.fusionDir,self.sampleID,self.fusionDir,self.sampleID)
			step6 = 'sed  s/^#//g %s/starfusion/%s/%s.fusion_predictions.final.abridged | cut -f  1-8 > %s/starfusion/%s/inal_fusion_genes/%s.fusion_predictions.final.abridged'%(self.fusionDir,self.sampleID,self.sampleID,self.fusionDir,self.sampleID,self.sampleID)
			step7 = 'mkdir %s/starfusion/%s/IGV'%(self.fusionDir,self.sampleID)
			step8 = 'cp {self.fusionDir}/starfusion/{self.sampleID}/{self.sampleID}.fa {self.fusionDir}/starfusion/{self.sampleID}/{self.sampleID}.gtf {self.fusionDir}/starfusion/{self.sampleID}/{self.sampleID}.bed {self.fusionDir}/starfusion/{self.sampleID}/{self.sampleID}.junction_reads.bam {self.fusionDir}/starfusion/{self.sampleID}/{self.sampleID}.junction.bam.bai {self.fusionDir}/starfusion/{self.sampleID}/{self.sampleID}.spanning_reads.bam {self.fusionDir}/starfusion/{self.sampleID}/{self.sampleID}.spanning_reads.bam  {self.fusionDir}/starfusion/{self.sampleID}/{self.sampleID}.spanning_reads.bam.bai {self.fusionDir}/starfusion/{self.sampleID}/IGV'.format(**locals())
			return ' && \\\n'.join([step1,step2,step3,step4,step5,step6,step7,step8]),order
	def fusion_stat(self,software='starfusion'):
		order = 'order fusion_stat_%s after call_fusion_%s' % (self.sampleID,self.sampleID)
		if software == 'starfusion':
			step1='mkdir -p  {self.fusionDir}/starfusion/{self.sampleID}/circos'.format(**locals())
			step2='\\\n\t'.join(['bin/Rcicos_starFusion_prepare',
				'--fusion {self.fusionDir}/starfusion/{self.sampleID}/{self.sampleID}.fusion_predictions.final.abridged'.format(**locals()),
				'--outdir {self.fusionDir}/starfusion/{self.sampleID}/circos'.format(**locals())])
			step3='sed -i s/rM/rMT/g {self.fusionDir}/starfusion/{self.sampleID}/circos/link.xls'.format(**locals())
			step4='sed -i s/rM/rMT/g {self.fusionDir}/starfusion/{self.sampleID}/circos/name.xls'.format(**locals())
			step5='\\\n\t'.join(['R/RCircos.R',
				'--cytoband {self.databases}.cytoband'.format(**locals()),
				'--link {self.fusionDir}/starfusion/{self.sampleID}/circos/link.xls'.format(**locals()),
				'--name {self.fusionDir}/starfusion/{self.sampleID}/circos/name.xls',
				'--prefix {self.fusionDir}/starfusion/{self.sampleID}/circos/{self.sampleID}.circos'])
			step6='mkdir -p {self.fusionDir}/starfusion/{self.sampleID}/Annotation'.format(**locals())
			step7='cp {self.fusionDir}/{self.sampleID}/{self.sampleID}.fa {self.fusionDir}/starfusion/{self.sampleID}/Annotation'.format(**locals())
			step8='/PUBLIC/source/RNA/med/Pipeline/medpipline1.2/bin/fusiongene_extact_id_seq  {self.fusionDir}/{self.sampleID}/{self.sampleID}.fa {self.fusionDir}/starfusion/{self.sampleID}/Annotation/fusiongene.seq'.format(**locals())
			step9="sed -i '1iFusion_name   Fusion_TranscriptSeq' {self.fusionDir}/starfusion/{self.sampleID}/Annotation/fusiongene.seq".format(**locals())
			step10='export PATH=software/jre1.8.0_51/bin:$PATH\nDatabase/Annotation/interproscan/interproscan-5.17-56.0/interproscan -t n -i   {self.fusionDir}/starfusion/{self.sampleID}/Annotation/{self.sampleID}.fa -f TSV -b {self.fusionDir}/starfusion/{self.sampleID}/Annotation/{self.sampleID} -T {self.fusionDir}/starfusion/{self.sampleID}/Annotation/temp -dp'.format(**locals())
			step11='bin/interproscandeal --infile {self.fusionDir}/starfusion/{self.sampleID}/Annotation/{self.sampleID}.tsv --outfile {self.fusionDir}/starfusion/{self.sampleID}/Annotation/{self.sampleID}.interproscan.xls'.format(**locals())
