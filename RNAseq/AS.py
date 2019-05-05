import os,sys
class AS:
	def __init__(self,CompareName,BamDir,ASDir,databases,softwares):
		self.CompareName = CompareName
		self.BamDir = BamDir
		self.ASDir = ASDir
		self.databases = databases
		self.softwares=softwares
	def run_MATS(self,soft='rMATS', pair='paired', length='150',libType='fr-unstranded'):
		if soft=='rMATS':
			order='oeder run_MATS_%s after mapping_rnaseq_%s\norder run_MATS_%s after mapping_rnaseq_%s)'%(self.CompareName,self.CompareName.strip().split("_vs_")[0],self.CompareName,self.CompareName.strip().split("_vs_")[1])
			step1='mkdir -p %s/rMATS/%s'%(self.ASDir,self.CompareName)
			step2='\\\n\t'.join([self.softwares['rMATS'],
				'-b1 %s/%s.bam'%(self.ASDir,self.CompareName.strip().split("_vs_")[0]),
				'-b2 %s.bam'%(os.path.join(self.BamDir,self.CompareName.strip().split("_vs_")[1])),
				'-gtf %s.gtf '%(self.databases),
				'-t %s -len %s -a 8 -r1 544 -r2 438 -sd1 512 -sd2 336 -c 0.01 -analysis U -libType %s'%(pair,length,libType),
				'-o %s'%(os.path.join(self.ASDir,'rMATS',self.CompareName))])
			return ' && \\\n'.join([step1,step2]),order
	def rMATS_result(self,soft='rMATS'):
		if soft=='rMATS':
			order='oeder rMATS_result%s after run_MATS_%s'%(self.CompareName,self.CompareName)
			AS_type=['SE','RI','MXE','A5SS','A3SS']
			AS_list=[]
			for each_type in AS_type:
				step1='\\\n\t'.join([self.softwares['rMATSresult'],
					'--infile %s/rMATS/%s/MATS_output//%s.MATS.ReadsOnTargetAndJunctionCounts.txt'%(self.ASDir,self.CompareName,each_type),
					'--cutoff 0.01 ',
					'--outfile %s/rMATS/%s/%s.significant.xls'%(self.ASDir,self.CompareName,each_type)])
				step2='cp %s/rMATS/%s/%s.significant.xls %s/rMATS/%s/%s.plot.txt'%(self.ASDir,self.CompareName,each_type,self.ASDir,self.CompareName,each_type)
				step3="number=$(wc -l %s/rMATS/%s/%s.significant.xls|cut -d ' ' -f 1)\nif [[ $number -lt 2 ]]; then\n "%(self.ASDir,self.CompareName,each_type)
				step4='head -n 6 %s/rMATS/%s/%s.MATS.ReadsOnTargetAndJunctionCounts.txt > %s/rMATS/%s/%s.plot.txt\nfi'%(self.ASDir,self.CompareName,each_type,self.ASDir,self.CompareName,each_type)
				step5="sed -i 's/NA/0/g' %s/rMATS/%s/%s.plot.txt"%(self.ASDir,self.CompareName,each_type)
				step6='''sed -i 's/"//g' %s/rMATS/%s/%s.MATS.ReadsOnTargetAndJunctionCounts.txt'''%(self.ASDir,self.CompareName,each_type)
				step7='''sed  -i 's/"//g' %s/rMATS/%s/%s.significant.xls'''%(self.ASDir,self.CompareName,each_type)
				step8='\\\n\t'.join([self.softwares['rmats2sashimiplot'],
					'-b1 %s.bam'%(os.path.join(self.BamDir,self.CompareName.strip().split("_vs_")[0])),
					'-b2 %s.bam'%(os.path.join(self.BamDir,self.CompareName.strip().split("_vs_")[1])),
					'-t %s'%(each_type),
					'-e %s/rMATS/%s/%s.plot.txt'%(self.ASDir,self.CompareName,each_type),
					'--l1 %s'%(self.CompareName.strip().split("_vs_")[0]),
					'--l2 %s'%(self.CompareName.strip().split("_vs_")[1]),
					'-o %s/rMATS/%s/%s'%(self.ASDir,self.CompareName,each_type)])
				AS_list.append(' && \\\n'.join([step1,step2,step3,step4,step5,step6,step7,step8]))
			return ' && \\\n'.join(AS_list),order
