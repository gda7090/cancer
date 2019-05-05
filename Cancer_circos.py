
# update      :		add sample name
#######################################################################
										  
#!/usr/bin/python
import os,sys
import argparse
import re,math

parser = argparse.ArgumentParser(description="Circos module for human cancer and mouse cancer.")
parser.add_argument('-b','--bam',help="Bam file",required=True)
parser.add_argument('-p','--sosnp',help="somatic snp vcf/bed",default=None)
parser.add_argument('-i','--soindel',help="somatic indel vcf/bed",default=None)
parser.add_argument('-c','--socnv',help="somatic cnv profile",default=None)
parser.add_argument('-v','--sosv',help="somatic sv gff",default=None)
parser.add_argument('-s','--sample',help="sample name",required=True)
parser.add_argument('-m','--mouse',help="Mouse circos",action='store_true')
parser.add_argument('-w','--window',help="window size, default 500000",type=int,default=500000)
parser.add_argument('-r','--rmchr',help="remove chromosome from circos, default: None, example: X or X,Y",default=None)
parser.add_argument('-o','--outdir',help="Output directory",required=True)
argv = parser.parse_args()


chrLen = {'chrs':[str(x) for x in range(1,23)]+['X','Y'],
	'len':{'1':249250621, '2':243199373, '3':198022430, '4':191154276, '5':180915260,
	  '6':171115067, '7':159138663, '8':146364022, '9':141213431, '10':135534747,
	  '11':135006516, '12':133851895, '13':115169878, '14':107349540, '15':102531392, 
	  '16':90354753, '17':81195210, '18':78077248, '19':59128983, '20':63025520, 
	  '21':48129895, '22':51304566, 'X':155270560, 'Y':59373566}}

mouse_chrLen = {'chrs':[str(x) for x in range(1,20)]+['X','Y'],
	"len":{'1':195471971, '2':182113224, '3':160039680, '4':156508116, '5':151834684,
      '6':149736546, '7':145441459, '8':129401213, '9':124595110, '10':130694993, '11':122082543,
      '12':120129022, '13':120421639, '14':124902244, '15':104043685, '16':98207768,
      '17':94987271, '18':90702639, '19':61431566, 'X':171031299, 'Y':91744698}}


def safe_open(file_name,mode='r'):
	try:
		if not file_name.endswith('.gz'):
			return open(file_name,mode)
		else:
			import gzip
			return gzip.open(file_name,mode)
	except IOError:
		print file_name + ' do not exist!'

def getInterval(chrLen,window=500000):
#	chrs = [str(x) for x in range(1,23)]+['X','Y']
	chrs = chrLen['chrs']
	intervals = []
	for each in chrs:
		interval = [(each,str(e+1),str(min(chrLen['len'][each],e+window))) for e in range(0,chrLen['len'][each],window)]
		intervals += interval      # for example, 1,50000,1000000
	return intervals


def catVcf(vcflist, cvcf):
	cvcf_tmp = safe_open(cvcf,'w')
	for line in safe_open(vcflist[0]):
		if line.startswith('#'):
			cvcf_tmp.write(line)
			continue
		array = line.strip().split('\t')
		array[0] = array[0].replace('chr','')
		cvcf_tmp.write('\t'.join(array)+'\n')
	if len(vcflist) == 1:
		cvcf_tmp.close()
		return True
	for eachvcf in vcflist[1:]:
		for line in safe_open(eachvcf):
			if line.startswith('#'):
				continue
			array = line.strip().split('\t')
			array[0] = array[0].replace('chr','')
			cvcf_tmp.write('\t'.join(array)+'\n')
	cvcf_tmp.close()
	return True

def bam2circos(bamfile,window,out, ifMouse=False):
	if ifMouse:
		shell = 'Circos/IGVTools/igvtools count -w %d %s %s.wig %s'%(window, bamfile, out, 'mm10')
	else:
		shell = 'Circos/IGVTools/igvtools count -w %d %s %s.wig %s'%(window, bamfile, out, 'novo37')
	print shell
	assert not os.system(shell)
	out_tmp = safe_open(out,'w')
	for line in safe_open(out+'.wig'):
		if line.startswith('track'):
			continue
		elif line.startswith('variableStep'):
			chrs = re.search('chrom=chr(\w+) ',line.strip())
			if chrs:
				if argv.mouse:
					chr = 'mm'+chrs.group(1)
				else:
					chr = 'hs'+chrs.group(1)
		else:
			array = line.strip().split()
			out_tmp.write('%s\t%s\t%s\t%f\n'%(chr,array[0],array[0],math.log(float(array[1])+1)))
	out_tmp.close()
	return True


def vcf2circos(vcf,window_bed,out):
	shell = 'bedtools intersect -a %s -b %s -c | awk \'$4>0{print "hs"$1"\\t"$2"\\t"$3"\\t"$4}\' > %s' % (window_bed, vcf, out)
	if argv.mouse:
		shell = 'bedtools intersect -a %s -b %s -c | awk \'$4>0{print "mm"$1"\\t"$2"\\t"$3"\\t"$4}\' > %s' % (window_bed, vcf, out)
		
	print shell
	assert not os.system(shell)

def cnv2circos(cnv,out,step=20):
	if not os.path.exists(cnv):
		cnv += '.gz'
	out_tmp = safe_open(out,'w')
	i = 0
	for line in safe_open(cnv):
		if i == 0:
			i += 1
			continue
		i += 1
		if not i % step == 2:
			continue
		# chr start ratio_raw ratio_norm cnv
		array = line.strip().split('\t')
		# color
		color = 'vvlgrey'
		cn = int(array[4])
		if cn >= 0 and cn<2:
			color = 'loss'
		elif cn == 2:
			color = 'normal'
		elif cn > 2:
			color = 'gain'
		# ratio/cn
		ratio = float(array[3])
		if array[3] == '-1':
			ratio = -0.1
			color = 'vvlgrey'
		else:
			ratio = ratio * 2
		if color == 'vvlgrey':
			continue
		if array[0].replace('chr','') in ['X','Y']:
			ratio = 2
			color = 'normal'
		if argv.mouse:
			out_tmp.write('mm%s\t%s\t%s\t%f\tcolor=%s\n'%(array[0].replace('chr',''), array[1], array[1], ratio, color))
		else:
			out_tmp.write('hs%s\t%s\t%s\t%f\tcolor=%s\n'%(array[0].replace('chr',''), array[1], array[1], ratio, color))
	out_tmp.close()
	return True

def sv2circos(svann,out):
	out_tmp = safe_open(out,'w')
	i = 1
	svid = 51   #?
	svtype = 53  #?
	svs = {}
	for line in safe_open(svann):
		array = line.strip().split('\t')
		if i == 1:
			i += 1
			svid = array.index('SVID')
			svtype = array.index('SVType')
			continue
		if array[svtype] == 'breakpoint':
			continue
		if argv.mouse:
			info1 = ['mm'+array[0].replace('chr',''),array[1],array[1]]
			info2 = ['mm'+array[0].replace('chr',''),array[2],array[2]]
		else:
			info1 = ['hs'+array[0].replace('chr',''),array[1],array[1]]
			info2 = ['hs'+array[0].replace('chr',''),array[2],array[2]]
		
		if array[svid] not in svs:
			svs[array[svid]] = {'type':array[svtype].lower(),'info':[]}
		if info1 == info2:
			svs[array[svid]]['info'] += info1
		else:
			svs[array[svid]]['info'] += info1+info2
		
	for each in svs:
		if not len(svs[each]['info']) == 6:
			continue
		if svs[each]['type'] == 'translocation':
			if svs[each]['info'][0] == svs[each]['info'][3]:
				out_tmp.write('\t'.join(svs[each]['info'])+'\tcolor=intra%s\n' % svs[each]['type'])  #color=inter?
			else:
				out_tmp.write('\t'.join(svs[each]['info'])+'\tcolor=inter%s\n' % svs[each]['type']) #color=intra ?
		else:
			out_tmp.write('\t'.join(svs[each]['info'])+'\tcolor=%s\n' % svs[each]['type'])
	out_tmp.close()
	return True


def create_dir(dir):
	if not os.path.exists(dir):
		assert not os.system('mkdir -p ' + dir)

outdir = argv.outdir
infodir = os.path.join(outdir,'infos')
create_dir(outdir)
create_dir(infodir)
bed_file = os.path.join(infodir,'window_%d.bed'%(argv.window))
bam_info = os.path.join(infodir,'%s.baminfo'%(argv.sample))
vcf_info = os.path.join(infodir,'%s.vcfinfo'%(argv.sample))
cnv_info = os.path.join(infodir,'%s.cnvinfo'%(argv.sample))
sv_info = os.path.join(infodir,'%s.svinfo'%(argv.sample))

if argv.mouse:
	intervals = getInterval(mouse_chrLen,window=argv.window)
	safe_open(bed_file,'w').write('\n'.join(['\t'.join(each) for each in intervals])+'\n')
else:
	intervals = getInterval(chrLen,window=argv.window)
	safe_open(bed_file,'w').write('\n'.join(['\t'.join(each) for each in intervals])+'\n')

info_files = []

ideogram = '''# The <ideogram> block defines the position, size, labels and other
# properties of the segments on which data are drawn. These segments
# are usually chromosomes, but can be any integer axis.
<ideogram>
# Spacing between ideograms. Suffix "r" denotes a relative value. It
# is relative to circle circumference (e.g. space is 0.5% of
# circumference).
<spacing>
default = 0.005r
# You can increase the spacing between specific ideograms.
#<pairwise hsY;hs1>
#spacing = 3r
#</pairwise>
</spacing>
<<include ideogram.position.conf>>
<<include ideogram.label.conf>>
<<include bands.conf>>
</ideogram>
'''

ideogram_label = '''show_label       = yes
label_font       = default
#label_radius     = dims(image,radius) - 60p
label_radius     = 1.075r
label_with_tag   = yes
label_size       = 35
label_parallel   = yes
label_case       = upper
#label_format     = eval(sprintf("chr%s",var(label)))
label_format     = eval(sprintf("%s",var(label)))
'''

ideogram_position='''# Ideogram position, thickness and fill. 
radius           = 0.80r
thickness        = 100p
fill             = yes
fill_color       = black
stroke_thickness = 2
stroke_color     = black
'''

bands='''show_bands            = yes
fill_bands            = yes
band_stroke_thickness = 2
band_stroke_color     = white
band_transparency     = 0
'''

ticks='''show_ticks        = yes
show_tick_labels  = yes
<ticks>
radius         = 1r
color          = black
thickness      = 2p
multiplier     = 1e-6
format         = %d

<tick>
spacing        = 10u
size           = 15p
</tick>

<tick>
spacing        = 50u
size           = 25p
show_label     = yes
label_size     = 20p
label_offset   = 10p
format         = %d
</tick>

</ticks>
'''

safe_open(os.path.join(outdir,'ideogram.conf'),'w').write(ideogram)
safe_open(os.path.join(outdir,'ideogram.label.conf'),'w').write(ideogram_label)
safe_open(os.path.join(outdir,'ideogram.position.conf'),'w').write(ideogram_position)
safe_open(os.path.join(outdir,'bands.conf'),'w').write(bands)
safe_open(os.path.join(outdir,'ticks.conf'),'w').write(ticks)

##config
if argv.mouse:
	config = 'karyotype = data/karyotype/karyotype.mouse.mm10.txt'+'\n'
else:
	config = 'karyotype = data/karyotype/karyotype.human.txt'+'\n'

if argv.rmchr:
	rmchrs = argv.rmchr.split(',')
	chromosomes=[]
	mouse_chrs=[str(x) for x in range(1,20)]+['X','Y']
	human_chrs=[str(x) for x in range(1,23)]+['X','Y']
	if argv.mouse:
		config+='chromosomes_display_default=no'+'\n'
		for each in rmchrs:
			if each in mouse_chrs:
				mouse_chrs.remove(each)
			else:
				continue
		config+='chromosomes=%s'%(';'.join(['mm%s'%each.replace('chr','') for each in mouse_chrs]))
		#config += 'chromosomes = %s'%(';'.join(['-mm%s'%each.replace('chr','') for each in rmchrs]))
	else:	
		config+='chromosomes_display_default=no'+'\n'
		for each in rmchrs:
			if each in human_chrs:
				human_chrs.remove(each)
			else:
				continue
		config+='chromosomes=%s'%(';'.join(['hs%s'%each.replace('chr','') for each in human_chrs]))
		#config += 'chromosomes = %s'%(';'.join(['-hs%s'%each.replace('chr','') for each in rmchrs]))
	

config += '''
chromosomes_units = 1000000
<colors>
loss = 0,0,205
gain = 205,38,38
normal = 162,205,90
deletion = gpos75
inversion = chr14
insertion = dgreen
intratranslocation = orange
intertranslocation = chr6
ctx = orange
del = gpos75
ins = dgreen
inv = chr14
itx = chr6
</colors>
<<include ideogram.conf>>
<<include ticks.conf>>
#The remaining content is standard and required.
<image>
<<include etc/image.conf>>
</image>
<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
data_out_of_range* = trim

<plots>
'''

r1 = 0.98
r0=0.90
### bam
bam2circos(argv.bam, argv.window, bam_info, argv.mouse)
info_files.append(bam_info)
config += '''
<plot>
type             = histogram
r1               = %0.2fr
r0               = %0.2fr
stroke_thickness = 0
file             = %s
orientation      = out
extend_bin       = no
fill_color       = blue
color       = blue
</plot>
'''%(r1,r0,bam_info)

## vcf
if argv.sosnp or argv.soindel:
	vcflist = []
	vcf_cmb = os.path.join(infodir,'%s.snv_indel.vcf'%(argv.sample))
	if argv.sosnp:
		vcflist.append(argv.sosnp)
	if argv.soindel:
		vcflist.append(argv.soindel)
	catVcf(vcflist,vcf_cmb)    #return the result of combining snp and indel.
	vcf2circos(vcf_cmb, bed_file, vcf_info)
	info_files.append(vcf_info)
	r1 = r0 - 0.01
	r0 = r1 - 0.1
	config += '''
<plot>
type             = scatter
r1               = %0.2fr
r0               = %0.2fr
stroke_thickness = 1
file             = %s
orientation      = out
fill_color       = dgreen
stroke_color     = dgreen
glyph            = circle
glyph_size       = 3
<backgrounds>
<background>
color            = vvlgrey
y1               = 1r
</background>
</backgrounds>
<axes>
<axis>
color            = lgrey
thickness        = 1
spacing          = 0.1r
y1               = 1r
</axis>
</axes>
</plot>
'''%(r1,r0,vcf_info)

## cnv
if argv.socnv:
	cnv2circos(argv.socnv, cnv_info, step=20)
	info_files.append(cnv_info)
	r1 = r0 - 0.01
	r0 = r1 - 0.1
	config += '''
<plot>
type         = scatter
r1           = %0.2fr
r0           = %0.2fr
file         = %s
min          = 0
max          = 6
glyph        = circle
glyph_size   = 3
thickness    =1


<backgrounds>
<background>
color        = vvlgrey
y1           = 1r
</background>
</backgrounds>
</plot>
'''%(r1,r0,cnv_info)
config += '</plots>'
## sv
if argv.sosv:
	sv2circos(argv.sosv, sv_info)
	info_files.append(sv_info)
	r1 = r0 - 0.01
	r0 = r1 - 0.02
	config += '''
<links>
z      = 50
radius = %0.2fr
crest  = 1
bezier_radius        = 0r
bezier_radius_purity = %0.2f
<link>
thickness = 3
file      = %s
</link>
</links>
''' % (r1,r0,sv_info)
safe_open(os.path.join(outdir,argv.sample+'.circos.conf'),'w').write(config)
outfile = os.path.join(outdir,argv.sample+'.circos')
os.chdir(outdir)
assert not os.system('circos -conf %s.conf -outputfile %s'%(outfile, outfile))

assert not os.system('sed -i \'s/<\/svg>/<g><text id="TextElement" x="1400" y="100" style="font-family:Verdana;font-size:25pt"> %s<\/text><\/g><\/svg>/\' %s.circos.svg' %(argv.sample,argv.sample))
assert not os.system('convert %s.circos.svg %s.circos.png'%(argv.sample, argv.sample))
