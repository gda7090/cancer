import os,sys

if len(sys.argv) < 2:
    exit('Usage: python %s input ouputbed' %sys.argv[0])

def safe_open(fileName,mode):
    try:
        if not fileName.endswith('.gz'):
            return open(fileName,mode)
        else:
            import gzip
            return gzip.open(fileName,mode)
    except IOError:
        print fileName + ' do not exist!'


input = safe_open(sys.argv[1],'r')
outbed = open(sys.argv[2],'w')

maxbed = ''
maxv = 0
for line in input:
    if line.startswith('N'):
        continue
    uline = line.strip().split('\t')
    bedline = '\t'.join([uline[1],uline[2],uline[3]])
    if 'Chr' in bedline:
        bedline = bedline.replace('Chr','')
    if uline[6] > maxv:
        maxv = uline[6]
        maxbed = bedline
outbed.write(maxbed+'\n')
