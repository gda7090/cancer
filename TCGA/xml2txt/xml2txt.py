#!Python2.7.12/bin/python
import xml.dom.minidom,sys,argparse
reload(sys)
sys.setdefaultencoding('utf-8')
parser=argparse.ArgumentParser(usage="xml2txt.py --in infile --out outfile\n",description="change XML 2 txt")
parser.add_argument('-i','--infile',help="the input file",required=True)
parser.add_argument('-o','--outfile',help="the output file",required=True)
argv = parser.parse_args()
infile=argv.infile
outfile=argv.outfile

def getnode(node,name):
    childs=node.childNodes
    num=len(childs)
    if num==1 :
        outfile.write(name+": "+childs[0].data+"\n")
    else :
        for e in childs:
            if e.nodeName!="#text":
                getnode(e,name+"/"+e.nodeName)

outfile=open(argv.outfile,"w")
DOMTree = xml.dom.minidom.parse(infile)
tcgaxml= DOMTree.documentElement
getnode(tcgaxml,tcgaxml.nodeName)
outfile.close()
