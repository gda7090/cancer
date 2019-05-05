#!Python2.7.12/bin/python
# -*- coding: UTF-8 -*-
import xml.dom.minidom,sys,argparse
parser=argparse.ArgumentParser(usage="TCGAXML.py --in infile --out outfile\n",description="change XML 2 txt")
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
	pass
    pass
pass

outfile=open(outfile,"w")
DOMTree = xml.dom.minidom.parse(infile)
tcgaxml= DOMTree.documentElement
getnode(tcgaxml,tcgaxml.nodeName)
outfile.close()
