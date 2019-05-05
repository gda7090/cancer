#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import gzip
import re
import time
import argparse
import getpass
parser = argparse.ArgumentParser(description="sum diff result to plot lanscape")
parser.add_argument('--compare',help="compare name ,for example:T15_vs_N15,T14_vs_N14,T18_vs_N18",required=True)
parser.add_argument('--diff_dir',help="diff dir",required=True)
parser.add_argument('--n',help="diff number",required=True)
parser.add_argument('--type',help="project type",required=True)
argv = parser.parse_args()
compare=argv.compare
n=int(argv.n)
type=(argv.type)
out=open("lanscape_data_result.xls",'w')
out.write("ID\t"+"\t".join(compare.split(","))+"\n")
###--------prepare_data_file----------###
compare_result_dic_list=[]
diff_id=[]
for each_compare in compare.split(","):
    diff_file=diff_dir+"/"+each_compare+"/"+each_compare+"."+type+".diffgene.xls"
    diffdic={}
    for eachline in open(diff_file).readlines()[1:]:
        each=eachline.strip().split("\t")
        diff_id.append(each[0])
        if float(each[5])>=0:
            diffdic[each[0]]='up'
        else:
            diffdic[each[0]]='down'
    compare_result_dic_list.append(diffdic)
for each_id in set(diff_id):
    out.write(each_id)
    for i in range(len(compare.split(","))):
        if each_id in compare_result_dic_list[i]:
            out.write("\t"+str(compare_result_dic_list[i][each_id]))
        else:
            out.write("\t")
    out.write("\n")

number_to_show=n
lanscape_sum_result=open("lanscape_sum_result.xls",'w')
lanscape_sum_result.write("ID\t"+"\t".join(compare.split(","))+"\n")
for each_id in set(diff_id):
    a=0
    for i in range(len(compare.split(","))):
        if each_id in compare_result_dic_list[i]:
            a+=1
        else:
            pass
    if a>=number_to_show:
        lanscape_sum_result.write(each_id)
        for i in range(len(compare.split(","))):
            if each_id in compare_result_dic_list[i]:
                lanscape_sum_result.write("\t"+str(compare_result_dic_list[i][each_id]))
            else:
                lanscape_sum_result.write("\t")
        lanscape_sum_result.write("\n")
os.system("cp Landscape_RNA/lanscape.R ./")
os.system("R-3.3.1/bin/R -f lanscape.R")
