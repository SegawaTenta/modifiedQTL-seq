import sys
import os
import re
import subprocess as sbp
from utils import time_stamp, clean_cmd, call_log


class Snpindex(object):

    def __init__(self,output_dir_aln,target_vcf,out_name):
        self.output_dir = output_dir_aln
        self.target_vcf = target_vcf
        self.out_name = out_name

    def run(self):
        print(time_stamp(),
              'start to calculate SNP-index.',
              flush=True)

        input_vcf = open(self.target_vcf, "r")
        output_file = open("{0}/{1}.snp_index.txt".format(self.output_dir,self.out_name), "w")
        for vcf_line in input_vcf :
            # print(vcf_line)
            vcf_line=vcf_line.replace('\n', '')

            if vcf_line.startswith('#'):
                a="AAA"
            else:
                colom = vcf_line.split()
                colom_num = len(colom)
                chr = str(colom[0])
                posi = str(colom[1])
                ref_base = str(colom[3])
                mut_base = str(colom[4])

                snpindex_list=[chr,posi,ref_base,mut_base]
                which_print="yes"

                # print(snpindex_list)

                if len(ref_base)==1 and  len(mut_base)==1:

                    for num in range(9,colom_num):

                        info = colom[num]
                        colom2 = info.split(':')
                        colom3 = colom2[3].split(',')
                        ref_depth=colom3[0]
                        mut_depth=colom3[1]

                        depth=int(ref_depth)+int(mut_depth)
                        snpindex=0

                        # print(snpindex_list,depth)

                        if depth==0:
                            which_print="no"
                            # print(vcf_line)
                        else:
                            snpindex=int(mut_depth)/int(depth)
                            snpindex=snpindex*1000000
                            snpindex=int(snpindex)
                            snpindex=snpindex/1000000
                            snpindex_list.append(str(depth))
                            snpindex_list.append(str(snpindex))

                    if which_print=="yes":
                        snpindex_list="\t".join(snpindex_list)
                        output_file.write(snpindex_list+"\n")

        input_vcf.close()

        print(time_stamp(),
              'calculation of SNP-index successfully finished.',
              flush=True)
