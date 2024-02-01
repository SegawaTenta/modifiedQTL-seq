import sys
import os
import re
import subprocess as sbp
import random
import pandas as pd
import scipy.stats as st
from utils import time_stamp, clean_cmd, call_log

class Filt(object):

    def __init__(self,mindepth,maxdepth,replicatipn,output_dir_aln,target_vcf,filt_pattern,filt_line,p1_name,p2_name):
        self.mindepth = mindepth
        self.maxdepth = maxdepth
        self.output_dir = output_dir_aln
        self.target_vcf = target_vcf
        self.filt_pattern = filt_pattern
        self.filt_line = filt_line
        self.p1_name = p1_name
        self.p2_name = p2_name
        self.replicatipn=replicatipn
        self.mydict = {}
        self.output_file2 = open("{0}/F1_sim.txt".format(self.output_dir), "w")

    def readsim_f1(self):

        if os.path.isfile("script/sim/F1_sim.txt"):
            sim_file = open("script/sim/F1_sim.txt", "r")
            for sim_file_line in sim_file:
                sim_file_line = sim_file_line.replace('\n','')
                array = sim_file_line.split('\t')
                depth=int(array[0])
                under=str(array[1])
                top=str(array[2])
                self.output_file2.write("{0}\t{1}\t{2}\n".format(depth,under,top))
                self.mydict[depth] = under+":"+top

    def sim_f1(self,depth):
        p95_d=int(int(self.replicatipn)*0.025)
        p95_u=int(int(self.replicatipn)*0.975)
        index_list=[]

        for rep in range(self.replicatipn):
            all_allele=[0,1]
            snp_count=0

            for ii in range(depth):
                allele=random.choice(all_allele)
                snp_count=snp_count+allele

            snp_index=snp_count/(ii+1)
            index_list.append(snp_index)

        index_list.sort()

        p95_d_index=index_list[p95_d]
        p95_d_index=p95_d_index*1000000
        p95_d_index=int(p95_d_index)
        p95_d_index=p95_d_index/1000000

        p95_u_index=index_list[p95_u]
        p95_u_index=p95_u_index*1000000
        p95_u_index=int(p95_u_index)
        p95_u_index=p95_u_index/1000000

        output2='{0}\t{1}\t{2}\n'.format(i,p95_d_index,p95_u_index)
        self.output_file2.write(output2)

        indexrange=str(p95_d_index)+":"+str(p95_u_index)
        self.mydict[i] = indexrange

    def run(self):

        print(time_stamp(),
              'start to filt SNP-index table.',
              flush=True)

        run_pattern = self.filt_pattern.split(':')
        run_line = self.filt_line.split(':')
        run_vcf = self.target_vcf
        run_num=len(run_line)

        for num in range(run_num):

            run_vcf_name=re.sub('.+/', '',run_vcf)

            filt_name=""
            filt_target=""

            if run_pattern[num] == "1" or run_pattern[num] == "2":
                filt_name="homo"
            elif run_pattern[num] == "3":
                filt_name="hetero"
                self.readsim_f1()
                # self.sim_f1()

            if run_line[num] == "1":
                filt_target=self.p1_name
            elif run_line[num] == "2":
                filt_target=self.p2_name
            elif run_line[num] == "3":
                filt_target="F1"

            name=filt_target+"_"+filt_name
            run_vcf_name=re.sub(filt_target,name,run_vcf_name)

            vcf = open(run_vcf, "r")
            output_file = open("{0}/{1}".format(self.output_dir,run_vcf_name), "w")
            run_vcf="{0}/{1}".format(self.output_dir,run_vcf_name)
            vcf_line = vcf.readline()

            while vcf_line:

                vcf_line=vcf_line.replace('\n', '')
                colom = vcf_line.split()

                depth_colom_num = 2 + 2 * (num+1)
                snpindex_colom_num = 3 + 2 * (num+1)

                depth = colom[depth_colom_num]
                snpindex = colom[snpindex_colom_num]

                if int(depth)>=self.mindepth and int(depth)<=self.maxdepth:

                    if run_pattern[num]=="1":
                        if float(snpindex) == 1:
                            output_file.write("{}\n".format(vcf_line))
                    elif run_pattern[num]=="2":
                        if float(snpindex) == 0:
                            output_file.write("{}\n".format(vcf_line))
                    elif run_pattern[num]=="3":

                        if int(depth) in self.mydict:
                            index_range = self.mydict[int(depth)].split(":")
                            under_index=index_range[0]
                            top_index=index_range[1]

                            if float(snpindex)>=float(under_index) and float(snpindex)<=float(top_index):

                                x_mutdepth=round(float(snpindex)*int(depth))
                                x_wtdepth=int(depth)-x_mutdepth
                                df = pd.DataFrame([[0,x_mutdepth], [colom[4],x_wtdepth]])
                                odds,p = st.fisher_exact(df)
                                if p < 0.05:
                                    output_file.write("{}\n".format(vcf_line))
                        else:
                            self.sim_f1(int(depth))
                            index_range = self.mydict[int(depth)].split(":")
                            under_index=index_range[0]
                            top_index=index_range[1]

                            if float(snpindex)>=float(under_index) and float(snpindex)<=float(top_index):

                                x_mutdepth=round(float(snpindex)*int(depth))
                                x_wtdepth=int(depth)-x_mutdepth
                                df = pd.DataFrame([[0,x_mutdepth], [colom[4],x_wtdepth]])
                                odds,p = st.fisher_exact(df)
                                if p < 0.05:
                                    output_file.write("{}\n".format(vcf_line))

                vcf_line = vcf.readline()

            vcf.close()
            output_file.close()

        print(time_stamp(),
              'filt SNP-index table successfully finished.',
              flush=True)
