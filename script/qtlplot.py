import sys
import os
import re
from utils import time_stamp, clean_cmd, call_log
import pandas as pd
import matplotlib.pyplot as plt
import math

class Qtlplot(object):

    def __init__(self,output_dir_aln,sim_file,target_vcf,mask_line,window_size,step_size,min_plot,bulk_info,mindepth,maxdepth,filt_bias):

        self.output_dir = output_dir_aln
        self.sim_file = sim_file
        self.target_vcf = target_vcf
        self.mask_line = mask_line
        self.window_size = window_size
        self.step_size = step_size
        self.min_plot = min_plot
        self.bulk_info = bulk_info
        self.mindepth = mindepth
        self.maxdepth = maxdepth
        self.filt_bias = filt_bias

    def read_sim(self):

        simfile = open(self.sim_file, "r")
        simfile_line = simfile.readline()
        self.mydict = {}

        while simfile_line:
            simfile_line=simfile_line.replace('\n', '')
            item = simfile_line.split()

            range_list=item[1]+":"+item[2]+":"+item[3]+":"+item[4]+":"+item[5]+":"+item[6]+":"+item[7]+":"+item[8]

            self.mydict[int(item[0])] = range_list
            simfile_line = simfile.readline()

    def add_sim1(self,name1,line):

        vcf = open(self.target_vcf, "r")
        vcf_line = vcf.readline()
        output_file1 = open("{0}/{1}.snp_index.sim.txt".format(self.output_dir,name1), "w")

        while vcf_line:

            vcf_line=vcf_line.replace('\n', '')
            item = vcf_line.split()

            target_depth1=item[int(line)*2+2]
            target_index1=item[int(line)*2+3]

            if int(target_depth1)>=int(self.mindepth) and int(target_depth1)<=int(self.maxdepth):
                sim_a=self.mydict[int(target_depth1)]
                sim_a=sim_a.replace(':', '\t')

                output1='{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(item[0],item[1],item[2],item[3],
                                                                     target_depth1,
                                                                     target_index1,
                                                                     sim_a)

                output_file1.write(output1)

            vcf_line = vcf.readline()

        vcf.close()
        output_file1.close()

    def add_sim2(self,name1,name2,line):

        vcf = open(self.target_vcf, "r")
        vcf_line = vcf.readline()
        output_file1 = open("{0}/{1}.snp_index.sim.txt".format(self.output_dir,name1), "w")
        output_file2 = open("{0}/{1}.snp_index.sim.txt".format(self.output_dir,name2), "w")
        output_file3 = open("{0}/{1}-{2}.snp_index.sim.txt".format(self.output_dir,name1,name2), "w")

        while vcf_line:

            vcf_line=vcf_line.replace('\n', '')
            item = vcf_line.split()
            target_depth1=item[int(line)*2+2]
            target_index1=item[int(line)*2+3]
            target_depth2=item[int(line)*2+4]
            target_index2=item[int(line)*2+5]

            target_depth3=0
            if int(target_depth1)>int(target_depth2):
                target_depth3=int(target_depth1)
            else:
                target_depth3=int(target_depth2)

            target_index3=float(target_index1)-float(target_index2)

            

            if int(target_depth1)>=int(self.mindepth) and int(target_depth2)>=int(self.mindepth) and int(target_depth1)<=int(self.maxdepth) and int(target_depth2)<=int(self.maxdepth) :
                if float(target_index1)<self.filt_bias and float(target_index2)<self.filt_bias:
                    a=1
                elif float(target_index1) > 1-self.filt_bias and float(target_index2) > 1-self.filt_bias:
                    a=1
                else:
                    sim_a=self.mydict[int(target_depth1)]
                    sim_a=sim_a.replace(':', '\t')
                    sim_b=self.mydict[int(target_depth2)]
                    sim_b=sim_b.replace(':', '\t')
                    sim_c=self.mydict[int(target_depth3)]
                    sim_c=sim_c.replace(':', '\t')

                    output1='{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(item[0],item[1],item[2],item[3],
                                                                        target_depth1,
                                                                        target_index1,
                                                                        sim_a)

                    output2='{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(item[0],item[1],item[2],item[3],
                                                                        target_depth2,
                                                                        target_index2,
                                                                        sim_b)

                    output3='{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(item[0],item[1],item[2],item[3],
                                                                        target_depth3,
                                                                        target_index3,
                                                                        sim_c)

                    output_file1.write(output1)
                    output_file2.write(output2)
                    output_file3.write(output3)

            vcf_line = vcf.readline()

        vcf.close()
        output_file1.close()
        output_file2.close()
        output_file3.close()

    def make_window(self,target_txt,name):

        window_mb=int(self.window_size)
        stepp_kb=int(self.step_size)

        vcf = open(target_txt, "r")
        output_file = open("{0}/{1}.{2}bp_window.{3}bp_step.txt".format(self.output_dir,name,int(window_mb),int(stepp_kb)), "w")

        df=pd.read_table(vcf,header=None)

        chr_item = df[0].unique()

        for chr in chr_item:

            df1 = df[df[0]==chr]
            max_posi=max(df1[1])

            step_posi=range(0, max_posi+int(self.step_size), int(self.step_size))

            for posi in step_posi:

                posi_t=posi+int(self.window_size)/2
                posi_u=posi-int(self.window_size)/2
                df2 = df1[df1[1]<=posi_t]
                df3 = df2[df2[1]>=posi_u]

                df3_line=len(df3)

                if df3_line>=self.min_plot:
                    ave1 = sum(df3[5]) / df3_line
                    ave1=ave1*1000000
                    ave1=int(ave1)
                    ave1=ave1/1000000

                    ave2 = sum(df3[6]) / df3_line
                    ave2=ave2*1000000
                    ave2=int(ave2)
                    ave2=ave2/1000000

                    ave3 = sum(df3[7]) / df3_line
                    ave3=ave3*1000000
                    ave3=int(ave3)
                    ave3=ave3/1000000

                    ave4 = sum(df3[8]) / df3_line
                    ave4=ave4*1000000
                    ave4=int(ave4)
                    ave4=ave4/1000000

                    ave5 = sum(df3[9]) / df3_line
                    ave5=ave5*1000000
                    ave5=int(ave5)
                    ave5=ave5/1000000

                    ave6 = sum(df3[10]) / df3_line
                    ave6=ave6*1000000
                    ave6=int(ave6)
                    ave6=ave6/1000000

                    ave7 = sum(df3[11]) / df3_line
                    ave7=ave7*1000000
                    ave7=int(ave7)
                    ave7=ave7/1000000

                    ave8 = sum(df3[12]) / df3_line
                    ave8=ave8*1000000
                    ave8=int(ave8)
                    ave8=ave8/1000000

                    ave9 = sum(df3[13]) / df3_line
                    ave9=ave9*1000000
                    ave9=int(ave9)
                    ave9=ave9/1000000

                    output="{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\n".format(chr,posi,ave1,ave2,ave3,ave4,ave5,ave6,ave7,ave8,ave9)
                    output_file.write(output)

        vcf.close()
        output_file.close()

    def singleplot(self,name,color):

        window_mb=int(self.window_size)
        stepp_kb=int(self.step_size)

        plotcolor=""
        if color==1:
            plotcolor="green"
        else:
            plotcolor="gold"


        window_txt = open("{0}/{1}.{2}bp_window.{3}bp_step.txt".format(self.output_dir,name,int(window_mb),int(stepp_kb)), "r")
        snp_txt = open("{0}/{1}.snp_index.sim.txt".format(self.output_dir,name,int(window_mb),int(stepp_kb)), "r")

        df1=pd.read_table(window_txt,header=None)
        df2=pd.read_table(snp_txt,header=None)

        maxposi1=max(df1[1])/1000000
        maxposi2=max(df2[1])/1000000
        maxposi=0

        if maxposi1>maxposi2:
            maxposi=maxposi1
        else:
            maxposi=maxposi2

        fig = plt.figure(figsize=(12,16),dpi=300)
        plt.subplots_adjust(wspace=0.4, hspace=0.6)

        chr_item = df1[0].unique()
        chr_num=len(chr_item)

        fig_col=chr_num/5
        fig_col=math.ceil(fig_col)
        fig_raw=5

        count1=0

        for chr in chr_item:

            count1=count1+1

            ax = fig.add_subplot(fig_raw, fig_col, count1)

            plt.title(chr,fontsize=9)
            plt.xlim(0,maxposi)
            plt.ylim(0,1)
            plt.xlabel("chr position(Mb)",fontsize=9)
            plt.ylabel("SNP-index",fontsize=9)
            plt.gca().spines['right'].set_visible(False)
            plt.gca().spines['top'].set_visible(False)

            df11=df1[df1[0]==chr]
            df22=df2[df2[0]==chr]

            ax.plot(df22[1]/1000000, df22[5], color=plotcolor,marker=".",markersize=4,linestyle='None')
            ax.plot(df11[1]/1000000, df11[3], color="orange")
            ax.plot(df11[1]/1000000, df11[4], color="orange")
            ax.plot(df11[1]/1000000, df11[5], color="lime")
            ax.plot(df11[1]/1000000, df11[6], color="lime")
            ax.plot(df11[1]/1000000, df11[2], color="red")
            ax.hlines([0.5], 0, maxposi, linestyles='dashed')

        fig.savefig("{0}/{1}.png".format(self.output_dir,name))

    def deltaplot(self,name):

        window_mb=int(self.window_size)
        stepp_kb=int(self.step_size)

        plotcolor="navy"

        window_txt = open("{0}/{1}.{2}bp_window.{3}bp_step.txt".format(self.output_dir,name,int(window_mb),int(stepp_kb)), "r")
        snp_txt = open("{0}/{1}.snp_index.sim.txt".format(self.output_dir,name,int(window_mb),int(stepp_kb)), "r")

        df1=pd.read_table(window_txt,header=None)
        df2=pd.read_table(snp_txt,header=None)

        maxposi1=max(df1[1])/1000000
        maxposi2=max(df2[1])/1000000
        maxposi=0

        if maxposi1>maxposi2:
            maxposi=maxposi1
        else:
            maxposi=maxposi2

        fig = plt.figure(figsize=(12,16),dpi=300)
        plt.subplots_adjust(wspace=0.4, hspace=0.6)

        chr_item = df1[0].unique()
        chr_num=len(chr_item)

        fig_col=chr_num/5
        fig_col=math.ceil(fig_col)
        fig_raw=5

        count1=0

        for chr in chr_item:

            count1=count1+1

            ax = fig.add_subplot(fig_raw, fig_col, count1)

            plt.title(chr,fontsize=9)
            plt.xlim(0,maxposi)
            plt.ylim(-1,1)
            plt.xlabel("chr position(Mb)",fontsize=9)
            plt.ylabel("SNP-index",fontsize=9)
            plt.gca().spines['right'].set_visible(False)
            plt.gca().spines['top'].set_visible(False)

            df11=df1[df1[0]==chr]
            df22=df2[df2[0]==chr]

            ax.plot(df22[1]/1000000, df22[5], color=plotcolor,marker=".",markersize=4,linestyle='None')
            ax.plot(df11[1]/1000000, df11[7], color="orange")
            ax.plot(df11[1]/1000000, df11[8], color="orange")
            ax.plot(df11[1]/1000000, df11[9], color="lime")
            ax.plot(df11[1]/1000000, df11[10], color="lime")
            ax.plot(df11[1]/1000000, df11[2], color="red")
            ax.hlines([0.0], 0, maxposi, linestyles='dashed')

        fig.savefig("{0}/{1}.png".format(self.output_dir,name))

    def run(self):

        print(time_stamp(),
              'start to plot QTL-seq graph.',
              flush=True)

        bulk_name = self.bulk_info.split('\t')
        bulk_count=len(bulk_name)

        self.read_sim()

        if bulk_count==1:

            bulk_line=int(self.mask_line)+1

            self.add_sim1(bulk_name[0],bulk_line)

            target_txt1="{0}/{1}.snp_index.sim.txt".format(self.output_dir,bulk_name[0])

            self.make_window(target_txt1,bulk_name[0])

            self.singleplot(bulk_name[0],1)


        elif bulk_count==2:

            bulk_line=int(self.mask_line)+1

            self.add_sim2(bulk_name[0],bulk_name[1],bulk_line)

            target_txt1="{0}/{1}.snp_index.sim.txt".format(self.output_dir,bulk_name[0])
            target_txt2="{0}/{1}.snp_index.sim.txt".format(self.output_dir,bulk_name[1])

            delta_name=bulk_name[0]+"-"+bulk_name[1]
            target_txt3="{0}/{1}.snp_index.sim.txt".format(self.output_dir,delta_name)

            self.make_window(target_txt1,bulk_name[0])
            self.make_window(target_txt2,bulk_name[1])
            self.make_window(target_txt3,delta_name)

            self.singleplot(bulk_name[0],1)
            self.singleplot(bulk_name[1],2)
            self.deltaplot(delta_name)


        print(time_stamp(),
              'ploting QTL-seq graph successfully finished.',
              flush=True)
