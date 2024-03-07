#!/usr/bin/env python3

from mqtlseq.__init__ import __version__
import os
import shutil
import sys
import argparse
from mqtlseq.utils import time_stamp, clean_cmd, call_log
from mqtlseq.refindex import RefIndex
from mqtlseq.alignment import Alignment
from mqtlseq.mpileup import Mpileup
from mqtlseq.snpindex import Snpindex
from mqtlseq.filt import Filt
from mqtlseq.makeref import Makeref
from mqtlseq.sim_range import Sim_range
from mqtlseq.qtlplot import Qtlplot


class Mqtlseq_mpileup(object):

    def __init__(self):

        parser = argparse.ArgumentParser(description='qtlseq2 pipeline 2022/11/15')

        parser.add_argument('-r','--ref_fa',
                            required=True,
                            type=str,
                            help='Referance fasta file.')

        parser.add_argument('-p1','--p1_name',
                            required=False,
                            type=str,
                            default="P1",
                            help='P1 name. [P1]')

        parser.add_argument('-p2','--p2_name',
                            required=False,
                            type=str,
                            default="P2",
                            help='P2 name. [P2]')

        parser.add_argument('-a','--bulkA_name',
                            required=False,
                            type=str,
                            default="A_bulk",
                            help='A bulk name. [A_bulk]')

        parser.add_argument('-b','--bulkB_name',
                            required=False,
                            type=str,
                            default="B_bulk",
                            help='B bulk name. [B_bulk]')

        parser.add_argument('-o','--output_dir',
                            required=False,
                            type=str,
                            default="QTL-seq",
                            help='Output directory. [QTL-seq]')

        parser.add_argument('-t','--thread',
                            required=False,
                            type=int,
                            default=1,
                            help='Thread. [1]')

        parser.add_argument('-mindepth','--mindepth',
                            required=False,
                            type=int,
                            default=10,
                            help='Minimum depth of SNP position used calculating SNP-index. [10]')

        parser.add_argument('-maxdepth','--maxdepth',
                            required=False,
                            type=int,
                            default=99,
                            help='Maximum depth of SNP position used calculating SNP-index. [99]')

        parser.add_argument('-sim_replication','--sim_replication',
                            required=False,
                            type=int,
                            default=10000,
                            help='Number of replication for simulating confidence interval of SNP-index. [10000]')

        parser.add_argument('-individual','--individual',
                            required=False,
                            type=int,
                            default=20,
                            help='Individual number of a bulk. [20]')

        parser.add_argument('-generation','--generation',
                            required=False,
                            type=str,
                            default="F2",
                            choices=['F2', 'RIL', 'BC1F1'],
                            help='Generation used for QTL-seq. You can chouse "F2" or "RIL" or "BC1F1". [F2]')

        parser.add_argument('-window_size','--window_size',
                            required=False,
                            type=int,
                            default=2000000,
                            help='Sliding window size. [2000000]')

        parser.add_argument('-step_size','--step_size',
                            required=False,
                            type=int,
                            default=50000,
                            help='Step size. [50000]')

        parser.add_argument('-min_plot','--min_plot',
                            required=False,
                            type=int,
                            default=10,
                            help='Minumum plot in window size. [10]')

        parser.add_argument('-bt','--bam_txt',
                            required=True,
                            type=str,
                            default="",
                            help='Path of reads_PATH.txt')

        parser.add_argument('-filt_bias','--filt_bias',
                            required=False,
                            type=int,
                            default=0.3,
                            help='Filt SNP-index bias. [0.3]')
        
        parser.add_argument('-bq','--base_quality',
                            required=False,
                            type=int,
                            default=13,
                            help='Minimum base quality. [13]')

        parser.add_argument('-mq','--mapping_quality',
                            required=False,
                            type=int,
                            default=40,
                            help='Minimum mapping quality. [40]')

        args = parser.parse_args()

        self.ref_fa=args.ref_fa
        self.p1_name=args.p1_name
        self.p2_name=args.p2_name
        self.bulkA_name=args.bulkA_name
        self.bulkB_name=args.bulkB_name
        self.output_dir=args.output_dir
        self.thread=args.thread
        self.mindepth=args.mindepth
        self.maxdepth=args.maxdepth
        self.sim_replication=args.sim_replication
        self.individual=args.individual
        self.generation=args.generation
        self.window_size=args.window_size
        self.step_size=args.step_size
        self.min_plot=args.min_plot
        self.bam_txt=args.bam_txt
        self.filt_bias=args.filt_bias
        self.base_quality=args.base_quality
        self.mapping_quality=args.mapping_quality

        self.filt_pattern="2"
        self.filt_line="1"

        self.sim_ind='{0}_{1}_sim.txt'.format(self.generation,self.individual)

        self.output_name=""
        self.bam_list=""

        self.p1_bam=""
        self.p2_bam=""
        self.f1_bam=""
        self.bulkA_bam=""
        self.bulkB_bam=""
        in_count=0
        self.filt_count=1
        self.bulk_info=self.bulkA_name
        self.plot_name=self.p1_name+"_homo"

        r_txt = open(self.bam_txt, "r")
        for r_txt_line in r_txt:
            r_txt_line=r_txt_line.replace('\n', '')
            if r_txt_line == "P1 bam":
                in_count=1
            elif r_txt_line == "P2 bam":
                in_count=2
            elif r_txt_line == "A bulk bam":
                in_count=3
            elif r_txt_line == "B bulk bam":
                in_count=4
            elif r_txt_line == "F1 bam":
                in_count=5
            elif in_count==1:
                self.p1_bam=r_txt_line
            elif in_count==2:
                self.p2_bam=r_txt_line
            elif in_count==3:
                self.bulkA_bam=r_txt_line
            elif in_count==4:
                self.bulkB_bam=r_txt_line
            elif in_count==5:
                self.f1_bam=r_txt_line
        r_txt.close()

        if self.p1_bam == "":
            print('Not setting P1 bam. Cheak {0}'.format(self.bam_txt))
            sys.exit()

        if self.bulkA_bam == "":
            print('Not setting A bulk bam. Cheak {0}'.format(self.bam_txt))
            sys.exit()

        self.output_name=self.p1_name
        self.bam_list=self.p1_bam

        if self.p2_bam != "":
            self.output_name=self.output_name+"."+self.p2_name
            self.bam_list=self.bam_list+" "+self.p2_bam
            self.filt_pattern=self.filt_pattern+":1"
            self.filt_line=self.filt_line+":2"
            self.filt_count=self.filt_count+1
            self.plot_name=self.plot_name+"."+self.p2_name+"_homo"

        if self.f1_bam != "":
            self.output_name=self.output_name+".F1"
            self.bam_list=self.bam_list+" "+self.f1_bam
            self.filt_pattern=self.filt_pattern+":3"
            self.filt_line=self.filt_line+":3"
            self.filt_count=self.filt_count+1
            self.plot_name=self.plot_name+".F1_hetero"

        self.output_name=self.output_name+"."+self.bulkA_name
        self.bam_list=self.bam_list+" "+self.bulkA_bam
        self.plot_name=self.plot_name+"."+self.bulkA_name

        if self.bulkB_bam != "":
            self.output_name=self.output_name+"."+self.bulkB_name
            self.bam_list=self.bam_list+" "+self.bulkB_bam
            self.bulk_info=self.bulk_info+"\t"+self.bulkB_name
            self.plot_name=self.plot_name+"."+self.bulkB_name

        self.plot_snpindex_txt="{0}/24_filt/{1}.snp_index.txt".format(self.output_dir,self.plot_name)

    def step22(self):
        os.makedirs('{}/22_vcf'.format(self.output_dir), exist_ok=True)
        output_dir_aln=self.output_dir+"/22_vcf"
        mpileup = Mpileup(output_dir_aln,self.bam_list,self.ref_fa,self.thread,self.output_name,self.base_quality,self.mapping_quality)
        mpileup.run()

    def step23(self):
        os.makedirs('{}/23_snpindex'.format(self.output_dir), exist_ok=True)
        output_dir_aln=self.output_dir+"/23_snpindex"
        target_vcf="{0}/22_vcf/{1}.vcf".format(self.output_dir,self.output_name)
        snpindex = Snpindex(output_dir_aln,target_vcf,self.output_name)
        snpindex.run()

    def step24(self):
        os.makedirs('{}/24_filt'.format(self.output_dir), exist_ok=True)
        output_dir_aln=self.output_dir+"/24_filt"
        target_snpindex_txt="{0}/23_snpindex/{1}.snp_index.txt".format(self.output_dir,self.output_name)
        filt = Filt(self.mindepth,self.maxdepth,self.sim_replication,output_dir_aln,target_snpindex_txt,self.filt_pattern,self.filt_line,self.p1_name,self.p2_name)
        filt.run()

    def step30(self):
        os.makedirs('{}/30_sim'.format(self.output_dir), exist_ok=True)

        if os.path.isfile("script/sim/{0}".format(self.sim_ind)) and self.maxdepth<100:
            if os.path.isfile("{0}/30_sim/{1}".format(self.output_dir,self.sim_ind)):
                print("{0}/30_sim/{1} already symlinked".format(self.output_dir,self.sim_ind))
            else:
                shutil.copyfile("script/sim/{0}".format(self.sim_ind),"{0}/30_sim/{1}".format(self.output_dir,self.sim_ind))
        else:
            output_dir_aln=self.output_dir+"/30_sim"
            sim_range = Sim_range(output_dir_aln,self.individual,self.sim_replication,self.generation,self.mindepth,self.maxdepth)
            sim_range.run()         

    def step40(self):
        os.makedirs('{}/40_plot'.format(self.output_dir), exist_ok=True)
        sim_file="{0}/30_sim/{1}".format(self.output_dir,self.sim_ind)
        output_dir_aln=self.output_dir+"/40_plot"
        qtlplot=Qtlplot(output_dir_aln,sim_file,self.plot_snpindex_txt,self.filt_count,self.window_size,self.step_size,self.min_plot,self.bulk_info,self.mindepth,self.maxdepth,self.filt_bias)
        qtlplot.run()

    def run(self):
        os.makedirs('{}'.format(self.output_dir), exist_ok=True)
        os.makedirs('{}/log'.format(self.output_dir), exist_ok=True)
        self.step22()
        self.step23()
        self.step24()
        self.step30()
        self.step40()

def main():
    Mqtlseq_mpileup().run()

if __name__ == '__main__':
    print(time_stamp(), 'start to run modified QTL-seq.', flush=True)
    Mqtlseq_mpileup().main()
    print(time_stamp(), 'modified QTL-seq successfully finished.\n', flush=True)

