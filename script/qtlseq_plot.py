#!/usr/bin/env python3

import os
import shutil
import sys
import argparse
from utils import time_stamp, clean_cmd, call_log
from refindex import RefIndex
from alignment import Alignment
from mpileup import Mpileup
from snpindex import Snpindex
from filt import Filt
from makeref import Makeref
from sim_range import Sim_range
from qtlplot import Qtlplot


class QTLseq(object):

    def __init__(self):

        parser = argparse.ArgumentParser(description='qtlseq2 pipeline 2022/11/15')

        parser.add_argument('-p2','--p2_name',
                            required=False,
                            default=False,
                            help='P2 name. [P2]')

        parser.add_argument('-f1','--f1',
                            required=False,
                            default=False,
                            help='F1 in vcf? [False]')

        parser.add_argument('-a','--bulkA_name',
                            required=False,
                            type=str,
                            default="A_bulk",
                            help='A bulk name. [A_bulk]')

        parser.add_argument('-b','--bulkB_name',
                            required=False,
                            type=str,
                            default="",
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
                            help='Sliding window size. [2000000000]')

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

        parser.add_argument('-txt','--txt',
                            required=True,
                            type=str,
                            default="",
                            help='Path of txt file')

        parser.add_argument('-filt_bias','--filt_bias',
                            required=False,
                            type=int,
                            default=0.3,
                            help='Filt SNP-index bias. [0.3]')

        args = parser.parse_args()

        self.p2_name=args.p2_name
        self.f1=args.f1
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
        self.plot_snpindex_txt=args.txt
        self.filt_bias=args.filt_bias

        self.sim_ind='{0}_{1}_sim.txt'.format(self.generation,self.individual)
        self.filt_count=1
        self.bulk_info=self.bulkA_name

        if self.p2_name :
            self.filt_count=self.filt_count+1
            
        if self.f1 :
            self.filt_count=self.filt_count+1

        if self.bulkB_name != "":
            self.bulk_info=self.bulk_info+"\t"+self.bulkB_name
            
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
        self.step30()
        self.step40()

if __name__ == '__main__':
    print(time_stamp(), 'start to run QTL-seq.', flush=True)
    QTLseq().run()
    print(time_stamp(), 'QTL-seq successfully finished.\n', flush=True)

