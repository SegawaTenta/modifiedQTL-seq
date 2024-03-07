import os
import random
import shutil
from mqtlseq.utils import time_stamp, clean_cmd, call_log

class Sim_range(object):

    def __init__(self,output_dir_aln,individural,sim_replication,generation,mindepth,maxdepth):

        self.output_dir = output_dir_aln
        self.individural = individural
        self.sim_replication = sim_replication
        self.generation = generation
        self.mindepth = mindepth
        self.maxdepth = maxdepth

    def run(self):

        print(time_stamp(),
              'start to simuration SNP-index value of confidence interval.',
              flush=True)

        p99_d=int(int(self.sim_replication)*0.005)
        p99_u=int(int(self.sim_replication)*0.995)
        p95_d=int(int(self.sim_replication)*0.025)
        p95_u=int(int(self.sim_replication)*0.975)

        if self.generation=="F2":

            if os.path.isfile("{0}/mqtlseq/{1}_{2}_sim.txt".format(os.path.dirname(__file__),self.output_dir,self.generation,self.individural)):
                shutil.copy("{0}/mqtlseq/{1}_{2}_sim.txt".format(os.path.dirname(__file__),self.output_dir,self.generation,self.individural), "{0}/{1}_{2}_sim.txt".format(self.output_dir,self.generation,self.individural))
            else:
                output_file = open("{0}/{1}_{2}_sim.txt".format(self.output_dir,self.generation,self.individural), "w")

                for i in range(int(self.mindepth),int(self.maxdepth)+1):

                    index_list=[]
                    dindex_list=[]

                    for rep in range(int(self.sim_replication)):

                        allele_pattern=[0,1]
                        all_allele1=[]
                        all_allele2=[]

                        for ind in range(int(self.individural)):

                            allele=random.choice(allele_pattern)
                            all_allele1.append(allele)
                            allele=random.choice(allele_pattern)
                            all_allele1.append(allele)

                            allele=random.choice(allele_pattern)
                            all_allele2.append(allele)
                            allele=random.choice(allele_pattern)
                            all_allele2.append(allele)

                        snp_count1=0
                        snp_count2=0

                        for ii in range(i):

                            allele1=random.choice(all_allele1)
                            allele2=random.choice(all_allele2)
                            snp_count1=snp_count1+allele1
                            snp_count2=snp_count2+allele2

                        snpindex1=snp_count1/(ii+1)
                        snpindex2=snp_count2/(ii+1)
                        deltasnpindex=snpindex1-snpindex2
                        index_list.append(snpindex1)
                        dindex_list.append(deltasnpindex)

                    index_list.sort()
                    dindex_list.sort()

                    p99_d_index=index_list[p99_d]
                    p99_d_index=p99_d_index*1000000
                    p99_d_index=int(p99_d_index)
                    p99_d_index=p99_d_index/1000000

                    p99_u_index=index_list[p99_u]
                    p99_u_index=p99_u_index*1000000
                    p99_u_index=int(p99_u_index)
                    p99_u_index=p99_u_index/1000000

                    p95_d_index=index_list[p95_d]
                    p95_d_index=p95_d_index*1000000
                    p95_d_index=int(p95_d_index)
                    p95_d_index=p95_d_index/1000000

                    p95_u_index=index_list[p95_u]
                    p95_u_index=p95_u_index*1000000
                    p95_u_index=int(p95_u_index)
                    p95_u_index=p95_u_index/1000000


                    dp99_d_index=dindex_list[p99_d]
                    dp99_d_index=dp99_d_index*1000000
                    dp99_d_index=int(dp99_d_index)
                    dp99_d_index=dp99_d_index/1000000

                    dp99_u_index=dindex_list[p99_u]
                    dp99_u_index=dp99_u_index*1000000
                    dp99_u_index=int(dp99_u_index)
                    dp99_u_index=dp99_u_index/1000000

                    dp95_d_index=dindex_list[p95_d]
                    dp95_d_index=dp95_d_index*1000000
                    dp95_d_index=int(dp95_d_index)
                    dp95_d_index=dp95_d_index/1000000

                    dp95_u_index=dindex_list[p95_u]
                    dp95_u_index=dp95_u_index*1000000
                    dp95_u_index=int(dp95_u_index)
                    dp95_u_index=dp95_u_index/1000000

                    output='{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(i,p99_d_index,p99_u_index,p95_d_index,p95_u_index,dp99_d_index,dp99_u_index,dp95_d_index,dp95_u_index)
                    output_file.write(output)

        if self.generation=="RIL":

            if os.path.isfile("{0}/mqtlseq/{1}_{2}_sim.txt".format(os.path.dirname(__file__),self.output_dir,self.generation,self.individural)):
                shutil.copy("{0}/mqtlseq/{1}_{2}_sim.txt".format(os.path.dirname(__file__),self.output_dir,self.generation,self.individural), "{0}/{1}_{2}_sim.txt".format(self.output_dir,self.generation,self.individural))
            else:

                output_file = open("{0}/{1}_{2}_sim.txt".format(self.output_dir,self.generation,self.individural), "w")

                for i in range(int(self.mindepth),int(self.maxdepth)+1):

                    index_list=[]
                    dindex_list=[]

                    for rep in range(int(self.sim_replication)):

                        allele_pattern=[0,1]
                        all_allele1=[]
                        all_allele2=[]

                        for ind in range(int(self.individural)):

                            allele=random.choice(allele_pattern)
                            all_allele1.append(allele)
                            all_allele1.append(allele)

                            allele=random.choice(allele_pattern)
                            all_allele2.append(allele)
                            all_allele2.append(allele)

                        snp_count1=0
                        snp_count2=0

                        for ii in range(i):

                            allele1=random.choice(all_allele1)
                            allele2=random.choice(all_allele2)
                            snp_count1=snp_count1+allele1
                            snp_count2=snp_count2+allele2

                        snpindex1=snp_count1/(ii+1)
                        snpindex2=snp_count2/(ii+1)
                        deltasnpindex=snpindex1-snpindex2
                        index_list.append(snpindex1)
                        dindex_list.append(deltasnpindex)

                    index_list.sort()
                    dindex_list.sort()

                    p99_d_index=index_list[p99_d]
                    p99_d_index=p99_d_index*1000000
                    p99_d_index=int(p99_d_index)
                    p99_d_index=p99_d_index/1000000

                    p99_u_index=index_list[p99_u]
                    p99_u_index=p99_u_index*1000000
                    p99_u_index=int(p99_u_index)
                    p99_u_index=p99_u_index/1000000

                    p95_d_index=index_list[p95_d]
                    p95_d_index=p95_d_index*1000000
                    p95_d_index=int(p95_d_index)
                    p95_d_index=p95_d_index/1000000

                    p95_u_index=index_list[p95_u]
                    p95_u_index=p95_u_index*1000000
                    p95_u_index=int(p95_u_index)
                    p95_u_index=p95_u_index/1000000

                    dp99_d_index=dindex_list[p99_d]
                    dp99_d_index=dp99_d_index*1000000
                    dp99_d_index=int(dp99_d_index)
                    dp99_d_index=dp99_d_index/1000000

                    dp99_u_index=dindex_list[p99_u]
                    dp99_u_index=dp99_u_index*1000000
                    dp99_u_index=int(dp99_u_index)
                    dp99_u_index=dp99_u_index/1000000

                    dp95_d_index=dindex_list[p95_d]
                    dp95_d_index=dp95_d_index*1000000
                    dp95_d_index=int(dp95_d_index)
                    dp95_d_index=dp95_d_index/1000000

                    dp95_u_index=dindex_list[p95_u]
                    dp95_u_index=dp95_u_index*1000000
                    dp95_u_index=int(dp95_u_index)
                    dp95_u_index=dp95_u_index/1000000

                    output='{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(i,p99_d_index,p99_u_index,p95_d_index,p95_u_index,dp99_d_index,dp99_u_index,dp95_d_index,dp95_u_index)
                    output_file.write(output)

        if self.generation=="BC1F1":

            if os.path.isfile("{0}/mqtlseq/{1}_{2}_sim.txt".format(os.path.dirname(__file__),self.output_dir,self.generation,self.individural)):
                shutil.copy("{0}/mqtlseq/{1}_{2}_sim.txt".format(os.path.dirname(__file__),self.output_dir,self.generation,self.individural), "{0}/{1}_{2}_sim.txt".format(self.output_dir,self.generation,self.individural))
            else:

                output_file = open("{0}/{1}_{2}_sim.txt".format(self.output_dir,self.generation,self.individural), "w")

                for i in range(int(self.mindepth),int(self.maxdepth)+1):

                    index_list=[]
                    dindex_list=[]

                    for rep in range(int(self.sim_replication)):

                        allele_pattern=[0,1]
                        all_allele1=[]
                        all_allele2=[]

                        for ind in range(int(self.individural)):

                            allele=random.choice(allele_pattern)
                            all_allele1.append(0)
                            all_allele1.append(allele)

                            allele=random.choice(allele_pattern)
                            all_allele2.append(0)
                            all_allele2.append(allele)

                        snp_count1=0
                        snp_count2=0

                        for ii in range(i):

                            allele1=random.choice(all_allele1)
                            allele2=random.choice(all_allele2)
                            snp_count1=snp_count1+allele1
                            snp_count2=snp_count2+allele2

                        snpindex1=snp_count1/(ii+1)
                        snpindex2=snp_count2/(ii+1)
                        deltasnpindex=snpindex1-snpindex2
                        index_list.append(snpindex1)
                        dindex_list.append(deltasnpindex)

                    index_list.sort()
                    dindex_list.sort()

                    p99_d_index=index_list[p99_d]
                    p99_d_index=p99_d_index*1000000
                    p99_d_index=int(p99_d_index)
                    p99_d_index=p99_d_index/1000000

                    p99_u_index=index_list[p99_u]
                    p99_u_index=p99_u_index*1000000
                    p99_u_index=int(p99_u_index)
                    p99_u_index=p99_u_index/1000000

                    p95_d_index=index_list[p95_d]
                    p95_d_index=p95_d_index*1000000
                    p95_d_index=int(p95_d_index)
                    p95_d_index=p95_d_index/1000000

                    p95_u_index=index_list[p95_u]
                    p95_u_index=p95_u_index*1000000
                    p95_u_index=int(p95_u_index)
                    p95_u_index=p95_u_index/1000000

                    dp99_d_index=dindex_list[p99_d]
                    dp99_d_index=dp99_d_index*1000000
                    dp99_d_index=int(dp99_d_index)
                    dp99_d_index=dp99_d_index/1000000

                    dp99_u_index=dindex_list[p99_u]
                    dp99_u_index=dp99_u_index*1000000
                    dp99_u_index=int(dp99_u_index)
                    dp99_u_index=dp99_u_index/1000000

                    dp95_d_index=dindex_list[p95_d]
                    dp95_d_index=dp95_d_index*1000000
                    dp95_d_index=int(dp95_d_index)
                    dp95_d_index=dp95_d_index/1000000

                    dp95_u_index=dindex_list[p95_u]
                    dp95_u_index=dp95_u_index*1000000
                    dp95_u_index=int(dp95_u_index)
                    dp95_u_index=dp95_u_index/1000000

                    output='{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(i,p99_d_index,p99_u_index,p95_d_index,p95_u_index,dp99_d_index,dp99_u_index,dp95_d_index,dp95_u_index)
                    output_file.write(output)

        print(time_stamp(),
              'simuration SNP-index value successfully finished.',
              flush=True)
