import sys
import os
import subprocess as sbp
from utils import time_stamp, clean_cmd, call_log



class Alignment(object):

    def __init__(self,ref,alignment_table,output_dir_aln,thread,bwa_type):
        self.output_dir = output_dir_aln
        self.ref = ref
        self.alignment_table = alignment_table
        self.thread = thread
        self.bwa_type = bwa_type

    def run(self):
        print(time_stamp(),
              'start to align reads by BWA.',
              flush=True)

        # print(self.alignment_table)

        aln_item = self.alignment_table.split('\n')
        aln_item.pop()
        # print(aln_item)

        count = 0
        keep_name = ""
        rm_list = ""
        
        command_mem = open(self.output_dir+"/mem.txt", "w")
        command_view = open(self.output_dir+"/view.txt", "w")
        command_merge = open(self.output_dir+"/merge.txt", "w")
        command_sort = open(self.output_dir+"/sort.txt", "w")
        command_rmdup = open(self.output_dir+"/rmdup.txt", "w")
        command_index = open(self.output_dir+"/index.txt", "w")
        command_coverage = open(self.output_dir+"/coverage.txt", "w")
        command_rm = open(self.output_dir+"/rm.txt", "w")

        for array_line in aln_item :
            # print(array_line)
            array = array_line.split('\t')
            name = array[0]
            fq1 = array[1]
            fq2 = array[2]

            sym_fq1 = '{0}/{1}.{2}.1.fq'.format(self.output_dir,name,count)
            sym_fq2 = '{0}/{1}.{2}.2.fq'.format(self.output_dir,name,count)

            if os.path.isfile(sym_fq1):
                print("{0} already symlinked".format(sym_fq1))
            else:
                os.symlink(fq1, sym_fq1)
                print(fq1+" -> "+sym_fq1)

            if os.path.isfile(sym_fq2):
                print("{0} already symlinked".format(sym_fq2))
            else:
                os.symlink(fq2, sym_fq2)
                print(fq2+" -> "+sym_fq2)

            rm_list=rm_list+" "+sym_fq1
            rm_list=rm_list+" "+sym_fq2

            mem_format = '{4} mem -t 10 {0} {3}/{1}.{2}.1.fq {3}/{1}.{2}.2.fq 1>{3}/{1}.{2}.sam 2>>{3}/../log/bwa.log\n'.format(self.ref,
                                                                                                name,
                                                                                                count,
                                                                                                self.output_dir,
                                                                                                self.bwa_type)

            view_format = 'samtools view -bS -F 4 {0}/{1}.{2}.sam 1>{0}/{1}.{2}.bam 2>>{0}/../log/samtools.log\n'.format(self.output_dir,
                                                                                            name,
                                                                                            count)
            command_mem.write(mem_format)
            command_view.write(view_format)

            rm_list1='rm {0}/{1}.{2}.sam\n'.format(self.output_dir,name,count)
            rm_list2='rm {0}/{1}.{2}.bam\n'.format(self.output_dir,name,count)
            rm_list3='rm {0}\n'.format(sym_fq1)
            rm_list4='rm {0}\n'.format(sym_fq2)


            command_rm.write(rm_list1)
            command_rm.write(rm_list2)
            command_rm.write(rm_list3)
            command_rm.write(rm_list4)


            if keep_name != name:
                keep_name = name

                merge_format = 'samtools merge {0}/{1}.merged.bam {0}/{1}.*.bam 2>>{0}/../log/samtools.log\n'.format(self.output_dir,name,count)

                sort_format = 'samtools sort -T {1}_sorting -@ 4 {0}/{1}.merged.bam -o {0}/{1}.sort.bam 2>>{0}/../log/samtools.log\n'.format(self.output_dir,name)

                rmdup_format = 'samtools rmdup {0}/{1}.sort.bam {0}/{1}.bam 2>>{0}/../log/samtools.log\n'.format(self.output_dir,name)

                index_format = 'samtools index {0}/{1}.bam 2>>{0}/../log/samtools.log\n'.format(self.output_dir,name)

                coverage_format = 'samtools coverage {0}/{1}.bam -o {0}/{1}.coverage.txt  2>>{0}/../log/samtools.log\n'.format(self.output_dir,name)

                command_merge.write(merge_format)
                command_sort.write(sort_format)
                command_rmdup.write(rmdup_format)
                command_index.write(index_format)
                command_coverage.write(coverage_format)

                rm_list5='rm {0}/{1}.merged.bam\n'.format(self.output_dir,name)
                rm_list6='rm {0}/{1}.sort.bam\n'.format(self.output_dir,name)

                command_rm.write(rm_list5)
                command_rm.write(rm_list6)

            count=count+1

        command_mem.close()
        command_view.close()
        command_merge.close()
        command_sort.close()
        command_rmdup.close()
        command_index.close()
        command_coverage.close()
        command_rm.close()
        

        cmd_mem = 'cat {0}/mem.txt|xargs -P{1} -I % sh -c %'.format(self.output_dir,self.thread)
        cmd_mem = clean_cmd(cmd_mem)

        cmd_view = 'cat {0}/view.txt|xargs -P{1} -I % sh -c %'.format(self.output_dir,self.thread)
        cmd_view = clean_cmd(cmd_view)

        cmd_merge = 'cat {0}/merge.txt|xargs -P{1} -I % sh -c %'.format(self.output_dir,self.thread)
        cmd_merge = clean_cmd(cmd_merge)

        cmd_sort = 'cat {0}/sort.txt|xargs -P{1} -I % sh -c %'.format(self.output_dir,self.thread)
        cmd_sort = clean_cmd(cmd_sort)

        cmd_rmdup = 'cat {0}/rmdup.txt|xargs -P{1} -I % sh -c %'.format(self.output_dir,self.thread)
        cmd_rmdup = clean_cmd(cmd_rmdup)

        cmd_index = 'cat {0}/index.txt|xargs -P{1} -I % sh -c %'.format(self.output_dir,self.thread)
        cmd_index = clean_cmd(cmd_index)

        cmd_coverage = 'cat {0}/coverage.txt|xargs -P{1} -I % sh -c %'.format(self.output_dir,self.thread)
        cmd_coverage = clean_cmd(cmd_coverage)

        cmd_rm = 'cat {0}/rm.txt|xargs -P{1} -I % sh -c %'.format(self.output_dir,self.thread)
        cmd_rm = clean_cmd(cmd_rm)
        
        sbp.run(cmd_mem,
                stdout=sbp.DEVNULL,
                stderr=sbp.DEVNULL,
                shell=True,
                check=True)

        sbp.run(cmd_view,
                stdout=sbp.DEVNULL,
                stderr=sbp.DEVNULL,
                shell=True,
                check=True)

        sbp.run(cmd_merge,
                stdout=sbp.DEVNULL,
                stderr=sbp.DEVNULL,
                shell=True,
                check=True)

        sbp.run(cmd_sort,
                stdout=sbp.DEVNULL,
                stderr=sbp.DEVNULL,
                shell=True,
                check=True)

        sbp.run(cmd_rmdup,
                stdout=sbp.DEVNULL,
                stderr=sbp.DEVNULL,
                shell=True,
                check=True)

        sbp.run(cmd_index,
                stdout=sbp.DEVNULL,
                stderr=sbp.DEVNULL,
                shell=True,
                check=True)

        sbp.run(cmd_coverage,
                stdout=sbp.DEVNULL,
                stderr=sbp.DEVNULL,
                shell=True,
                check=True)

        sbp.run(cmd_rm,
                stdout=sbp.DEVNULL,
                stderr=sbp.DEVNULL,
                shell=True,
                check=True)

        print(time_stamp(),
              'alignment successfully finished.',
              flush=True)
