import subprocess as sbp
from mqtlseq.utils import time_stamp, clean_cmd, call_log
import pandas as pd

class Mpileup(object):

    def __init__(self,output_dir_aln,bam_list,ref,thread,output_name,base_quality,mapping_quality):
        self.output_dir = output_dir_aln
        self.ref = ref
        self.bam_table = bam_list
        self.thread = thread
        self.output_name = output_name
        self.base_quality = base_quality
        self.mapping_quality = mapping_quality

    def run(self):
        print(time_stamp(),
              'start to mpileup by bcftools.',
              flush=True)

        print(time_stamp(),
              self.bam_table,
              flush=True)

        output_name=self.output_name
        rm_list = ""

        fai = open("{0}.fai".format(self.ref), "r")
        df=pd.read_table(fai,header=None)
        chr_item = df[0]

        command_mpile = open(self.output_dir+"/mpileup.txt", "w")
        command_call = open(self.output_dir+"/call.txt", "w")
        command_fil = open(self.output_dir+"/filter.txt", "w")
        command_rm = open(self.output_dir+"/rm.txt", "w")


        for chr in chr_item:

            mpile_format = 'bcftools mpileup -r {0} -a AD,DP -B -q {5} -Q {6} -O u -f {1} --ignore-RG {2} 1>{3}/{0}.{4}.bcf 2>>{3}/../log/bcftools.log\n'.format(chr,
                                                                                                                                                         self.ref,
                                                                                                                                                         self.bam_table,
                                                                                                                                                         self.output_dir,
                                                                                                                                                         output_name,
                                                                                                                                                         self.mapping_quality,
                                                                                                                                                         self.base_quality)

            call_format = 'bcftools call -vm -O u {0}/{2}.{1}.bcf 1>{0}/{2}.{1}.call.bcf 2>>{0}/../log/bcftools.log\n'.format(self.output_dir,output_name,chr)

            filter_format = 'bcftools filter -O v -o {0}/{2}.{1}.vcf {0}/{2}.{1}.call.bcf 2>>{0}/../log/bcftools.log\n'.format(self.output_dir,output_name,chr)

            rm1_format = 'rm {0}/{2}.{1}.bcf\n'.format(self.output_dir,output_name,chr)
            rm2_format = 'rm {0}/{2}.{1}.call.bcf\n'.format(self.output_dir,output_name,chr)
            rm3_format = 'rm {0}/{2}.{1}.vcf\n'.format(self.output_dir,output_name,chr)

            command_mpile.write(mpile_format)
            command_call.write(call_format)
            command_fil.write(filter_format)
            command_rm.write(rm1_format)
            command_rm.write(rm2_format)
            command_rm.write(rm3_format)

        command_mpile.close()
        command_call.close()
        command_fil.close()
        command_rm.close()

        cmd_mpile = 'cat {0}/mpileup.txt|xargs -P{1} -I % sh -c %'.format(self.output_dir,self.thread)
        cmd_mpile = clean_cmd(cmd_mpile)

        cmd_call = 'cat {0}/call.txt|xargs -P{1} -I % sh -c %'.format(self.output_dir,self.thread)
        cmd_call = clean_cmd(cmd_call)

        cmd_fil = 'cat {0}/filter.txt|xargs -P{1} -I % sh -c %'.format(self.output_dir,self.thread)
        cmd_fil = clean_cmd(cmd_fil)

        cmd_rm = 'cat {0}/rm.txt|xargs -P{1} -I % sh -c %'.format(self.output_dir,self.thread)
        cmd_rm = clean_cmd(cmd_rm)


        sbp.run(cmd_mpile,
                shell=True,
                check=True)

        sbp.run(cmd_call,
                shell=True,
                check=True)

        sbp.run(cmd_fil,
                shell=True,
                check=True)

        all_vcf = open('{0}/{1}.vcf'.format(self.output_dir,output_name), "w")

        for chr in chr_item:

            chr_vcf = open("{0}/{2}.{1}.vcf".format(self.output_dir,output_name,chr), "r")
            chr_vcf_line = chr_vcf.readline()
            # print("{0}/{2}.{1}.vcf".format(self.output_dir,output_name,chr))

            while chr_vcf_line:

                chr_vcf_line=chr_vcf_line.replace('\n', '')

                if chr_vcf_line.startswith('#'):
                    aaa="aaa"
                else:
                    all_vcf.write('{0}\n'.format(chr_vcf_line))

                chr_vcf_line = chr_vcf.readline()

            chr_vcf.close()

        sbp.run(cmd_rm,
                shell=True,
                check=True)

        print(time_stamp(),
              'mpileup successfully finished.',
              flush=True)
