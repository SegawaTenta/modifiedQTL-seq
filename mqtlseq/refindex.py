import sys
import subprocess as sbp
from mqtlseq.utils import time_stamp, clean_cmd, call_log


class RefIndex(object):

    def __init__(self, output_dir, ref_fa,bwa_type):
        self.output_dir = output_dir
        self.ref_fa = ref_fa
        self.bwa_type = bwa_type

    def run(self):
        print(time_stamp(),
              'start to index reference fasta.',
              flush=True)

        cmd1 = '{} index {} \
                >> {}/log/bwa.log \
                2>&1'.format(self.bwa_type,self.ref_fa, self.output_dir)

        cmd2 = 'samtools faidx {} \
                >> {}/log/samtools.log \
                2>&1'.format(self.ref_fa, self.output_dir)

        cmd1 = clean_cmd(cmd1)
        cmd2 = clean_cmd(cmd2)

        try:
            sbp.run(cmd1,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)
        except sbp.CalledProcessError:
            call_log(self.out, 'bwa', cmd1)
            sys.exit(1)

        try:
            sbp.run(cmd2,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)
        except sbp.CalledProcessError:
            call_log(self.output_dir, 'samtools', cmd2)
            sys.exit(1)

        print(time_stamp(),
              'indexing of reference successfully finished.',
              flush=True)
