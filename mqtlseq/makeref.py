from mqtlseq.utils import time_stamp, clean_cmd, call_log

class Makeref(object):

    def __init__(self,output_dir_aln,target_txt,p1_name,sym_ref):
        self.output_dir = output_dir_aln
        self.target_txt = target_txt
        self.p1_name = p1_name
        self.sym_ref = sym_ref

    def run(self):

        print(time_stamp(),
              'start to make {0} referance sequence.'.format(self.p1_name),
              flush=True)

        mydict = {}
        txt = open(self.target_txt, "r")
        txt_line = txt.readline()

        while txt_line:

            txt_line=txt_line.replace('\n', '')
            colom = txt_line.split()

            chr=colom[0]
            posi=colom[1]
            ref_base=colom[2]
            mut_base=colom[3]

            chr_posi=chr+":::"+posi
            base=ref_base+":::"+mut_base
            # print(chr_posi,base)
            mydict[chr_posi] = base

            txt_line = txt.readline()

        txt.close()

        ref = open(self.sym_ref, "r")
        output_file = open("{0}/{1}.fa".format(self.output_dir,self.p1_name), "w")
        ref_line = ref.readline()

        base_count=0
        stock_chr=""

        while ref_line:

            ref_line=ref_line.replace('\n', '')

            if ref_line.startswith('>'):

                output_file.write("{}\n".format(ref_line))
                colom = ref_line.split()
                chr=colom[0].replace('>', '')
                base_count=0
                stock_chr=chr

            else:

                colom=list(ref_line)
                array=""

                for item in colom:
                    base_count=base_count+1
                    chr_posi=stock_chr+":::"+str(base_count)

                    # print(chr_posi)
                    if chr_posi in mydict:
                        base = mydict[chr_posi].split(":::")
                        array=array+base[1]
                        # print("{0}\t{1}\t{2}\t{3}".format(stock_chr,base_count,item,base[0],base[1]))
                    else:
                        array=array+item

                output_file.write("{}\n".format(array))

            ref_line = ref.readline()

        ref.close()

        print(time_stamp(),
              'make {0} referance sequence successfully finished.'.format(self.p1_name),
              flush=True)
