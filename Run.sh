#! /bin/sh
#$ -S /bin/sh
#$ -cwd


export PATH=/lustre7/home/lustre3/segawa-tenta/anaconda3/envs/QTL-seq/bin:$PATH:$HOME/bin:/lustre7/home/lustre3/segawa-tenta/tools/Mark_maker


python script/qtlseq.py -r /lustre7/home/lustre3/segawa-tenta/fasta/m_Brapa_chiifu_v41_genome20230413.fasta \
                        -rt reads_PATH.txt\
                        -p1 AKA\
                        -p2 KAN\
                        -a F2whole\
                        -o QTL-seq\
                        -t 10\
                        -window_size 2000000\
                        -step_size 50000
<<COMMENTOUT
################################### example1 ###################################
python script/qtlseq.py -r /lustre7/home/lustre3/segawa-tenta/Other/221110_QTL-seq_test_data/test_refarence.fa \
                        -rt reads_PATH.txt\
                        -p1 AKA \
                        -p2 KAN \
                        -a R_bulk\
                        -b P_bulk\
                        -o Run1\
                        -t 1\
                        -window_size 5000\
                        -step_size 1000\
                        -min_plot 1\
                        -t 10\
                        
                        
                        

################################### option1 ###################################
# 必須
# -r : Referance fasta のパス
# -rt : reads_PATH.txtのパス

# 任意
# -p1 : P1 の名前 : P1
# -p2 : P2 の名前 : P2
# -a : A bulk の名前 : A_bulk
# -b : B bulk の名前 : B_bulk
# -o : OUTPUT したいディレクトリのパス
# -t : Thread : 1
# -bwa : アライナー "bwa" or "bwa-mem2" : bwa
# -mindepth : 設定値 <= depth 場合、SNP-index 計算に使用 : 10
# -maxdepth : 設定値 >= depth 場合、SNP-index 計算に使用 : 99
# -sim_replication : シミュレーションの反復数 : 10000
# -individual : bulk した個体数 : 20
# -generation : bulk した世代 "F2" or "RIL" or "BC1F1"。 "BC1F1" の場合、P1を戻し親にする : F2
# -window_size : Sliding window の window サイズ (bp) : 2000000
# -step_size : Sliding window の step サイズ (bp) : 50000
# -min_plot : Window 内の最小 SNP 数。設定値以下では、Sliding window 計算時にスキップされる : 10
# -filt_bias : 2バルク使用時のみ使用される。2バルクともで(設定値)以下、または、(1-設定値)以上のSNP-indexをとるポジションをフィルターする : 0.3

################################### example2 ###################################
python script/qtlseq_mpileup.py -r /lustre7/home/lustre3/segawa-tenta/Other/221109_QTL-seq_pipline/new/Run1/20_ref/AKA.fa \
                                -bt bams_PATH.txt\
                                -p1 AKA \
                                -p2 KAN \
                                -a R_bulk\
                                -b P_bulk\
                                -o Run2\
                                -t 1\
                                -window_size 5000\
                                -step_size 1000\
                                -min_plot 1\
                                -t 10\
                        

################################### option2 ###################################
# 必須
# -r : Referance fasta のパス fai ファイルも同じディレクトリに必須
# -bt : bams_PATH.txtのパス

# 任意
# -p1 : P1 の名前 : P1
# -p2 : P2 の名前 : P2
# -a : A bulk の名前 : A_bulk
# -b : B bulk の名前 : B_bulk
# -o : OUTPUT したいディレクトリのパス
# -t : Thread : 1
# -mindepth : 設定値 <= depth 場合、SNP-index 計算に使用 : 10
# -maxdepth : 設定値 >= depth 場合、SNP-index 計算に使用 : 99
# -sim_replication : シミュレーションの反復数 : 10000
# -individual : bulk した個体数 : 20
# -generation : bulk した世代 "F2" or "RIL" or "BC1F1"。 "BC1F1" の場合、P1を戻し親にする : F2
# -window_size : Sliding window の window サイズ (bp) : 2000000
# -step_size : Sliding window の step サイズ (bp) : 50000
# -min_plot : Window 内の最小 SNP 数。設定値以下では、Sliding window 計算時にスキップされる : 10
# -filt_bias : 2バルク使用時のみ使用される。2バルクともで(設定値)以下、または、(1-設定値)以上のSNP-indexをとるポジションをフィルターする : 0.3

################################### example3 ###################################
python script/qtlseq_snpindex.py -vcf /lustre7/home/lustre3/segawa-tenta/Other/221109_QTL-seq_pipline/new/Run2/22_vcf/AKA.KAN.F1.R_bulk.P_bulk.vcf\
                                 -p1 AKA \
                                 -p2 KAN \
                                 -f1 True \
                                 -a R_bulk\
                                 -b P_bulk\
                                 -o Run3\
                                 -t 1\
                                 -window_size 5000\
                                 -step_size 1000\
                                 -min_plot 1\
                                 -t 10\
                        

################################### option3 ###################################
# 必須
# -vcf : vcfのパス : "GT:PL:DP:AD"をP1,P2,F1,Abulk,Bbulkの順でコールしたvcfの絶対パス

# 任意
# -p1 : P1 の名前 : P1
# -p2 : P2 の名前 : P2があるなら必須
# -a : A bulk の名前 : A_bulk
# -b : B bulk の名前 : B_bulkがあるなら必須
# -f1 : f1があるならTrue
# -o : OUTPUT したいディレクトリのパス
# -t : Thread : 1
# -mindepth : 設定値 <= depth 場合、SNP-index 計算に使用 : 10
# -maxdepth : 設定値 >= depth 場合、SNP-index 計算に使用 : 99
# -sim_replication : シミュレーションの反復数 : 10000
# -individual : bulk した個体数 : 20
# -generation : bulk した世代 "F2" or "RIL" or "BC1F1"。 "BC1F1" の場合、P1を戻し親にする : F2
# -window_size : Sliding window の window サイズ (bp) : 2000000
# -step_size : Sliding window の step サイズ (bp) : 50000
# -min_plot : Window 内の最小 SNP 数。設定値以下では、Sliding window 計算時にスキップされる : 10
# -filt_bias : 2バルク使用時のみ使用される。2バルクともで(設定値)以下、または、(1-設定値)以上のSNP-indexをとるポジションをフィルターする : 0.3

################################### example4 ###################################
python script/qtlseq_plot.py -txt /lustre7/home/lustre3/segawa-tenta/Other/221109_QTL-seq_pipline/new/Run3/24_filt/AKA_homo.KAN_homo.F1_hetero.R_bulk.P_bulk.snp_index.txt\
                             -p2 True \
                             -f1 True \
                             -a R_bulk\
                             -b P_bulk\
                             -o Run4\
                             -t 1\
                             -window_size 5000\
                             -step_size 1000\
                             -min_plot 1\
                             -t 10\
                        

################################### option4 ###################################
# 必須
# -txt : txtのパス : snp-indexとdepthがフォーマット通りに並ぶ必要あり

# 任意
# -p2 : P2があるならTrue
# -a : A bulk の名前 : A_bulk
# -b : B bulk の名前 : B_bulkがあるなら必須
# -f1 : f1があるならTrue
# -o : OUTPUT したいディレクトリのパス
# -t : Thread : 1
# -mindepth : 設定値 <= depth 場合、SNP-index 計算に使用 : 10
# -maxdepth : 設定値 >= depth 場合、SNP-index 計算に使用 : 99
# -sim_replication : シミュレーションの反復数 : 10000
# -individual : bulk した個体数 : 20
# -generation : bulk した世代 "F2" or "RIL" or "BC1F1"。 "BC1F1" の場合、P1を戻し親にする : F2
# -window_size : Sliding window の window サイズ (bp) : 2000000
# -step_size : Sliding window の step サイズ (bp) : 50000
# -min_plot : Window 内の最小 SNP 数。設定値以下では、Sliding window 計算時にスキップされる : 10
# -filt_bias : 2バルク使用時のみ使用される。2バルクともで(設定値)以下、または、(1-設定値)以上のSNP-indexをとるポジションをフィルターする : 0.3


COMMENTOUT