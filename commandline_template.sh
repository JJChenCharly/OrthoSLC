# cml template
# install ----
chmod a+x install.sh

time ./install.sh src/ bins/

# bin directory ----
bin_dir="./bins"
# dir for output
wd="./test_op"
mkdir $wd
# thread number
cpu="36"

# Step 1 ----
time $bin_dir/Step1_preparation \
-i test_inputs/ \
-o $wd"/Step1_op.txt"

# Step 2 ----
time $bin_dir/Step2_simple_derep \
-i $wd"/Step1_op.txt" \
-o $wd"/S2_op_dereped" \
-u $cpu
# remove very small files output by step 2 if necessary !!!!!!!!!!!!!!!!!!!
grep -c ">" $wd/S2_op_dereped/* | sort -t: -k2 # 小到大
# grep -c ">" $wd/S2_op_dereped/* | sort -t: -k2rn # 大到小

# Step 3 ----
time $bin_dir/Step3_pre_cluster \
-i $wd"/S2_op_dereped" \
-d $wd"/S3_op_dereped_cated.fasta" \
-n $wd"/S3_op_nr_genomes" \
-l $wd"/S3_op_seq_len.txt" \
-m $wd"/S3_op_id_info.txt" \
-p $wd"/S3_op_pre_cluster.txt" \
-u $cpu

# Step 4 ----
time python3 $bin_dir/Step4_makeblastdb.py \
-c makeblastdb \
-i $wd"/S3_op_nr_genomes" \
-o $wd"/S4_op_dbs" \
-u $cpu \
-t nucl

time python3 $bin_dir/Step4_reciprocal_blast.py \
-c blastn \
-i $wd"/S3_op_dereped_cated.fasta" \
-d $wd"/S4_op_dbs" \
-o $wd"/S4_op_blast_op" \
-e 1e-5 \
-u $cpu \
-m off

# Step 5 ----
time $bin_dir/Step5_query_binning \
-i $wd"/S4_op_blast_op" \
-o $wd"/S5_op" \
-u $cpu \
-L 64 \
-k off

# Step 6 ----
time $bin_dir/Step6_filter_n_bin \
-i $wd"/S5_op" \
-o $wd"/S6_op" \
-s $wd"/S3_op_seq_len.txt" \
-p $wd"/S3_op_pre_cluster.txt" \
-L 48 -r 0.3 -u $cpu -k off

# Step 7 ----
time $bin_dir/Step7_RBF \
-i $wd"/S6_op" \
-o $wd"/S7_op" \
-u $cpu -L 36 -k off

# Step 8 ----
time $bin_dir/Step8_SLC \
-i $wd"/S7_op" \
-o $wd"/SLC_1" \
-u $cpu \
-S 1

time $bin_dir/Step8_SLC \
-i $wd"/SLC_1" \
-o $wd"/SLC_final" \
-p $wd"/S3_op_pre_cluster.txt" \
-u $cpu \
-S all

# Step 10, here is an example for case where 10 genomes participates analysis, the amount left in $wd"/S2_op_dereped" ----
# check how many participated
ls $wd"/S2_op_dereped/" | cat -n | tail

time $bin_dir/Step9_write_clusters \
-i $wd"/SLC_final/0.txt" \
-o $wd"/S9_write_fasta" \
-f $wd"/S3_op_dereped_cated.fasta" \
-m $wd"/S3_op_id_info.txt" \
-p $wd"/S3_op_pre_cluster.txt" \
-c 16 \
-t strict,surplus,accessory \
-u $cpu
