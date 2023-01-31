# cml template
# install ----
chmod a+x install.sh

time ./install.sh src/ bins/

# bin directory ----
bin_dir=""
# dir for output
cp=""

# Step 1 ----
time $bin_dir/Step1_preparation \
-i input/dir \
-o $cp"/Step1_op.txt"

# Step 2 ----
time $bin_dir/Step2_simple_derep \
-i $cp"/Step1_op.txt" \
-o $cp"/S2_op_dereped" \
-u 16

# Step 3 ----
time $bin_dir/Step3_seq_concat \
-i $cp"/S2_op_dereped/" \
-c $cp"/S3_op_cated.fasta"

# Step4 ----
time $bin_dir/Step4_pre_cluster \
-i $cp"/S3_op_cated.fasta" \
-d $cp"/S4_op_dereped_cated.fasta" \
-n $cp"/S4_op_nr_genomes" \
-l $cp"/S4_op_seq_len.txt" \
-p $cp"/S4_op_pre_cluster.txt" \
-u 16

# Step 5 ----
time python3 $bin_dir/Step5_makeblastdb.py \
-c makeblastdb \
-i $cp"/S4_op_nr_genomes" \
-o $cp"/S5_op_dbs" \
-u 16 \
-t nucl

time python3 $bin_dir/Step5_reciprocal_blast.py \
-c blastn \
-i $cp"/S4_op_dereped_cated.fasta" \
-d $cp"/S5_op_dbs" \
-o $cp"/S5_op_blast_op" \
-e 1e-5 \
-u 16 \
-m off

# Step 6 ----
time $bin_dir/Step6_query_binning \
-i $cp"/S5_op_blast_op" \
-o $cp"/S6_op" \
-u 16 \
-L 64 \
-k off

# Step 7 ----
time $bin_dir/Step7_filter_n_bin \
-i $cp"/S6_op" \
-o $cp"/S7_op" \
-s $cp"/S4_op_seq_len.txt" \
-p $cp"/S4_op_pre_cluster.txt" \
-L 48 -r 0.3 -u 16 -k off

# Step 8 ----
time $bin_dir/Step8_RBF \
-i $cp"/S7_op" \
-o $cp"/S8_op" \
-u 16 -L 32 -k off

# Step 9 ----
time $bin_dir/Step9_SLC \
-i $cp"/S8_op" \
-o $cp"/SLC_1" \
-u 16 \
-S 1

time $bin_dir/Step9_SLC \
-i $cp"/SLC_1" \
-o $cp"/SLC_final" \
-p $cp"/S4_op_pre_cluster.txt" \
-u 16 \
-S all

# Step 10, here is an example for case where 200 genomes participates analysis, the amount left in $cp"/S2_op_dereped" ----
time $bin_dir/Step10_write_clusters \
-i $cp"/SLC_final/0.txt" \
-o $cp"/S10_write_fasta" \
-f $cp"/S3_op_cated.fasta" \
-c 200 \
-t strict,surplus,accessory \
-u 16
