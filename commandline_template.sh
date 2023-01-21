# cml template
# install ----
chmod a+x install.sh

time ./install.sh src/ bins/

# working directory ----
cp="/working_dir"

# Step 1 ----
time ./bins/Step1_preparation \
-i $cp"/dir_input_annotated/" \
-o $cp"/test_op/Step1_pre_res.txt"

# Step 2 ----
time ./bins/Step2_simple_derep \
-i $cp"/test_op/Step1_pre_res.txt" \
-o $cp"/test_op/dereped" \
-u 16

# Step 3 ----
time ./bins/Step3_seq_concat \
-i $cp"/test_op/dereped/" \
-c $cp"/test_op/cated.fasta"

# Step4 ----
time ./bins/Step4_pre_cluster \
-i $cp"/test_op/cated.fasta" \
-d $cp"/test_op/dereped_cated.fasta" \
-n $cp"/test_op/nr_genomes" \
-l $cp"/test_op/seq_len.txt" \
-p $cp"/test_op/pre_cluster.txt" \
-u 16

# Step 5 ----
time python3 ./bins/Step5_makeblastdb.py \
-c makeblastdb \
-i $cp"/test_op/nr_genomes" \
-o $cp"/test_op/dbs" \
-u 16 \
-t nucl

time python3 ./bins/Step5_reciprocal_blast.py \
-c blastn \
-i $cp"/test_op/dereped_cated.fasta" \
-d $cp"/test_op/dbs" \
-o $cp"/test_op/blast_op" \
-e 1e-5 \
-u 14 \
-m off

# Step 6 ----
time ./bins/Step6_query_binning \
-i $cp"/test_op/blast_op" \
-o $cp"/test_op/query_binning" \
-u 16 \
-L 64 \
-k off

# Step 7 ----
time ./bins/Step7_filter_n_bin \
-i $cp"/test_op/query_binning" \
-o $cp"/test_op/bin_op" \
-s $cp"/test_op/seq_len.txt" \
-p $cp"/test_op/pre_cluster.txt" \
-L 48 -r 0.3 -u 16 -k off

# Step 8 ----
time ./bins/Step8_RBF \
-i $cp"/test_op/bin_op" \
-o $cp"/test_op/RBB_op" \
-u 16 -L 32 -k off

# Step 9 ----
time ./bins/Step9_SLC \
-i $cp"/test_op/RBB_op" \
-o $cp"/test_op/SLC_1" \
-u 16 \
-S 1

time ./bins/Step9_SLC \
-i $cp"/test_op/SLC_1" \
-o $cp"/test_op/SLC_final" \
-p $cp"/test_op/pre_cluster.txt" \
-u 16 \
-S all

# Step 10, 200 genomes participates analysis ----
time ./bins/Step10_write_clusters \
-i $cp"/test_op/SLC_final/0.txt" \
-o $cp"/test_op/write_fasta" \
-f $cp"/test_op/cated.fasta" \
-c 200 \
-t strict,surplus,accessory \
-u 16