# cml template----

# by bash ----
bin_dir=$1
wd=$2
in_dir=$3
mkdir -p $wd
cpu=$4
bin_num=$5
in_ext=$6

# by manual ----
# bin_dir="./bins"
# wd="./test_op"
# in_dir="./test_inputs"
# mkdir $wd
# cpu="16"
# bin_num='64'

# Step 1 ----
time python3 $bin_dir/Step1_preWalk.py \
-i $in_dir \
-o $wd"/Step1_op.txt" \
-f $in_ext

# Step 2 ----
time $bin_dir/Step2_simple_derep \
-i $wd"/Step1_op.txt" \
-o $wd"/S2_op_dereped" \
-c $wd"/S2_op_copy_info.txt" \
-u $cpu
# remove very small files output by step 2 if necessary !!!!!!!!!!!!!!!!!!!
# grep -c ">" $wd/S2_op_dereped/* | sort -t: -k2 # 小到大
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
-f "6 qseqid sseqid pident score evalue" \
-u $cpu

# Step 5 ----
time $bin_dir/Step5_query_binning \
-i $wd"/S4_op_blast_op" \
-o $wd"/S5_op" \
-u $cpu \
-L $bin_num \
-k off

# Step 6 ----
time $bin_dir/Step6_filter_n_bin \
-i $wd"/S5_op" \
-o $wd"/S6_op" \
-s $wd"/S3_op_seq_len.txt" \
-p $wd"/S3_op_pre_cluster.txt" \
-L $bin_num -r 0.8 -m 80 -w bitscore -u $cpu -k off

# Step 7 ----
time $bin_dir/Step7_RBF \
-i $wd"/S6_op" \
-o $wd"/S7_op" \
-u $cpu -L $bin_num -k off

# Step 8 ----
########### option 1: single linkage clustering (select -w 'none' in Step 6) ###########
# time $bin_dir/Step8_SLC \
# -i $wd"/S7_op" \
# -o $wd"/SLC_1" \
# -u $cpu \
# -S 5

# time $bin_dir/Step8_SLC \
# -i $wd"/SLC_1" \
# -o $wd"/Final_cluster" \
# -p $wd"/S3_op_pre_cluster.txt" \
# -u $cpu \
# -S all

########### option 2: MCL ###########
time python3 $bin_dir/Step8_MCL.py \
-c mcl \
-i $wd"/S7_op" \
-p $wd"/S3_op_pre_cluster.txt" \
-o $wd"/Final_cluster" \
-u $cpu \
-I 1.5 --abc

# Step 9 ----
time $bin_dir/Step9_write_clusters \
-i $wd"/Final_cluster/0.txt" \
-o $wd"/S9_write_fasta" \
-f $wd"/S3_op_dereped_cated.fasta" \
-m $wd"/S3_op_id_info.txt" \
-p $wd"/S3_op_pre_cluster.txt" \
-a 0 \
-c `find $wd"/S2_op_dereped" -type f | wc -l` \
-t strict,surplus,accessory \
-u $cpu

# clean unnecessaries (easy re-run with parameter change) if you want ----
for dir in "$wd/S5_op" "$wd/S6_op" "$wd/S7_op" "$wd/S2_op_dereped" "$wd/S3_op_nr_genomes"; do
    if [ -d "$dir" ]; then
        rm -r "$dir"
    fi
done

for slc_dir in "$wd"/SLC_*; do
    if [ -d "$slc_dir" ]; then
        rm -r "$slc_dir"
    fi
done
