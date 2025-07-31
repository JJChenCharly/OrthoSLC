#!/bin/bash

# Check if g++ is installed
if command -v g++ > /dev/null; then
    echo "g++ is installed"
    echo "Compiling..."
else
    echo "g++ is not installed"
    echo "For ubuntu, you may try: "
    echo "'apt install build-essential'"
    echo "'apt-get update'"
    echo "'apt-get install g++'"
    echo "For centos, you can try:"
    echo "'yum update'"
    echo "'sudo yum install gcc-c++'"
    exit
fi

util_path=$1"/Utils.cpp"

cp $1"/Step1_preWalk.py" $2"/Step1_preWalk.py" && \
g++ -std=c++17 -pthread \
$1"/Step2_simple_derep.cpp" \
$util_path \
-o $2"/Step2_simple_derep" -O2 && \
g++ -std=c++17 -pthread \
$1"/Step3_pre_cluster.cpp" \
$util_path \
-o $2"/Step3_pre_cluster" -O2 && \
cp $1"/OrthoSLC.py" $2"/OrthoSLC.py" && \
cp $1"/Step4_makeblastdb.py" $2"/Step4_makeblastdb.py" && \
cp $1"/Step4_reciprocal_blast.py" $2"/Step4_reciprocal_blast.py" && \
g++ -std=c++17 -pthread \
$1"/Step5_query_binning.cpp" \
-o $2"/Step5_query_binning" -O2 && \
g++ -std=c++17 -pthread \
$1"/Step6_filter_n_bin.cpp" \
-o $2"/Step6_filter_n_bin" -O2 && \
g++ -std=c++17 -pthread \
$1"/Step7_RBF.cpp" \
-o $2"/Step7_RBF" -O2 && \
cp $1"/Step8_MCL.py" $2"/Step8_MCL.py" && \
g++ -std=c++17 \
$1/"cluster_fusion.cpp" \
-o $2"/cluster_fusion" -O2 && \
g++ -std=c++17 -pthread \
$1"/Step8_SLC.cpp" \
-o $2"/Step8_SLC" -O2 && \
g++ -std=c++17 -pthread \
$1"/Step9_write_clusters.cpp" \
$util_path \
-o $2"/Step9_write_clusters" -O2 && \
cp $1"/TK_kalign.py" $2"/TK_kalign.py" && \
cp $1"/TK_mafft.py" $2"/TK_mafft.py" && \
cp $1"/TK_AlnConcat.py" $2"/TK_AlnConcat.py" && \
cp $1"/TK_SNPmat.py" $2"/TK_SNPmat.py" && \
echo "Done"