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

g++ -std=c++17 \
$1"/Step1_preparation.cpp" \
-o $2"/Step1_preparation" && \
g++ -std=c++17 -pthread \
$1"/Step2_simple_derep.cpp" \
-o $2"/Step2_simple_derep" && \
g++ -std=c++17 -pthread \
$1"/Step3_pre_cluster.cpp" \
-o $2"/Step3_pre_cluster" && \
cp $1"/Blast.py" $2"/Blast.py" && \
cp $1"/Step4_makeblastdb.py" $2"/Step4_makeblastdb.py" && \
cp $1"/Step4_reciprocal_blast.py" $2"/Step4_reciprocal_blast.py" && \
g++ -std=c++17 -pthread \
$1"/Step5_query_binning.cpp" \
-o $2"/Step5_query_binning" && \
g++ -std=c++17 -pthread \
$1"/Step6_filter_n_bin.cpp" \
-o $2"/Step6_filter_n_bin" && \
g++ -std=c++17 -pthread \
$1"/Step7_RBF.cpp" \
-o $2"/Step7_RBF" && \
g++ -std=c++17 -pthread \
$1"/Step8_SLC.cpp" \
-o $2"/Step8_SLC" && \
g++ -std=c++17 -pthread \
$1"/Step9_write_clusters.cpp" \
-o $2"/Step9_write_clusters" && \
echo "Done"