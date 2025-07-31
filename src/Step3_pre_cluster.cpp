#define WRITER_BUF 2 // MB
#define KEEP_SHORT_ID false // MB
#define MAP_RESRVE_COEF 1250  //5000 * 0.25 * amont of files

#include "ThreadPool.h"
#include "Utils.hpp"

#include <iostream>
#include <string>
#include <memory>
#include <filesystem>

namespace fs = std::filesystem;
using namespace fasta;

int main(int argc, char** argv) {
    // parameter parsing
    std::string dereped_dir;
    std::string dereped_cated_fasta_pth;
    std::string nr_dir;
    std::string len_info_path;
    std::string id_info_path;
    std::string pre_cluster_path;
    int process_num = 1;
    bool keyInclude = KEEP_SHORT_ID;

    for (int i = 1; i < argc; ++i)
    {
        if (std::string(argv[i]) == "--input_path" || std::string(argv[i]) == "-i")
        {
            dereped_dir = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--dereped" || std::string(argv[i]) == "-d")
        {
            dereped_cated_fasta_pth = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--nr_genomes" || std::string(argv[i]) == "-n")
        {
            nr_dir = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--seq_len_info" || std::string(argv[i]) == "-l")
        {
            len_info_path = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--id_info" || std::string(argv[i]) == "-m")
        {
            id_info_path = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--pre_cluster" || std::string(argv[i]) == "-p")
        {
            pre_cluster_path = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--thread_number" || std::string(argv[i]) == "-u")
        {
            process_num = std::stoi(argv[i + 1]);
        }
        else if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h")
        {
            std::cout << "Thanks for using OrthoSLC! (version: " << __version__ << ")\n\n";
            std::cout << "Usage: Step3_pre_cluster -i concatenated.fasta -d dereped/ -n nr_genome/ -l seq_len.txt -p pre_cluster.txt [options...]\n\n";
            std::cout << "  -i or --input_path -----> <dir> path/to/directory of dereplicated FASTA from Step2\n";
            std::cout << "  -d or --derep_fasta ----> <fasta> path/to/output/dereplicated concatenated FASTA\n";
            std::cout << "  -n or --nr_genomes -----> <dir> path/to/directory/of/output/non-reundant/genomes\n";
            std::cout << "  -l or --seq_len_info ---> <txt> path/to/output/sequence_length_table\n";
            std::cout << "  -m or --id_info --------> <txt> path/to/output/id_info_table\n";
            std::cout << "  -p or --pre_cluster ----> <txt> path/to/output/pre_clustered_file\n";
            std::cout << "  -u or --thread_number --> <int> thread number, default: 1\n";
            std::cout << "  -h or --help -----------> display this information\n";
            exit(0);
        }
    }
    // if file path exist
    if (!(fs::exists(dereped_dir))) {
        std::cerr << "Error: path provided to '-i or --input_path' does not exist. 路径不存在\n";
        exit(0);
    }

    fs::path parent_dereped_cated_fasta_pth = fs::path(dereped_cated_fasta_pth).parent_path();
    fs::path parent_nr_dir = fs::path(nr_dir).parent_path();
    fs::path parent_len_info_path = fs::path(len_info_path).parent_path();
    fs::path parent_pre_cluster_path = fs::path(pre_cluster_path).parent_path();
    fs::path parent_id_info_path = fs::path(id_info_path).parent_path();

    if (!(fs::exists(parent_dereped_cated_fasta_pth))) {
        std::cerr << "Error: parent path provided to '-d or --concatenated_fasta' does not exist. 路径不存在\n";
        exit(0);
    } 
    else if (!(fs::exists(parent_nr_dir))) {
        std::cerr << "Error: parent path provided to '-n or --nr_genomes' does not exist. 路径不存在\n";
        exit(0);
    }
    else if (!(fs::exists(parent_len_info_path))) {
        std::cerr << "Error: parent path provided to '-l or --seq_len_info' does not exist. 路径不存在\n";
        exit(0);
    }
    else if (!(fs::exists(parent_pre_cluster_path))) {
        std::cerr << "Error: parent path provided to '-p or --pre_cluster' does not exist. 路径不存在\n";
        exit(0);
    }
    else if (!(fs::exists(parent_id_info_path))) {
        std::cerr << "Error: parent path provided to '-m or --id_info' does not exist. 路径不存在\n";
        exit(0);
    }

    // read in each fasta and dereplicate all, get seq_len, id_info ----

    // 统计文件数量
    int file_count = 0;
    for (const auto& entry : fs::directory_iterator(dereped_dir)) {
        if (entry.is_regular_file()) {
            file_count++;
        }
    }
    
    Writer cated_fasta_writer(dereped_cated_fasta_pth, WRITER_BUF * 1024 * 1024);
    DeduplicatorExact dedup;  // or DeduplicatorExact for perfect identity
    dedup.reserve(MAP_RESRVE_COEF*file_count);

    std::unordered_map<std::string, std::unordered_map<std::string, const std::string*>> each_spe_d;
    
    std::string current_path;
    for (const auto& entry : fs::directory_iterator(dereped_dir)) {
        
        current_path = entry.path();
        current_path = fs::absolute(current_path);
        
        Reader reader(current_path);
        std::string header, seq;
        while (reader.next(header, seq)) {
            size_t fisrt_space = header.find(' ');
            const std::string id = header.substr(0, fisrt_space);

            // write original info for each input
            tsv::id_info_TSVwriter(id_info_path, id, header.substr(fisrt_space+1), 1 * 1024 * 1024);

            auto [ seq_ptr, is_new ] = dedup.add_and_get(seq, id);
            if (is_new) {
                // the dereplicated concatenated fasta
                cated_fasta_writer.write(id, seq); 

                // write len for dereplicated ones
                tsv::len_info_TSVwriter(len_info_path, id, seq.length(), 1 * 1024 * 1024);

                std::string species = id.substr(0, id.find('-'));
                each_spe_d[species].emplace(id, seq_ptr);

            }
        }
    }
    cated_fasta_writer.flush();

    // nr genomes ----
    // if op path exit
    if (!(fs::exists(nr_dir))){
        fs::create_directory(nr_dir);
    }
    
    // mt ----
    ThreadPool pool(process_num);
    
    std::vector<std::future<void>> futures;
    
    // write pre-cluster ----
    void (*TSVwriterFp)(const std::string&, const tsv::Map&) 
    = keyInclude 
        ? &tsv::writeTSVWithKey 
        : &tsv::writeTSV;

    const auto& hashMap = dedup.map();
    futures.emplace_back(
        pool.enqueue(
            [&pre_cluster_path, &hashMap, &TSVwriterFp](){
                TSVwriterFp(pre_cluster_path, hashMap);
                return;
            }
        )
    );
    
    std::string out_dir = nr_dir;
    for (const auto& [species, id_map] : each_spe_d) {
        futures.emplace_back(
            pool.enqueue(
                [ species, id_map, out_dir ]() {
                    // build the filename: <nr_dir>/<species>.fasta
                    std::string out_path = fs::path(out_dir) / fs::path( species + ".fasta");
                    // open a fasta::Writer with the same buffer size macro
                    Writer writer(out_path, WRITER_BUF * 1024 * 1024);

                    // write each <id, seq> pair
                    for (const auto& [id, seq_ptr] : id_map) {
                        // seq_ptr is const std::string* from your dedup map
                        writer.write(id, *seq_ptr);
                    }

                    // flush before destructor (optional—destructor will also flush)
                    writer.flush();
                }
            )
        );
        
    }

    for (auto& future : futures) {
        future.wait();
    }
    return 0;
}