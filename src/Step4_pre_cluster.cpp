#include "ThreadPool.h"

#include <iostream>
#include <unordered_map>
#include <fstream>
#include <string>
#include <filesystem>

namespace fs = std::filesystem;

int main(int argc, char** argv) {
    // parameter parsing
    std::string cated_fasta_pth;
    std::string dereped_cated_fasta_pth;
    std::string nr_dir;
    std::string len_info_path;
    std::string pre_cluster_path;
    int process_num = 1;

    for (int i = 1; i < argc; ++i)
    {
        if (std::string(argv[i]) == "--input_path" || std::string(argv[i]) == "-i")
        {
            cated_fasta_pth = argv[i + 1];
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
            std::cout << "Thanks for using OrthoSLC! (version: 0.1)\n\n";
            std::cout << "Usage: Step4_pre_cluster -i concatenated.fasta -d dereped.fasta -n nr_genome/ -l seq_len.txt -p pre_cluster.txt [options...]\n\n";
            std::cout << "  -i or --input_path ---------> path/to/input.FASTA\n";
            std::cout << "  -d or --concatenated_fasta -> path/to/output/dereplicated FASTA\n";
            std::cout << "  -n or --nr_genomes ---------> path/to/directory/of/output/non-reundant/genomes\n";
            std::cout << "  -l or --seq_len_info -------> path/to/output/sequence_length_table\n";
            std::cout << "  -p or --pre_cluster --------> path/to/output/pre-clustered_file\n";
            std::cout << "  -u or --thread_number ------> thread number, default: 1\n";
            std::cout << "  -h or --help ---------------> display this information\n";
            exit(0);
        }
    }
    // if file path exist
    if (!(fs::exists(cated_fasta_pth))) {
        std::cerr << "Error: path provided to '-i or --input_path' do not exist. 路径不存在\n";
        exit(0);
    }

    // read cated fasta ----
    std::ifstream cated_fasta(cated_fasta_pth);
    std::ofstream len_info(len_info_path, std::ios::trunc);
    std::ofstream dereped_cated_fasta(dereped_cated_fasta_pth, std::ios::trunc);
    
    typedef std::unordered_map<std::string, std::vector<std::string>> s_v;
    typedef std::unordered_map<std::string, std::string> s_s;
    s_v pre_cluster, each_d; // seq: short_id

    s_s dereped_cated_dict;

    std::string a_line, seq_id, original_id;
    std::string a_seq = "";

    getline (cated_fasta, a_line);
    original_id = a_line.substr(1);
    seq_id = a_line.substr(1, 11);

    while (getline (cated_fasta, a_line)) {
                
        if (a_line[0] == '>') { // if the row of id 

            // operate on previous SeqRecord using previous seq_id
            if (pre_cluster.count(a_seq) == 0) { // if no prescence before, need comparison with find
                pre_cluster[a_seq].push_back(seq_id);

                len_info << seq_id + "\t" << a_seq.length() << "\n";

                std::string rec_id = ">" + original_id + "\n";
            
                dereped_cated_fasta << rec_id;

                std::string a_seq_to_write;
                for (int saver = 0; saver < a_seq.length(); saver += 60) {
                    a_seq_to_write = a_seq_to_write + a_seq.substr(saver, 60) + "\n";
                }

                dereped_cated_fasta << a_seq_to_write;

                dereped_cated_dict[seq_id] = rec_id + a_seq_to_write;

                each_d[seq_id.substr(0, 5)].push_back(seq_id);

            } else { // if presence before
                (pre_cluster.at(a_seq)).push_back(seq_id);
                
            }
            
            // cover previous seq_id
            original_id = a_line.substr(1);
            seq_id = a_line.substr(1, 11);
            
            // cover previous a_seq
            a_seq = ""; // reset a_seq
        } else {
            a_seq = a_seq + a_line;
        }
    }

    // last rec
    if (pre_cluster.count(a_seq) == 0) {
        pre_cluster[a_seq].push_back(seq_id);

        len_info << seq_id + "\t" << a_seq.length() << "\n";

        std::string rec_id = ">" + original_id + "\n";
    
        dereped_cated_fasta << rec_id;

        std::string a_seq_to_write;
        for (int saver = 0; saver < a_seq.length(); saver += 60) {
            a_seq_to_write = a_seq_to_write + a_seq.substr(saver, 60) + "\n";
        }

        dereped_cated_fasta << a_seq_to_write;

        dereped_cated_dict[seq_id] = rec_id + a_seq_to_write;

        each_d[seq_id.substr(0, 5)].push_back(seq_id);

    } else {
        (pre_cluster.at(a_seq)).push_back(seq_id);
        
    }

    cated_fasta.close();
    len_info.close();
    dereped_cated_fasta.close();

    // write pre-cluster ----
    std::ofstream pre_cluster_f(pre_cluster_path, std::ios::trunc);

    for (const auto& [key, value] : pre_cluster) {

        int LEN = static_cast<int>(value.size());
        int LEN_1 = LEN - 1;

        for (int i = 0; i <= LEN_1; i++) {
            if (i == LEN_1) {
                pre_cluster_f << value[i] + "\n";
            } else {
                pre_cluster_f << value[i] + "\t";
            }
        }
    }

    pre_cluster_f.close();

    // nr genomes ----
    // if op path exit
    if (!(fs::exists(nr_dir))){
        fs::create_directory(nr_dir);
    }
    
    // mt ----
    ThreadPool pool(process_num);

    for (const auto& [key, value] : each_d) {
        s_s* dereped_cated_dict_ptr = &dereped_cated_dict;
        std::string* op = &nr_dir;

        pool.enqueue([key, op, dereped_cated_dict_ptr, value] {
            std::string op_full_path = fs::path(*op) / fs::path(key + ".fasta");
            std::ofstream a_nr_genome(op_full_path, std::ios::trunc);

            for(const auto& v : value)
            {
                a_nr_genome << dereped_cated_dict_ptr->at(v);
            }

            a_nr_genome.close();
        });

    }


    return 0;
}