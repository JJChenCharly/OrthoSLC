#include "ThreadPool.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_set>
#include <filesystem>

namespace fs = std::filesystem;

int main(int argc, char** argv) {
    // parameter parsing
    std::string strain_info_tsv_pth;
    std::string rr_ed_dir;
    int process_num = 1;

    for (int i = 1; i < argc; ++i)
    {
        if (std::string(argv[i]) == "--input_path" || std::string(argv[i]) == "-i")
        {
            strain_info_tsv_pth = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--output_path" || std::string(argv[i]) == "-o")
        {
            rr_ed_dir = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--thread_number" || std::string(argv[i]) == "-u")
        {
            process_num = std::stoi(argv[i + 1]);
        }
        else if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h")
        {
            std::cout << "Thanks for using OrthoSLC! (version: " << __version__ << ")\n\n";
            std::cout << "Usage: Step2_simple_derep -i input_file -o output/ [options...]\n\n";
            std::cout << "  -i or --input_path -----> <txt> path/to/file/output/by/Step1\n";
            std::cout << "  -o or --output_path ----> <dir> path/to/output/directory\n";
            std::cout << "  -u or --thread_number --> <int> thread number, default: 1\n";
            std::cout << "  -h or --help -----------> display this information\n";
            exit(0);
        }
    }
    // if file path exist
    if (!(fs::exists(strain_info_tsv_pth))) {
        std::cerr << "Error: path provided to '-i or --input_path' does not exist. 路径不存在\n";
        exit(0);
    }

    fs::path parent_rr_ed_dir = fs::path(rr_ed_dir).parent_path();
    if (!(fs::exists(parent_rr_ed_dir))) {
        std::cerr << "Error: parent path provided to '-o or --output_path' does not exist. 路径不存在\n";
        exit(0);
    }

    // if op path exit
    if (!(fs::exists(rr_ed_dir))){
        fs::create_directory(rr_ed_dir);
    }
    
    // get abs pth and short id ----
    std::string a_row;

    std::vector<std::string> pth_v, s_id_v;

    std::ifstream strain_in_put(strain_info_tsv_pth);

    std::string sid, naam, abs_p;
    while (getline(strain_in_put, a_row)) {

        std::stringstream ss(a_row);
        std::string key, value;
        std::getline(ss, sid, '\t');
        std::getline(ss, naam, '\t');
        std::getline(ss, abs_p);

        pth_v.push_back(abs_p);
        s_id_v.push_back(sid);
    }
    strain_in_put.close();

    int task_len = pth_v.size();

    // mt ----
    ThreadPool pool(process_num);

    for(int i = 0; i < task_len; ++i) {
        std::string in_p = pth_v[i];
        std::string s_id = s_id_v[i];
        std::string* op = &rr_ed_dir;
        pool.enqueue([in_p, s_id, op] { // regradless of number of input parameter
            // read in
            std::ifstream fasta_in_put(in_p);

            // for file to save
            std::string file_naam = fs::path(in_p).filename();
            std::string op_full_path = fs::path(*op) / fs::path(file_naam);
            std::ofstream fasta_out_put(op_full_path, std::ios::trunc);

            std::string a_line, info;
            std::string a_seq = "";
            int seq_id = 0; // start

            std::unordered_set<std::string> added_seq = {"added_seq"};

            // skip the first void seq write-in \n
            getline (fasta_in_put, a_line);
            info = a_line.substr(1);

            int saver;
            while (getline (fasta_in_put, a_line)) {
                
                if (a_line[0] == '>') { // if the row of id 

                    // operate on previous SeqRecord using previous seq_id
                    if (added_seq.count(a_seq) == 0) { // if no prescence before, need comparison with find
                        added_seq.insert(a_seq);

                        fasta_out_put << ">" + s_id + "-" << seq_id << " " + info + "\n"; // using previour seq_id
                        seq_id++;

                        // save 60 nchar a row
                        for (saver = 0; saver < a_seq.length(); saver += 60) {
                            fasta_out_put << a_seq.substr(saver, 60) + "\n";
                        }
                        // cover previous seq_id
                    }
                    
                    info = a_line.substr(1);
                    // cover previous a_seq
                    a_seq = ""; // reset a_seq
                } else {
                    a_seq = a_seq + a_line;
                }
            }
            
            // save last seq rec
            if (added_seq.count(a_seq) == 0) {// if no prescence before
                // added_seq.insert(a_seq); no need to insert last one

                fasta_out_put << ">" + s_id + "-" << seq_id << " " + info + "\n"; // using previour seq_id
                
                // save 60 nchar a row
                for (saver = 0; saver < a_seq.length(); saver += 60) {
                    fasta_out_put << a_seq.substr(saver, 60) + "\n";
                }
            } 

            fasta_in_put.close();
            fasta_out_put.close();
            return 0;
        });
    }

    return 0;
}