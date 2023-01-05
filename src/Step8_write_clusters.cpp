#include "ThreadPool.h"

#include <filesystem>
#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <utility>

namespace fs = std::filesystem;

int main(int argc, char** argv) {
    // parameter parsing ----
    std::string final_cluster_path;
    std::string output_path;
    std::string cated_fasta_path;
    int total_count;
    std::string cluster_type = "accessory,strict,surplus";
    int process_num = 1;

    for (int i = 1; i < argc; ++i)
    {
        if (std::string(argv[i]) == "--input_path" || std::string(argv[i]) == "-i")
        {
            final_cluster_path = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--output_path" || std::string(argv[i]) == "-o")
        {
            output_path = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--thread_number" || std::string(argv[i]) == "-u")
        {
            process_num = std::stoi(argv[i + 1]);
        }
        else if (std::string(argv[i]) == "--cluster_type" || std::string(argv[i]) == "-t")
        {
            cluster_type = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--fasta_path" || std::string(argv[i]) == "-f")
        {
            cated_fasta_path = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--total_count" || std::string(argv[i]) == "-c")
        {
            total_count = std::stoi(argv[i + 1]);
        }
        else if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h")
        {
            std::cout << "Thanks for using OrthoSLC! (version: 0.1Alpha)\n\n";
            std::cout << "Usage: Step8_write_clusters -i input_path -o output/ -f concatenated.fasta [options...]\n\n";
            std::cout << "options:\n";
            std::cout << "  -i or --input_path ------> path/to/input/cluster_file\n";
            std::cout << "  -o or --output_path -----> path/to/output/directory\n";
            std::cout << "  -f or --fasta_path ------> path to concatenated FASTA file\n";
            std::cout << "  -c or --total_count -----> amonut of genomes to analyze\n";
            std::cout << "  -t or --cluster_type ----> select from < accessory / strict / surplus >, separate by ',', all types if not specified\n";
            std::cout << "  -u or --thread_number ---> thread number, default: 1\n";
            std::cout << "  -h or --help ------------> display this information\n";
            exit(0);
        }
    }

    // if op path exit ----

    std::unordered_map<std::string, bool> types_to_write{
        {"accessory_cluster", false},
        {"strict_core", false},
        {"surplus_core", false}
    };

    if (!(fs::exists(output_path))) {
        fs::create_directory(output_path);
    }

    // types to write ----
    std::vector<std::string> to_write;

    std::stringstream ss(cluster_type);
    std::string token;
    while (std::getline(ss, token, ',')) {
        to_write.push_back(token);
    }

    for (const auto& str : to_write) {
        if (str == "accessory") {
            types_to_write.at("accessory_cluster") = true;
            std::filesystem::path dir_path = fs::path(output_path) / fs::path("accessory_cluster");
            fs::create_directory(dir_path);
        } else if (str == "strict") {
            types_to_write.at("strict_core") = true;
            std::filesystem::path dir_path = fs::path(output_path) / fs::path("strict_core");
            fs::create_directory(dir_path);
        } else if (str == "surplus") {
            types_to_write.at("surplus_core") = true;
            std::filesystem::path dir_path = fs::path(output_path) / fs::path("surplus_core");
            fs::create_directory(dir_path);
        }
    }
    
    // mission list ----
    std::ifstream CLUSTER_file(final_cluster_path);
    std::vector<std::vector<std::string>> CLUSTERS;
    std::string a_line;

    while (getline(CLUSTER_file, a_line)) {
        std::vector<std::string> v;
        std::stringstream iss(a_line);
        std::string element;

        while (std::getline(iss, element, '\t')) {
            v.push_back(element);
        }
    
        CLUSTERS.push_back(v);
    }
    
    CLUSTER_file.close();

    // cat fasta ----
    typedef std::unordered_map<std::string, std::string> FASTA_DICT;
    FASTA_DICT fasta_dict;

    std::ifstream cat_file(cated_fasta_path);
    
    // first one
    getline(cat_file, a_line);
    std::string key = a_line.substr(1, 11);
    std::string info = a_line + "\n";

    while (getline(cat_file, a_line)) {
        if (a_line[0] == '>') {
            // save previous one
            fasta_dict[key] = info;

            key = a_line.substr(1, 11);
            info = a_line + "\n";
        } else {
            info = info + a_line + "\n";
        }
    }
    // last one
    fasta_dict[key] = info;

    cat_file.close();

    // mt
    ThreadPool pool(process_num);

    for (const std::vector<std::string>& a_c : CLUSTERS) {
        std::vector<std::string> the_cluster = a_c;
        int* total_amount = &total_count;
        std::unordered_map<std::string, bool>* which_to_write = &types_to_write;
        std::string op_p = output_path;
        FASTA_DICT* ptr_fasta_dict = &fasta_dict;

        pool.enqueue([the_cluster, total_amount, which_to_write, op_p, ptr_fasta_dict] {
            std::unordered_set<std::string> spe_count;

            for (const std::string& s : the_cluster) {
                spe_count.insert(s.substr(5));
            }

            fs::path file_naam;
            bool ending;

            if (spe_count.size() < *total_amount) {
                if(which_to_write->at("accessory_cluster")) {
                    file_naam = fs::path(op_p) / fs::path("accessory_cluster") / fs::path(the_cluster[0] + ".fasta");
                } else{
                    ending = true;
                }
            } else {
                if (the_cluster.size() == *total_amount) {
                    if (which_to_write->at("strict_core")) {
                        file_naam = fs::path(op_p) / fs::path("strict_core") / fs::path(the_cluster[0] + ".fasta");
                    } else {
                        ending = true;
                    }
                } else {
                    if (which_to_write->at("strict_core")) {
                        file_naam = fs::path(op_p) / fs::path("surplus_core") / fs::path(the_cluster[0] + ".fasta");
                    } else {
                        ending = true;
                    }
                }
            }
            
            if (ending) {
                ;
            } else {
                std::ofstream fasta_out_put(file_naam, std::ios::trunc);

                for (const std::string& a_g : the_cluster) {
                    fasta_out_put << ptr_fasta_dict->at(a_g);
                }
                
                fasta_out_put.close();
            }
            
        });
        
    }


    return 0;
}