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
    std::string id_info_path;
    std::string pre_cluster_path;
    std::string dereped_cated_fasta_pth;
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
            dereped_cated_fasta_pth = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--pre_cluster_path" || std::string(argv[i]) == "-p")
        {
            pre_cluster_path = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--id_info_path" || std::string(argv[i]) == "-m")
        {
            id_info_path = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--total_count" || std::string(argv[i]) == "-c")
        {
            total_count = std::stoi(argv[i + 1]);
        }
        else if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h")
        {
            std::cout << "Thanks for using OrthoSLC! (version: " << __version__ << ")\n\n";
            std::cout << "Usage: Step9_write_clusters -i input_path -o output/ -f concatenated.fasta [options...]\n\n";
            std::cout << "options:\n";
            std::cout << "  -i or --input_path --------> <dir> path/to/input/final_cluster_file from Step 8\n";
            std::cout << "  -o or --output_path -------> <dir> path/to/output/directory\n";
            std::cout << "  -f or --fasta_path --------> <fasta> path/to/dereped_cated_fasta from Step 3\n";
            std::cout << "  -m or --id_info_path ------> <txt> path/to/output/id_info_table from Step 3\n";
            std::cout << "  -p or --pre_cluster_path --> <txt> path/to/pre_clustered_file from Step 3\n";
            std::cout << "  -c or --total_count -------> <int> amonut of genomes to analyze\n";
            std::cout << "  -t or --cluster_type ------> <txt> select from < accessory / strict / surplus >, separate by ',', all types if not specified\n";
            std::cout << "  -u or --thread_number -----> <int> thread number, default: 1\n";
            std::cout << "  -h or --help --------------> display this information\n";
            exit(0);
        }
    }
    // if file path exist
    if (!(fs::exists(final_cluster_path))) {
        std::cerr << "Error: path provided to '-i or --input_path' does not exist. 路径不存在\n";
        exit(0);
    }
    else if (!(fs::exists(dereped_cated_fasta_pth))) {
        std::cerr << "Error: path provided to '-f or --fasta_path' does not exist. 路径不存在\n";
        exit(0);
    }
    else if (!(fs::exists(pre_cluster_path))) {
        std::cerr << "Error: path provided to '-p or --pre_cluster_path' does not exist. 路径不存在\n";
        exit(0);
    }
    else if (!(fs::exists(id_info_path))) {
        std::cerr << "Error: path provided to '-m or --id_info_path' does not exist. 路径不存在\n";
        exit(0);
    }

    fs::path parent_output_path = fs::path(output_path).parent_path();
    if (!(fs::exists(parent_output_path))) {
        std::cerr << "Error: parent path provided to '-o or --output_path' does not exist. 路径不存在\n";
        exit(0);
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
            if (!(fs::exists(dir_path))) {
                fs::create_directory(dir_path);
            }
        } else if (str == "strict") {
            types_to_write.at("strict_core") = true;
            std::filesystem::path dir_path = fs::path(output_path) / fs::path("strict_core");
            if (!(fs::exists(dir_path))) {
                fs::create_directory(dir_path);
            }
        } else if (str == "surplus") {
            types_to_write.at("surplus_core") = true;
            std::filesystem::path dir_path = fs::path(output_path) / fs::path("surplus_core");
            if (!(fs::exists(dir_path))) {
                fs::create_directory(dir_path);
            }
        }
    }
    
// mission list and final cluster ----
    std::ifstream final_CLUSTER_file(final_cluster_path);
    std::vector<std::vector<std::string>> CLUSTERS;
    // [[id1,id2], [id3, id4]]
    
    std::string a_line;

    while (getline(final_CLUSTER_file, a_line)) {
        std::vector<std::string> v;
        std::stringstream iss(a_line);
        std::string element;

        while (std::getline(iss, element, '\t')) {
            v.push_back(element);
        }
    
        CLUSTERS.push_back(v);
    }
    final_CLUSTER_file.close();

// id pre cluster dcit ----
    std::ifstream pre_CLUSTER_file(pre_cluster_path);

    typedef std::unordered_map<std::string, std::string> s_s;
    s_s id_cluster_dict;
    // {id1: id1,
    // id2: id1}

    while (getline(pre_CLUSTER_file, a_line)) {
        std::vector<std::string> v;
        std::stringstream iss(a_line);
        std::string element;

        while (std::getline(iss, element, '\t')) {
            v.push_back(element);

            id_cluster_dict[element] = v[0];
        }
    }
    pre_CLUSTER_file.close();

// id info dict ----
    std::ifstream ID_INFO_file(id_info_path);
    typedef std::unordered_map<std::string, std::string> s_s;
    s_s id_info_dict;

    while (getline(ID_INFO_file, a_line)) {
        std::vector<std::string> v;
        std::stringstream iss(a_line);
        std::string element;

        while (std::getline(iss, element, '\t')) {
            v.push_back(element);
        }

        id_info_dict[v[0]] = v[1];
    }
    ID_INFO_file.close();
    // {id: description}

// dereped cat fasta ----
    s_s fasta_dict;

    std::ifstream dereped_cat_file(dereped_cated_fasta_pth);
    
    // first one
    getline(dereped_cat_file, a_line);
    std::string key = a_line.substr(1, 11);
    std::string seq = "";

    while (getline(dereped_cat_file, a_line)) {
        if (a_line[0] == '>') {
            // save previous one
            fasta_dict[key] = seq;

            key = a_line.substr(1, 11);
            seq = "";
        } else {
            seq = seq + a_line + "\n";
        }
    }
    // last one
    fasta_dict[key] = seq;

    dereped_cat_file.close();

// mt ----
    ThreadPool pool(process_num);

    for (const std::vector<std::string>& a_c : CLUSTERS) {
        std::vector<std::string> the_cluster = a_c;
        int* total_amount = &total_count;
        std::unordered_map<std::string, bool>* which_to_write = &types_to_write;
        std::string op_p = output_path;
        s_s* ptr_fasta_dict = &fasta_dict;
        s_s* ptr_id_info_dict = &id_info_dict;
        s_s* ptr_id_cluster_dict = &id_cluster_dict;

        pool.enqueue([the_cluster, total_amount, which_to_write, op_p, ptr_fasta_dict, ptr_id_info_dict, ptr_id_cluster_dict] {
            
            std::unordered_set<std::string> spe_count;

            for (const std::string& s : the_cluster) {
                spe_count.insert(s.substr(0, 5));
            }

            fs::path file_naam;
            bool ending = false;

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
                    fasta_out_put << ">" + a_g + " " + ptr_id_info_dict->at(a_g) + "\n" + ptr_fasta_dict->at(ptr_id_cluster_dict->at(a_g)); 
                }
                
                fasta_out_put.close();
            } 
        });   
    }

    return 0;
}