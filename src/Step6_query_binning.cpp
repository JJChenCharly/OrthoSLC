#include "ThreadPool.h"

#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <filesystem>
#include <mutex>

namespace fs = std::filesystem;

int main(int argc, char** argv) {
    // parameter parsing
    std::string blast_op_pth;
    std::string bin_op_pth;
    int bin_level = 10;
    int process_num = 1;
    bool no_lock_mode = false;

    for (int i = 1; i < argc; ++i)
    {
        if (std::string(argv[i]) == "--input_path" || std::string(argv[i]) == "-i")
        {
            blast_op_pth = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--output_path" || std::string(argv[i]) == "-o")
        {
            bin_op_pth = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--thread_number" || std::string(argv[i]) == "-u")
        {
            process_num = std::stoi(argv[i + 1]);
        }
        else if (std::string(argv[i]) == "--no_lock_mode" || std::string(argv[i]) == "-k")
        {
            std::string no_lock_mode_arg = argv[i + 1];

            if (no_lock_mode_arg == "on") {
                no_lock_mode = true;
            } else if (no_lock_mode_arg == "off") {
                no_lock_mode = false;
            } else {
                std::cerr << "Error: no lock mode option can only be <on> or <off>\n";
                exit(0);
            }
        }
        else if (std::string(argv[i]) == "--bin_level" || std::string(argv[i]) == "-L")
        {
            bin_level = std::stoi(argv[i + 1]);
            if((bin_level < 0) | (bin_level > 9999)) {
                std::cerr << "Error: -L or --bin_level must be integer 0 < L <= 9999\n";
                exit(0);
            }
        }
        else if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h")
        {
            std::cout << "Thanks for using OrthoSLC! (version: 0.1.1)\n\n";
            std::cout << "Usage: Step6_query_binning -i input/ -o output/ [options...]\n\n";
            std::cout << "  -i or --input_path -------> path/to/input/directory\n";
            std::cout << "  -o or --output_path ------> path/to/output/directory\n";
            std::cout << "  -u or --thread_number ----> thread number, default: 1\n";
            std::cout << "  -L or --bin_level --------> binning level, an intger 0 < L <= 9999, default: 10\n";
            std::cout << "  -k or --no_lock_mode -----> select to turn no lock mode <on> or <off>, default: off\n";
            std::cout << "  -h or --help -------------> display this information\n";
            exit(0);
        }
    }
    // if file path exist
    if (!(fs::exists(blast_op_pth))) {
        std::cerr << "Error: path provided to '-i or --input_path' do not exist. 路径不存在\n";
        exit(0);
    }

    fs::path parent_bin_op_pth = fs::path(bin_op_pth).parent_path();
    if (!(fs::exists(parent_bin_op_pth))) {
        std::cerr << "Error: parent path provided to '-o or --output_path' do not exist. 路径不存在\n";
        exit(0);
    }

    // mission lst ----
    std::vector<std::string> missions;
    std::string current_path;
    for (const auto& entry : fs::directory_iterator(blast_op_pth)) {
        current_path = entry.path();
        current_path = fs::absolute(current_path);
        missions.push_back(current_path);
    }

    // if op path exit
    if (!(fs::exists(bin_op_pth))){
        fs::create_directory(bin_op_pth);
    }

    // mt ----
    std::mutex op_mutex[bin_level];

    ThreadPool pool(process_num);

    for(int i = 0; i < missions.size(); ++i) {
        std::string* in_p = &missions[i];
        std::string* bin_op = &bin_op_pth;
        int b_level = bin_level;
        std::mutex* LOCK = op_mutex;
        bool lock_or_not = no_lock_mode;

        pool.enqueue([in_p, bin_op, b_level, lock_or_not, LOCK] {
            std::ifstream a_blast_op(*in_p);

            std::hash<std::string> hash_fn;

            std::unordered_map<int, std::vector<std::string>> bins_to_save;
            std::string a_line;

            while (std::getline(a_blast_op, a_line)) {
                a_line = a_line + "\n";

                std::string query = a_line.substr(0, 11);
                
                size_t h = hash_fn(query);
                int b = h % b_level;
                
                bins_to_save[b].push_back(a_line);
            }
            a_blast_op.close();

            // saver ----
            for (auto it = bins_to_save.begin(); it != bins_to_save.end(); ++it) {

                if (lock_or_not) {
                    std::ofstream save_into_a_bin(fs::path(*bin_op) / fs::path(std::to_string(it->first) + ".txt"), std::ios::app);

                    for (auto const &item : it->second) {
                        save_into_a_bin << item;
                    }

                    save_into_a_bin.close();
                    
                } else {
                    std::unique_lock<std::mutex> lock(LOCK[it->first]);
                    
                    std::ofstream save_into_a_bin(fs::path(*bin_op) / fs::path(std::to_string(it->first) + ".txt"), std::ios::app);

                    for (auto const &item : it->second) {
                        save_into_a_bin << item;
                    }

                    save_into_a_bin.close();
                    lock.unlock();
                }
            }

        });
    }

    return 0;
}