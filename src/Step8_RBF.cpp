#include "ThreadPool.h"

#include <iostream>
#include <fstream>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <filesystem>
#include <vector>
#include <sstream>
#include <cmath>

namespace fs = std::filesystem;

int main(int argc, char** argv) {
    // parameter parsing ----
    std::string bin_op_dir_pth;
    std::string RBB_op_dir_pth;
    int process_num = 1;
    int bin_level = 10;
    bool no_lock_mode = false;

    for (int i = 1; i < argc; ++i)
    {
        if (std::string(argv[i]) == "--input_path" || std::string(argv[i]) == "-i")
        {
            bin_op_dir_pth = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--output_path" || std::string(argv[i]) == "-o")
        {
            RBB_op_dir_pth = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--thread_number" || std::string(argv[i]) == "-u")
        {
            process_num = std::stoi(argv[i + 1]);
        }
        else if (std::string(argv[i]) == "--bin_level" || std::string(argv[i]) == "-L")
        {
            bin_level = std::stoi(argv[i + 1]);
            if((bin_level < 0) | (bin_level > 9999)) {
                std::cerr << "-L or --bin_level must be integer 0 < L <= 9999";
                exit(0);
            }
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
        else if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h")
        {
            std::cout << "Thanks for using OrthoSLC! (version: 0.1)\n\n";
            std::cout << "Usage: Step8_RBF -i input/ -o output/ [options...]\n\n";
            std::cout << "  -i or --input_path -------> path/to/input/directory\n";
            std::cout << "  -o or --output_path ------> path/to/output/directory\n";
            std::cout << "  -u or --thread_number ----> thread number, default: 1\n";
            std::cout << "  -k or --no_lock_mode -----> select to turn no lock mode <on> or <off>, default: off\n";
            std::cout << "  -L or --bin_level --------> binning level, an intger 0 < L <= 9999 , default: 10\n";
            std::cout << "  -h or --help -------------> display this information\n";
            exit(0);
        }
    }
    // if file path exist
    if (!(fs::exists(bin_op_dir_pth))) {
        std::cerr << "Error: path provided to '-i or --input_path' do not exist. 路径不存在\n";
        exit(0);
    }

    // if op path exit ----
    if (!(fs::exists(RBB_op_dir_pth))){
        fs::create_directory(RBB_op_dir_pth);
    }

    // mission list ----
    std::vector<std::string> mission_lst;
    std::string current_path;
    for (const auto& entry : fs::directory_iterator(bin_op_dir_pth)) {
        current_path = entry.path();
        current_path = fs::absolute(current_path);
        mission_lst.push_back(current_path);
    }
    
    // mt ----
    ThreadPool pool(process_num);

    std::mutex op_mutex[bin_level];

    for(int i = 0; i < mission_lst.size(); ++i) {
        std::string* op = &RBB_op_dir_pth;
        std::string* in_p = &mission_lst[i];
        std::mutex* LOCK = op_mutex;
        int b_level = bin_level;
        bool lock_or_not = no_lock_mode;
        
        pool.enqueue([op, in_p, LOCK, b_level, lock_or_not]{

            std::ifstream a_bin_op(*in_p);
            std::unordered_set<std::string> rec;
            std::unordered_map<int, std::vector<std::string>> bins_to_save;

            std::hash<std::string> hash_fn;

            std::string a_line;
            while(getline(a_bin_op, a_line)) {
                if (rec.count(a_line) != 0) { // if prescence 
                    std::string first;
                    std::getline(std::stringstream(a_line), first, '\t');

                    size_t h = hash_fn(first);
                    int b = h % b_level;

                    bins_to_save[b].push_back(a_line + "\n");
                    
                    rec.erase(a_line);
                } else {
                    rec.insert(a_line);
                }
            }

            a_bin_op.close();

            // saving ----
            // saving ----
            for (auto it = bins_to_save.begin(); it != bins_to_save.end(); ++it) {

                if (lock_or_not) {
                    std::ofstream save_into_a_bin(fs::path(*op) / fs::path(std::to_string(it->first) + ".txt"), std::ios::app);

                    for (auto const &item : it->second) {
                        save_into_a_bin << item;
                    }

                    save_into_a_bin.close();
                    
                } else {
                    std::unique_lock<std::mutex> lock(LOCK[it->first]);
                    
                    std::ofstream save_into_a_bin(fs::path(*op) / fs::path(std::to_string(it->first) + ".txt"), std::ios::app);

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