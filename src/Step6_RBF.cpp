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
    int bin_level = 2;

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
            if(bin_level > 4) {
                std::cerr << "-L or --bin_level must be less than 4";
                exit(0);
            }
        }
        else if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h")
        {
            std::cout << "Thanks for using OrthoSLC! (version: 0.1Alpha)\n\n";
            std::cout << "Usage: Step6_RBF -i input/ -o output/ [options...]\n\n";
            std::cout << "  -i or --input_path -------> path/to/input/directory\n";
            std::cout << "  -o or --output_path ------> path/to/output/directory\n";
            std::cout << "  -u or --thread_number ----> thread number, default: 1\n";
            std::cout << "  -L or --bin_level --------> binning level, an intger 0 < L <= 4 , default: 2, means 100 bins\n";
            std::cout << "  -h or --help -------------> display this information\n";
            exit(0);
        }
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

    double bin_nums = pow(10.0, bin_level);
    std::mutex op_mutex[static_cast<int>(bin_nums)];

    for(int i = 0; i < mission_lst.size(); ++i) {
        std::string* op = &RBB_op_dir_pth;
        std::string* in_p = &mission_lst[i];
        std::mutex* LOCK = op_mutex;
        int b_level = bin_level;
        
        pool.enqueue([op, in_p, LOCK, b_level]{

            std::ifstream a_bin_op(*in_p);
            std::unordered_set<std::string> rec;
            std::unordered_map<std::string, std::vector<std::string>> bins_to_save;

            std::hash<std::string> hash_fn;

            std::string a_line;
            while(getline(a_bin_op, a_line)) {
                if (rec.count(a_line) != 0) { // if prescence 
                    std::string first;
                    std::getline(std::stringstream(a_line), first, '\t');

                    size_t h = hash_fn(first);
                    std::stringstream h_to_str;
                    h_to_str << h;
                    std::string h_str = h_to_str.str();

                    std::string its_bin = h_str.substr(h_str.length() - b_level);

                    bins_to_save[its_bin].push_back(a_line);
                    
                    rec.erase(a_line);
                } else {
                    rec.insert(a_line);
                }
            }

            a_bin_op.close();

            // saving ----
            for (auto it3 = bins_to_save.begin(); it3 != bins_to_save.end(); ++it3) {
                const std::string& the_bin = it3->first;
                const std::vector<std::string> vec_ = it3->second;

                std::unique_lock<std::mutex> lock(LOCK[std::stoi(the_bin)]);
                std::ofstream RBB_out(*op + "/" + the_bin + ".txt", std::ios::app);

                for (std::vector<std::string>::const_iterator it4 = vec_.begin(); 
                    it4 != vec_.end(); 
                    ++it4) {
                    RBB_out << *it4 <<"\n";
                }
                RBB_out.close();
                lock.unlock();
            }            
        });
    }

    return 0;
}