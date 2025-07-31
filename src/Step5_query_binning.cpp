#include "ThreadPool.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <filesystem>

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
            std::cout << "Thanks for using OrthoSLC! (version: " << __version__ << ")\n\n";
            std::cout << "Usage: Step5_query_binning -i input/ -o output/ [options...]\n\n";
            std::cout << "  -i or --input_path -----> <dir> path/to/input/directory of blast outputs\n";
            std::cout << "  -o or --output_path ----> <dir> path/to/output/directory\n";
            std::cout << "  -u or --thread_number --> <int> thread number, default: 1\n";
            std::cout << "  -L or --bin_level ------> <int> binning level, an intger 0 < L <= 9999, default: 10\n";
            std::cout << "  -k or --no_lock_mode ---> <on/off> select to turn no lock mode <on> or <off>, default: off\n";
            std::cout << "  -h or --help -----------> display this information\n";
            exit(0);
        }
    }
    // if file path exist
    if (!(fs::exists(blast_op_pth))) {
        std::cerr << "Error: path provided to '-i or --input_path' does not exist. 路径不存在\n";
        exit(0);
    }

    fs::path parent_bin_op_pth = fs::path(bin_op_pth).parent_path();
    if (!(fs::exists(parent_bin_op_pth))) {
        std::cerr << "Error: parent path provided to '-o or --output_path' does not exist. 路径不存在\n";
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
    std::vector<std::future<void>> futures;

    for(int i = 0; i < missions.size(); ++i) {
        std::string* in_p = &missions[i];
        std::string* bin_op = &bin_op_pth;
        int b_level = bin_level;
        std::mutex* LOCK = op_mutex;
        bool lock_or_not = no_lock_mode;

        futures.emplace_back(
            pool.enqueue([in_p, bin_op, b_level, lock_or_not, LOCK] {
                std::ifstream a_blast_op(*in_p);

                std::hash<std::string> hash_fn;

                std::unordered_map<int, std::vector<std::string>> bins_to_save;
                std::string a_line;
                std::string query, subject;

                while (std::getline(a_blast_op, a_line)) {
                    std::stringstream ss_(a_line);

                    // std::string query = a_line.substr(0, a_line.find("\t"));
                    std::getline(ss_, query, '\t');
                    std::getline(ss_, subject, '\t');
                    if (query.substr(0, query.find("-")) == subject.substr(0, subject.find("-"))) continue;
                    
                    size_t h = hash_fn(query);
                    int b = h % b_level;
                    
                    bins_to_save[b].push_back(a_line + "\n");
                }
                a_blast_op.close();

                // saver ----
                bool use_lock = !lock_or_not;  // rename for clarity

                auto save_bin = [&](int bin_id, const std::vector<std::string>& items) {
                    // build the path correctly:
                    auto filename  = std::to_string(bin_id) + ".txt";
                    auto full_path = fs::path(*bin_op) / filename;

                    if (use_lock) {
                        std::unique_lock<std::mutex> lock(LOCK[bin_id]);
                        std::ofstream ofs(full_path, std::ios::app);
                        for (auto const &line : items) ofs << line;
                    } else {
                        std::ofstream ofs(full_path, std::ios::app);
                        for (auto const &line : items) ofs << line;
                    }
                };

                for (auto const &kv : bins_to_save) {
                    save_bin(kv.first, kv.second);
                }
            })
        );


    }
    
    for (auto& future : futures) {
        future.wait();
    }

    return 0;
}