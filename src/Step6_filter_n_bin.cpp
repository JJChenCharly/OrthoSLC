#include "ThreadPool.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <filesystem>
#include <utility>
#include <cmath>
#include <functional>
#include <iomanip>
#include <mutex>

namespace fs = std::filesystem;

int main(int argc, char** argv) {
    // parameter parsing
    std::string q_bin_op;
    std::string bin_op_pth;
    std::string seq_len_pth;
    std::string pre_cluster_path;
    int bin_level = 10;
    double r_ = 0.3;
    int process_num = 1;
    bool no_lock_mode = false;

    for (int i = 1; i < argc; ++i)
    {
        if (std::string(argv[i]) == "--input_path" || std::string(argv[i]) == "-i")
        {
            q_bin_op = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--output_path" || std::string(argv[i]) == "-o")
        {
            bin_op_pth = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--seq_len_path" || std::string(argv[i]) == "-s")
        {
            seq_len_pth = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--pre_cluster_path" || std::string(argv[i]) == "-p")
        {
            pre_cluster_path = argv[i + 1];
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
            if((bin_level < 0) || (bin_level > 9999)) {
                std::cerr << "-L or --bin_level must be integer 0 < L <= 9999";
                exit(0);
            }
        }
        else if (std::string(argv[i]) == "--length_limit" || std::string(argv[i]) == "-r")
        {
            r_ = std::stod(argv[i + 1]);
            if((r_ < 0) || (r_ > 1)) {
                std::cerr << "-r or --length_limit must be integer 0 < r <= 1";
                exit(0);
            }
        }
        else if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h")
        {
            std::cout << "Thanks for using OrthoSLC! (version: " << __version__ << ")\n\n";
            std::cout << "Usage: Step6_filter_n_bin -i input/ -o output/ -s seq_len_info.txt [options...]\n\n";
            std::cout << "  -i or --input_path --------> <dir> path/to/input/directory from Step 5\n";
            std::cout << "  -o or --output_path -------> <dir> path/to/output/directory\n";
            std::cout << "  -s or --seq_len_path ------> <txt> path/to/seq_len_info.txt from Step 3\n";
            std::cout << "  -p or --pre_cluster_path --> <txt> path/to/pre_cluster.txt from Step 3\n";
            std::cout << "  -L or --bin_level ---------> <int> binning level, an intger 0 < L <= 9999 , default: 10\n";
            std::cout << "  -r or --length_limit ------> <float> length difference limit, 0 < r <= 1, default: 0.3\n";
            std::cout << "  -k or --no_lock_mode ------> <on/off> select to turn no lock mode <on> or <off>, default: off\n";
            std::cout << "  -u or --thread_number -----> <int> thread number, default: 1\n";
            std::cout << "  -h or --help --------------> display this information\n";
            exit(0);
        }
    }
    // if file path exist
    if (!(fs::exists(q_bin_op))) {
        std::cerr << "Error: path provided to '-i or --input_path' does not exist. 路径不存在\n";
        exit(0);
    } else if (!(fs::exists(seq_len_pth))) {
        std::cerr << "Error: path provided to '-s or --seq_len_path' does not exist. 路径不存在\n";
        exit(0);
    } else if (!(fs::exists(pre_cluster_path))) {
        std::cerr << "Error: path provided to '-p or --pre_cluster_path' does not exist. 路径不存在\n";
        exit(0);
    }

    fs::path parent_bin_op_pth = fs::path(bin_op_pth).parent_path();
    if (!(fs::exists(parent_bin_op_pth))) {
        std::cerr << "Error: parent path provided to '-o or --output_path' do not exist. 路径不存在\n";
        exit(0);
    }

    // if op path exit
    if (!(fs::exists(bin_op_pth))){
        fs::create_directory(bin_op_pth);
    }

    // read seq len ----
    std::unordered_map<std::string, double> seq_len;

    std::ifstream seq_len_file(seq_len_pth);

    std::string a_line;
    std::string key, value;
    while (std::getline(seq_len_file, a_line)) {
        // Split the line into two parts using a stringstream and getline
        std::stringstream ss(a_line);
        
        std::getline(ss, key, '\t');
        std::getline(ss, value);

        // Insert the key-value pair into the map
        seq_len[key] = std::stod(value);
    }

    seq_len_file.close();

    // read pre_cluster ----
    std::ifstream pre_cluster(pre_cluster_path);
    typedef std::unordered_map<std::string, std::unordered_set<std::string>> um_s_uss;
    um_s_uss bin_spe;

    while(getline(pre_cluster, a_line)) {
        std::stringstream ss(a_line);
        std::string element;

        int key_ = 0;
        std::string key;
        while (std::getline(ss, element, '\t')) {
            if (key_ == 0) {
                key = element;
            }
            key_++;
            bin_spe[key].insert(element.substr(0, 5));
        }
        
    }

    pre_cluster.close();

    // mission VEC ----
    std::vector<std::string> missions;
    std::string current_path;
    for (const auto& entry : fs::directory_iterator(q_bin_op)) {
        current_path = entry.path();
        current_path = fs::absolute(current_path);
        missions.push_back(current_path);
    }

    // mt ----
    // mutex each op file
    std::mutex op_mutex[bin_level];

    ThreadPool pool(process_num);

    for(int i = 0; i < missions.size(); ++i) {
        std::string* in_p = &missions[i];
        std::string* bin_op = &bin_op_pth;
        std::unordered_map<std::string, double>* seq_len_ = &seq_len;
        int b_level = bin_level;
        double r_limit = r_;
        um_s_uss* bin_spe_ = &bin_spe;
        std::mutex* LOCK = op_mutex;
        bool lock_or_not = no_lock_mode;

        pool.enqueue([in_p, bin_op, seq_len_, b_level, r_limit, bin_spe_, lock_or_not, LOCK] { // regradless of number of input parameter
                
            std::ifstream a_qb_op(*in_p);

            std::unordered_map<std::string, std::unordered_map<std::string, int>> q_s_d;
            std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::string>>> sspe_m;
            
            std::string a_hit;
            std::string query, subject, score;

            // data formatting ----

            while (std::getline(a_qb_op, a_hit)) {
                std::stringstream ss_(a_hit);
    
                std::getline(ss_, query, '\t');
                std::getline(ss_, subject, '\t');
                std::getline(ss_, score);

                if (query.substr(0, 5) == subject.substr(0, 5)) {
                    continue;
                } else {
                    int score_ = stoi(score);
                    q_s_d[query][subject] = score_;

                    for (const std::string& a_spe_of_clus : bin_spe_->at(subject)) {
                        sspe_m[query][a_spe_of_clus].push_back(subject);
                    }
                }

            }

            a_qb_op.close();

            // filter ----
            std::unordered_map<int, std::vector<std::string>> bins_to_save;
            std::hash<std::string> hash_fn;

            for (const auto& [qs, nested_um] : q_s_d) {
                double q_len = seq_len_->at(qs);

                for (const auto& [sub, its_score] : nested_um) {

                    std::string s_spe = sub.substr(0, 5);

                    if (bin_spe_->at(qs).find(s_spe) != bin_spe_->at(qs).end()) { // if qs pre cluster has that spe
                        continue;
                    } else {

                        bool best_hit = true;

                        std::vector<std::string> same_score;                            

                        if (sspe_m.at(qs).at(s_spe).size() > 1) { // if multiple hit
                            for (const auto &id_ : sspe_m.at(qs).at(s_spe)) {
                                if (sub == id_) {
                                    continue; // skip self 
                                }

                                if (its_score < q_s_d.at(qs).at(id_)) { // if a better hit
                                    best_hit = false;
                                    break;
                                } else if (its_score == q_s_d.at(qs).at(id_)) {
                                    same_score.push_back(id_);
                                }
                            }
                        }

                        if (same_score.size() > 0 && best_hit) {
                            double qs_len_diff = std::fabs(q_len - seq_len_->at(sub));

                            for (const auto &id_ : same_score) {
                                double q_id_len_diff = std::fabs(q_len - seq_len_->at(id_));

                                if (q_id_len_diff < qs_len_diff) {
                                    best_hit = false; // there is a seq with more similar length
                                    break;
                                }
                            }
                        }
                        
                        if (best_hit) {
                            double len_ratio = q_len/(seq_len_->at(sub));

                            bool c1 = (len_ratio <= r_limit);
                            bool c2 = (len_ratio >= (1/r_limit));

                            if (c1 || c2) { // if len too diff
                                continue;
                            } else {
                                std::string p_to_s;
                                if (qs < sub) {
                                    p_to_s = qs + "\t" + sub + "\n";
                                } else {
                                    p_to_s = sub + "\t" + qs + "\n";
                                }

                                size_t h = hash_fn(p_to_s);
                                int b = h % b_level;

                                bins_to_save[b].push_back(p_to_s);
                            }
                        }
                    }
                }
            }
            
            // saving ----
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