#include "ThreadPool.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <filesystem>
#include <utility>
#include <cmath>
#include <functional>
#include <iomanip>
#include <mutex>

namespace fs = std::filesystem;

int main(int argc, char** argv) {
    // parameter parsing
    std::string blast_op_pth;
    std::string bin_op_pth;
    std::string seq_len_pth;
    int bin_level = 2;
    double r_ = 0.3;
    int process_num = 1;

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
        else if (std::string(argv[i]) == "--seq_len_path" || std::string(argv[i]) == "-s")
        {
            seq_len_pth = argv[i + 1];
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
        else if (std::string(argv[i]) == "--length_limit" || std::string(argv[i]) == "-r")
        {
            r_ = std::stod(argv[i + 1]);
        }
        else if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h")
        {
            std::cout << "Thanks for using OrthoSLC! (version: 0.1Alpha)\n\n";
            std::cout << "Usage: Step5_filter_n_bin -i input/ -o output/ -s seq_len_info.txt [options...]\n\n";
            std::cout << "  -i or --input_path -------> path/to/input/directory\n";
            std::cout << "  -o or --output_path ------> path/to/output/directory\n";
            std::cout << "  -s or --seq_len_path -----> path/to/output/seq_len_info.txt\n";
            std::cout << "  -u or --thread_number ----> thread number, default: 1\n";
            std::cout << "  -L or --bin_level --------> binning level, an intger 0 < L <= 4 , default: 2, means 100 bins\n";
            std::cout << "  -r or --length_limit -----> length difference limit, default: 0.3\n";
            std::cout << "  -h or --help -------------> display this information\n";
            exit(0);
        }
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

    // mission VEC ----
    std::vector<std::string> missions;
    std::string current_path;
    for (const auto& entry : fs::directory_iterator(blast_op_pth)) {
        current_path = entry.path();
        current_path = fs::absolute(current_path);
        missions.push_back(current_path);
    }

    // mt ----
    // mutex each op file
    double bin_nums = pow(10.0, bin_level);
    std::mutex op_mutex[static_cast<int>(bin_nums)];

    ThreadPool pool(process_num);

    for(int i = 0; i < missions.size(); ++i) {
        std::string* in_p = &missions[i];
        std::string bin_op = bin_op_pth;
        std::unordered_map<std::string, double>* seq_len_ = &seq_len;
        int b_level = bin_level;
        double r_limit = r_;
        std::mutex* LOCK = op_mutex;

        pool.enqueue([in_p, bin_op, seq_len_, b_level, r_limit, LOCK] { // regradless of number of input parameter
                
            std::ifstream a_blast_op(*in_p);

            typedef std::unordered_map<std::string, std::pair<std::string, double>> how_to_rec;
            how_to_rec rec;

            std::string a_hit;
            std::string query, subject, score;

            // filter ----
            while (std::getline(a_blast_op, a_hit)) {
                std::stringstream ss_(a_hit);
    
                std::getline(ss_, query, '\t');
                std::getline(ss_, subject, '\t');
                std::getline(ss_, score);

                if (query.substr(0, 5) == subject.substr(0, 5)) { // if paralog
                    continue;
                } else {
                    double q_len = seq_len_->at(query);
                    double s_len = seq_len_->at(subject);
                    double len_ratio = q_len/s_len;
                    
                    bool c1 = (len_ratio < r_limit);
                    bool c2 = (len_ratio > (1/r_limit));

                    if (c1 || c2) { // if len too diff
                        continue;
                    } else {
                        double score_ = std::stod(score);
                        // if already hit
                        how_to_rec::iterator it = rec.find(query);
                        if (it != rec.end()) { // if already hit
                            if (score_ > rec[query].second) { // if a better hit
                                std::pair<std::string, double> p = std::make_pair(subject, score_);
                                rec[query] = p;

                            } else if (score_ == rec[query].second) { // if same score

                                double query_len = seq_len_->at(query);
                                double len1 = seq_len_->at(subject) - query_len;// double len2 = seq_len_.at(rec[query].first) - query_len; // len of existing sub of the query
                                double len2 = seq_len_->at(rec[query].first) - query_len;
                                if (fabs(len1) < fabs(len2)) { // if more similar
                                    std::pair<std::string, double> p = std::make_pair(subject, score_);
                                    rec[query] = p;

                                }
                            }
                        } else { // if no hit bofore
                            std::pair<std::string, double> p = std::make_pair(subject, score_);
                            rec[query] = p;
                        }
                    }
                }
            }
            a_blast_op.close();

            // binning ----
            std::unordered_map<std::string, std::vector<std::string>> bins_to_save;
            std::hash<std::string> hash_fn;

            for (auto it2 = rec.begin(); it2 != rec.end(); ++it2) {
                const std::string& qu = it2->first;
                const std::pair<std::string, double>& value2 = it2->second;
                const std::string& su = value2.first;

                std::string p_to_s;
                if (qu < su) {
                    p_to_s = qu + "\t" + su;
                } else {
                    p_to_s = su + "\t" + qu;
                }
                
                size_t h = hash_fn(p_to_s);
                std::stringstream h_to_str;
                h_to_str << h;
                std::string h_str = h_to_str.str();

                std::string its_bin = h_str.substr(h_str.length() - b_level);

                bins_to_save[its_bin].push_back(p_to_s);

            }
            
            // saving ----
            for (auto it3 = bins_to_save.begin(); it3 != bins_to_save.end(); ++it3) {
                const std::string& the_bin = it3->first;
                const std::vector<std::string> vec_ = it3->second;

                std::unique_lock<std::mutex> lock(LOCK[std::stoi(the_bin)]);
                std::ofstream save_into_a_bin(bin_op + '/' + the_bin + ".txt", std::ios::app);

                for (std::vector<std::string>::const_iterator it4 = vec_.begin(); 
                    it4 != vec_.end(); 
                    ++it4) {
                    save_into_a_bin << *it4 <<"\n";
                }
                save_into_a_bin.close();
                lock.unlock();
            }
        });
    }

    return 0;
}