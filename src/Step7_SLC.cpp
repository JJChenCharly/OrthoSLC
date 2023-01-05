#include "ThreadPool.h"

#include <filesystem>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <sstream>

namespace fs = std::filesystem;

int main(int argc, char** argv) {
    // parameter parsing ----
    std::string in_dir;
    std::string out_dir;
    int process_num = 1;
    int compression_size = 10;
    bool one_step = false;

    for (int i = 1; i < argc; ++i)
    {
        if (std::string(argv[i]) == "--input_path" || std::string(argv[i]) == "-i")
        {
            in_dir = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--output_path" || std::string(argv[i]) == "-o")
        {
            out_dir = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--thread_number" || std::string(argv[i]) == "-u")
        {
            process_num = std::stoi(argv[i + 1]);
        }
        else if (std::string(argv[i]) == "--compression_size" || std::string(argv[i]) == "-S")
        {
            if (std::string(argv[i + 1]) == "all") {
                one_step = true;
            } else {
                compression_size = std::stoi(argv[i + 1]);
            }
        }
        else if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h")
        {
            std::cout << "Thanks for using OrthoSLC! (version: 0.1Alpha)\n\n";
            std::cout << "Usage: Step7_SLC -i input/ -o output/ [options...]\n\n";
            std::cout << "options: " << "\n";
            std::cout << "  -i or --input_path -------> path/to/input/directory\n";
            std::cout << "  -o or --output_path ------> path/to/output/directory\n";
            std::cout << "  -u or --thread_number ----> thread number, default: 1\n";
            std::cout << "  -S or --compression_size -> compression size, default: 10, 'all' means one-step-to-final\n";
            std::cout << "  -h or --help -------------> display this information\n";
            exit(0);
        }
    }

    // if op path exit ----
    if (!(fs::exists(out_dir))){
        fs::create_directory(out_dir);
    }

    //mission lst and spliter ----
    std::vector<std::string> mission_lst;
    std::string current_path;
    for (const auto& entry : fs::directory_iterator(in_dir)) {
        current_path = entry.path();
        current_path = fs::absolute(current_path);
        mission_lst.push_back(current_path);
    }

    if (one_step) {
        compression_size = mission_lst.size();
    }

    int task_len = static_cast<int>(mission_lst.size());
    // // if wrong compression_size
    if (compression_size > mission_lst.size()) {
        std::cerr << "Error: -S or --compression_size should be smaller than " << task_len << "\n";
        exit(0);
    }

    // // m_lst_
    typedef std::unordered_map<int, std::vector<std::string*>> m_lst;
    m_lst mission_lst_;

    int remainer = task_len%compression_size;

    int rem = 0;
    int f = 0;
    for (int ind = 0; ind < task_len; ind+0) {
        if (rem < remainer) {
            for (int step_len = 0; step_len < (compression_size + 1) && ind < task_len; step_len++) {
                mission_lst_[f].push_back(&mission_lst[ind]);
                ind++;
            }
            rem++;
        } else {
            for (int step_len = 0; step_len < compression_size && ind < task_len; step_len++) {
                mission_lst_[f].push_back(&mission_lst[ind]);
                ind++;
            }
        }
        f++;
         
    }
    //  mt ----
    if (mission_lst_.size() < process_num) {
        process_num = static_cast<int>(mission_lst_.size());
    }

    ThreadPool pool(process_num);

    m_lst::iterator it = mission_lst_.begin();

    while (it != mission_lst_.end()) {
        std::string op_ind = std::to_string(it->first);
        std::vector<std::string*> fls = it->second;
        std::string* op = &out_dir;

        pool.enqueue([op_ind, fls, op]{
            
            //clustering
            std::unordered_map<std::string, long int> ptr_coll;
            std::vector<long int> bin_connector;
            std::vector<std::unordered_set<std::string>> to_save;

            long int _ind_every = 0;
            std::string a_line;

            for (std::size_t i = 0; i < fls.size(); ++i) {
                
                std::ifstream a_in(*fls[i]);

                while(getline(a_in, a_line)) {
                    
                    std::unordered_set<std::string> elements;
                    std::stringstream ss(a_line);
                    std::string element;

                    while (std::getline(ss, element, '\t')) {
                        elements.insert(element);
                    }
                    to_save.push_back(elements);

                    long int ind_to_go = _ind_every;

                    bin_connector.push_back(_ind_every);

                    // traverse set just made
                    std::unordered_set<std::string>::iterator _it_ = elements.begin();

                    while (_it_ != elements.end()) {
                        if (ptr_coll.count(*_it_) > 0) { // if prescence before
                            long int ind_ = ptr_coll.at(*_it_);

                            while (bin_connector[ind_] != ind_) { // to the lowest ind
                                ind_ = bin_connector[ind_];
                            }

                            long int its_lowest_bin_ind = ind_;

                            if (its_lowest_bin_ind == ind_to_go) {
                                ++_it_;
                                continue;
                            } 

                            if (its_lowest_bin_ind > ind_to_go) {
                                
                                long int tmp = ind_to_go;
                                
                                ind_to_go = its_lowest_bin_ind;
                                its_lowest_bin_ind = tmp;
                            }  

                            to_save[its_lowest_bin_ind].insert(to_save[ind_to_go].begin(), to_save[ind_to_go].end());

                            to_save[ind_to_go].clear();

                            bin_connector[ind_to_go] = its_lowest_bin_ind;
                            
                            ind_to_go = its_lowest_bin_ind;

                        }  else { // if no prescence
                            ptr_coll[*_it_] = ind_to_go;
                        }
                        ++_it_;
                        
                    }
                    _ind_every++;
                }
                a_in.close();
            }

            // saving
            std::ofstream a_out(*op + "/" + op_ind + ".txt", std::ios::trunc);
            
            std::vector<std::unordered_set<std::string>>::iterator a_set = to_save.begin();
            while (a_set != to_save.end()) {
                if (a_set->empty()) {
                    ++a_set;
                    continue;
                } else {
                    std::unordered_set<std::string>::iterator an_id = a_set->begin();
                    size_t last_or_not = 0;
                    size_t size = a_set->size();
                    while (an_id != a_set->end()) {
                        last_or_not++; 
                        if (last_or_not == size) {
                            a_out << *an_id << "\n";
                        } else {
                            a_out << *an_id << "\t";
                        }
                        an_id++;
                    }

                }
                ++a_set;
            }

            a_out.close();

        });
        it++;
    }

    return 0;
}