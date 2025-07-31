#define WRITER_BUF 1 // MB
#define KEEP_SHORT_ID false // MB
#define MAP_RESERVE 5000 // average of EC

#include "ThreadPool.h"
#include "Utils.hpp"

#include <iostream>
#include <filesystem>

namespace fs = std::filesystem;
using namespace fasta;


int main(int argc, char** argv) {
    // parameter parsing
    std::string strain_info_tsv_pth;
    std::string copy_info_pth;
    std::string rr_ed_dir;
    bool keyInclude = KEEP_SHORT_ID;
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
        else if (std::string(argv[i]) == "--copy_info_path" || std::string(argv[i]) == "-c")
        {
            copy_info_pth = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--thread_number" || std::string(argv[i]) == "-u")
        {
            process_num = std::stoi(argv[i + 1]);
        }
        else if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h")
        {
            std::cout << "Thanks for using OrthoSLC! (version: " << __version__ << ")\n\n";
            std::cout << "Usage: Step2_simple_derep -i input_file -o output/ [options...]\n\n";
            std::cout << "  -i or --input_path ------> <txt> path/to/file/output/by/Step1\n";
            std::cout << "  -o or --output_path -----> <dir> path/to/output/directory\n";
            std::cout << "  -c or --copy_info_path --> <txt> path/to/output/copy_info.txt\n";
            std::cout << "  -u or --thread_number ---> <int> thread number, default: 1\n";
            std::cout << "  -h or --help ------------> display this information\n";
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

    fs::path parent_copy_info_pth = fs::path(copy_info_pth).parent_path();
    if (!(fs::exists(parent_copy_info_pth))) {
        std::cerr << "Error: parent path provided to '-c or --copy_info_path' does not exist. 路径不存在\n";
        exit(0);
    }

    // if op path exit
    if (!(fs::exists(rr_ed_dir))){
        fs::create_directory(rr_ed_dir);
    }
    
    // get abs pth and short id ----
    tsv::Map pth_m;
    tsv::readTSV(strain_info_tsv_pth, pth_m);

    // mt ----
    // Thread-safe data structures for copy info
    tsv::Map copy_info;
    std::mutex copy_info_mutex;
    void (*TSVwriterFp)(const std::string&, const tsv::Map&) 
        = keyInclude 
            ? &tsv::writeTSVWithKey 
            : &tsv::writeTSV;


    ThreadPool pool(process_num);
    // 存储所有任务的 future 对象
    std::vector<std::future<void>> futures;

    for (const auto& kv : pth_m) {
        const std::string& in_p = kv.second[2];
        const std::string& save_name = kv.second[1];
        const std::string& s_id = kv.first;
        const std::string* op = &rr_ed_dir;
        futures.emplace_back(pool.enqueue([in_p, save_name, s_id, op, &copy_info, &copy_info_mutex] { // regradless of number of input parameter
            // read in file is 
            Reader reader(in_p);
            // for file to save
            // std::string file_naam = fs::path(in_p).filename();
            std::string op_full_path = fs::path(*op) / fs::path(save_name + ".fasta");
            Writer writer(op_full_path, WRITER_BUF * 1024 * 1024);

            DeduplicatorExact dedup;  // or DeduplicatorExact for perfect identity
            dedup.reserve(MAP_RESERVE);

            int seq_id = 0; // start
            std::string header, seq;
             while (reader.next(header, seq)) {
                const std::string id = header.substr(0, header.find(' '));

                bool is_new = dedup.add_and_get(seq, id).second;
                if (is_new) {
                    std::string ortho_id = s_id + "-" + std::to_string(seq_id);
                    seq_id++;
                    writer.write(ortho_id + " " + header, seq);
                }
            }
            writer.flush();

            std::unique_lock<std::mutex> lock(copy_info_mutex);
            const auto& hashMap = dedup.map();
            for (const auto& [key, vec] : hashMap) {
                if (vec.size() > 1) {
                    copy_info.emplace(vec[0], vec);
                }
            }
            lock.unlock();

            return;
        }));
    }

    // COPY INFO ----
    // 等待所有任务完成
    for (auto& future : futures) {
        future.wait();
    }
    TSVwriterFp(copy_info_pth, copy_info);

    return 0;
}