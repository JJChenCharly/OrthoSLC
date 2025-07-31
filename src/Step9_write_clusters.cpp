#include "ThreadPool.h"
#include "Utils.hpp"

#include <filesystem>
#include <unordered_set>
#include <utility>
#include <cmath>

namespace fs = std::filesystem;
using namespace fasta;

constexpr std::size_t WRITER_BUF = 1 * 1024 * 1024; // 2 MB buffer for fasta::Writer

int main(int argc, char** argv) {
    // parameter parsing ----
    std::string final_cluster_path;
    std::string output_path;
    std::string id_info_path;
    std::string pre_cluster_path;
    std::string derep_fasta_path;
    int total_count = 0;
    std::string cluster_type = "accessory,strict,surplus";
    int process_num = 1;
    int pct_threshold = 0;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-i" || arg == "--input_path") {
            final_cluster_path = argv[++i];
        } else if (arg == "-o" || arg == "--output_path") {
            output_path = argv[++i];
        } else if (arg == "-f" || arg == "--fasta_path") {
            derep_fasta_path = argv[++i];
        } else if (arg == "-m" || arg == "--id_info_path") {
            id_info_path = argv[++i];
        } else if (arg == "-p" || arg == "--pre_cluster_path") {
            pre_cluster_path = argv[++i];
        } else if (arg == "-c" || arg == "--total_count") {
            total_count = std::stoi(argv[++i]);
        } else if (arg == "-u" || arg == "--thread_number") {
            process_num = std::stoi(argv[++i]);
        } else if (arg == "-t" || arg == "--cluster_type") {
            cluster_type = argv[++i];
        } else if (arg == "-a" || arg == "--pct_threshold") {
            pct_threshold = std::stoi(argv[++i]);
        }
        else if (arg == "-h" || arg == "--help") {        
            std::cout << "Thanks for using OrthoSLC! (version: " << __version__ << ")\n\n";
            std::cout << "Usage: Step9_write_clusters -i input_path -o output/ -f concatenated.fasta [options...]\n\n";
            std::cout << "options:\n";
            std::cout << "  -i or --input_path --------> <dir> path/to/input/final_cluster_file from Step 8\n";
            std::cout << "  -o or --output_path -------> <dir> path/to/output/directory\n";
            std::cout << "  -f or --fasta_path --------> <fasta> path/to/dereped_cated_fasta from Step 3\n";
            std::cout << "  -m or --id_info_path ------> <txt> path/to/output/id_info_table from Step 3\n";
            std::cout << "  -p or --pre_cluster_path --> <txt> path/to/pre_clustered_file from Step 3\n";
            std::cout << "  -c or --total_count -------> <int> amonut of genomes to analyze\n";
            std::cout << "  -a or --pct_threshold -----> <float> only write accessory clusters shared by >=n% (0<=n<100) of genomes, default: 0\n";
            std::cout << "  -t or --cluster_type ------> <txt> select from < accessory / strict / surplus >, separate by ',', all types if not specified\n";
            std::cout << "  -u or --thread_number -----> <int> thread number, default: 1\n";
            std::cout << "  -h or --help --------------> display this information\n";
            exit(0);
        }
    }

    // sanity checks
    if (final_cluster_path.empty() || output_path.empty() || derep_fasta_path.empty() ||
        id_info_path.empty() || pre_cluster_path.empty() || total_count <= 0) {
        return 1;
    }
    if (!fs::exists(final_cluster_path) || !fs::exists(derep_fasta_path) ||
        !fs::exists(id_info_path) || !fs::exists(pre_cluster_path) ||
        (pct_threshold < 0 || pct_threshold >= 100)) {
        std::cerr << "Error: invalid argument or file not found.\n";
        return 1;
    }
    fs::create_directories(output_path);

    // setup cluster type directories
    std::unordered_map<std::string,bool> types_to_write{
        {"accessory_cluster", false},
        {"strict_core", false},
        {"surplus_core", false}
    };
    std::stringstream ss(cluster_type);
    std::string token;
    while (std::getline(ss, token, ',')) {
        if (token == "accessory") { types_to_write["accessory_cluster"] = true;
            fs::create_directory(output_path + "/accessory_cluster"); }
        else if (token == "strict")  { types_to_write["strict_core"] = true;
            fs::create_directory(output_path + "/strict_core"); }
        else if (token == "surplus") { types_to_write["surplus_core"] = true;
            fs::create_directory(output_path + "/surplus_core"); }
    }



    // parallel load helper data
    ThreadPool pool(process_num);
    std::vector<std::future<void>> futures;

    std::unordered_map<std::string,std::string> id_cluster_dict;
    std::unordered_map<std::string,std::string> id_info_dict;

    // read final clusters
    std::vector<std::vector<std::string>> CLUSTERS;
    futures.emplace_back(
        pool.enqueue([&final_cluster_path, &CLUSTERS](){
            std::ifstream in(final_cluster_path);
            std::string line;
            while (std::getline(in, line)) {
                if (line.empty()) continue;
                std::vector<std::string> cluster;
                std::stringstream ls(line);
                std::string id;
                while (std::getline(ls, id, '\t')) cluster.push_back(id);
                CLUSTERS.push_back(std::move(cluster));
            }
        })
    );

    // load pre-cluster
    futures.emplace_back(pool.enqueue([&](){
        std::ifstream in(pre_cluster_path);
        std::string a_line;
        while (std::getline(in, a_line)) {
            std::vector<std::string> v;
            std::stringstream iss(a_line);
            std::string element;

            while (std::getline(iss, element, '\t')) {
                v.push_back(element);

                id_cluster_dict[element] = v[0];
            }
        }
    }));
    // load id info
    futures.emplace_back(pool.enqueue([&id_info_path, &id_info_dict](){
        std::ifstream in(id_info_path);
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            std::istringstream iss(line);
            std::string id, info;
            if (std::getline(iss, id, '\t') &&
                std::getline(iss, info, '\t')) {
                id_info_dict[id] = info;
            }
        }
    }));

    // load fasta
    std::unordered_map<std::string,std::string> fasta_dict;
    futures.emplace_back(pool.enqueue([&fasta_dict, &derep_fasta_path](){
        Reader reader(derep_fasta_path);
        std::string header, seq;
        while (reader.next(header, seq)) {
            
            auto pos = header.find(' ');
            std::string gid = (pos == std::string::npos ? header : header.substr(0, pos));
            
            // std::string gid = header;
            fasta_dict[gid] = std::move(seq);
        }
    }));

    for (auto & f : futures) f.wait();
    futures.clear();
    

    int threshold = static_cast<int>(std::ceil(pct_threshold / 100.0 * total_count));

    // write clusters
    for (auto & cluster : CLUSTERS) {
        futures.emplace_back(pool.enqueue([&, cluster](){
            std::unordered_set<std::string> species_set;
            for (auto & gid : cluster) {
                species_set.insert(gid.substr(0, gid.find('-')));
            }
            if ((int)species_set.size() < threshold) return;

            bool do_write = false;
            fs::path out_file;
            size_t sc = species_set.size();
            if (sc < (size_t)total_count && types_to_write["accessory_cluster"]) {
                out_file = fs::path(output_path) / "accessory_cluster" / (cluster[0] + ".fasta");
                do_write = true;
            } else if (cluster.size() == (size_t)total_count && types_to_write["strict_core"]) {
                out_file = fs::path(output_path) / "strict_core" / (cluster[0] + ".fasta");
                do_write = true;
            } else if (cluster.size() > (size_t)total_count && types_to_write["surplus_core"]) {
                out_file = fs::path(output_path) / "surplus_core" / (cluster[0] + ".fasta");
                do_write = true;
            }
            if (!do_write) return;

            Writer writer(out_file.string(), WRITER_BUF);
            for (auto & gid : cluster) {
                // std::cout << gid << id_info_dict[gid] << std::endl;
                // std::cout << gid << fasta_dict[id_cluster_dict[gid]] << std::endl;
                writer.write(gid + " " + id_info_dict[gid], fasta_dict[id_cluster_dict[gid]]);
            }
            writer.flush();
        }));
    }
    for (auto & f : futures) f.wait();

    return 0;
}
