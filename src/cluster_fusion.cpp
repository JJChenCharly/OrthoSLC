#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <filesystem>

namespace fs = std::filesystem;

// Split a line by tabs
static std::vector<std::string> split_tab(const std::string &line) {
    std::vector<std::string> toks;
    std::istringstream ss(line);
    std::string token;
    while (std::getline(ss, token, '\t')) {
        toks.push_back(std::move(token));
    }
    return toks;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: cluster_fusion <pre_cluster.txt> <final_cluster_dir>\n";
        return 1;
    }

    const fs::path pre_path = argv[1];
    const fs::path out_dir  = argv[2];
    const fs::path out_file = out_dir / "0.txt";

    // 1) Read pre‑cluster file into an unordered_map
    std::ifstream pre_in(pre_path);
    if (!pre_in) {
        std::cerr << "Error: cannot open pre‑cluster file: " << pre_path << "\n";
        return 1;
    }

    std::unordered_map<std::string, std::unordered_set<std::string>> pre_cluster;
    std::string line;
    while (std::getline(pre_in, line)) {
        if (line.empty()) continue;
        auto toks = split_tab(line);
        auto &bucket = pre_cluster[toks[0]];            // creates empty set if not exists
        for (auto &g : toks) bucket.insert(std::move(g));
    }
    pre_in.close();

    // 2) Read MCL output from stdin and fuse with pre‑clusters
    std::unordered_map<std::string, std::unordered_set<std::string>> fused;
    std::unordered_set<std::string> seen_in_mcl;

    while (std::getline(std::cin, line)) {
        if (line.empty()) continue;
        auto toks = split_tab(line);
        std::unordered_set<std::string> merged(toks.begin(), toks.end());
        for (auto &g : toks) {
            seen_in_mcl.insert(g);
            auto it = pre_cluster.find(g);
            if (it != pre_cluster.end()) {
                merged.insert(it->second.begin(), it->second.end());
            }
        }
        fused.emplace(toks[0], std::move(merged));
    }

    // 3) Carry over any genes never seen by MCL
    for (auto & [gene, bucket] : pre_cluster) {
        if (seen_in_mcl.find(gene) == seen_in_mcl.end()) {
            fused.emplace(gene, bucket);
        }
    }

    // 4) Make sure output directory exists
    if (!fs::exists(out_dir)) {
        fs::create_directories(out_dir);
    }

    // 5) Write fused clusters to `0.txt`
    std::ofstream out(out_file);
    if (!out) {
        std::cerr << "Error: cannot create output file: " << out_file << "\n";
        return 1;
    }

    for (auto & [cid, members] : fused) {
        bool first = true;
        for (auto &g : members) {
            if (!first) out << '\t';
            out << g;
            first = false;
        }
        out << '\n';
    }
    out.close();

    return 0;
}
