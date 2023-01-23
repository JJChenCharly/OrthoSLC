#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>

namespace fs = std::filesystem;

int main(int argc, char** argv) {
    // parameter parsing
    std::string rr_ed_dir;
    std::string cated_fasta_pth;

    for (int i = 1; i < argc; ++i)
    {
        if (std::string(argv[i]) == "--input_path" || std::string(argv[i]) == "-i")
        {
            rr_ed_dir = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--concatenated_fasta" || std::string(argv[i]) == "-c")
        {
            cated_fasta_pth = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h")
        {
            std::cout << "Thanks for using OrthoSLC! (version: 0.1.1)\n\n";
            std::cout << "Usage: Step3_seq_concat -i input/ -c concatenated.fasta\n\n";
            std::cout << "  -i or --input_path ---------> path/to/input/directory\n";
            std::cout << "  -c or --concatenated_fasta -> path/to/output/concatenated.fasta\n";
            std::cout << "  -h or --help ---------------> display this information\n";
            exit(0);
        }
    }
    // if file path exist
    if (!(fs::exists(rr_ed_dir))) {
        std::cerr << "Error: path provided to '-i or --input_path' do not exist. 路径不存在\n";
        exit(0);
    }

    fs::path parent_cated_fasta_pth = fs::path(cated_fasta_pth).parent_path();
    if (!(fs::exists(parent_cated_fasta_pth))) {
        std::cerr << "Error: parent path provided to '-c or --concatenated_fasta' do not exist. 路径不存在\n";
        exit(0);
    }

    //abs pth
    rr_ed_dir = fs::absolute(fs::path(rr_ed_dir));
    cated_fasta_pth = fs::absolute(fs::path(cated_fasta_pth));

    std::ofstream cated_fasta(cated_fasta_pth, std::ios::trunc);

    std::string a_line;
    std::string a_seq = "";
    
    std::string in_p;
    // Iterate over the entries in the rr_ed_dir
    for (const auto& entry : fs::directory_iterator(rr_ed_dir)) {
        in_p = entry.path();
        std::ifstream fasta_in_put(in_p);

        // first line
        getline(fasta_in_put, a_line);
        cated_fasta << a_line << "\n";
        
        while (getline(fasta_in_put, a_line)) {
            cated_fasta << a_line + "\n";
        }

        fasta_in_put.close();

    }

    cated_fasta.close();

    return 0;
}