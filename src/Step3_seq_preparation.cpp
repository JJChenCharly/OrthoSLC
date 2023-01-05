#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>

namespace fs = std::filesystem;

int main(int argc, char** argv) {
    // parameter parsing
    std::string rr_ed_dir;
    std::string cated_fasta_pth;
    std::string seq_len_pth;

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
        else if (std::string(argv[i]) == "--seq_len_tbl" || std::string(argv[i]) == "-l")
        {
            seq_len_pth = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h")
        {
            std::cout << "Thanks for using OrthoSLC! (version: 0.1Beta)\n\n";
            std::cout << "Usage: Step3_seq_preparation -i input/ -c concatenated.fasta -l seq_len.txt\n\n";
            std::cout << "  -i or --input_path ---------> path/to/input/directory\n";
            std::cout << "  -c or --concatenated_fasta -> path/to/output/concatenated.fasta\n";
            std::cout << "  -l or --seq_len_tbl --------> path/to/output/sequence_length_table\n";
            std::cout << "  -h or --help ---------------> display this information\n";
            exit(0);
        }
    }

    //abs pth
    rr_ed_dir = fs::absolute(fs::path(rr_ed_dir));
    cated_fasta_pth = fs::absolute(fs::path(cated_fasta_pth));
    seq_len_pth = fs::absolute(fs::path(seq_len_pth));

    std::ofstream cated_fasta(cated_fasta_pth, std::ios::trunc);
    std::ofstream seq_len(seq_len_pth, std::ios::trunc);

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
        seq_len << a_line.substr(1, 11) << "\t";
        
        while (getline(fasta_in_put, a_line)) {
            cated_fasta << a_line << "\n";

            if (a_line[0] == '>') { // if the row of id 
                // save previous seq len
                seq_len << a_seq.length() << "\n";
                a_seq = "";

                seq_len << a_line.substr(1, 11) << "\t";
            } else {
                a_seq = a_seq + a_line;
            }
        }
        
        // last length info
        seq_len << a_seq.length() << "\n";

        fasta_in_put.close();

    }

    cated_fasta.close();
    seq_len.close();

    return 0;
}