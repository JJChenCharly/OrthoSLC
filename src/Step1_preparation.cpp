#include <iostream>
#include <filesystem>
#include <string>
#include <fstream>

#include "ThreadPool.h"

namespace fs = std::filesystem;

int main(int argc, char** argv) {   
    // parameter parsing ----
    std::string raw_annotated_dir_path;
    std::string pre_res_path;

    for (int i = 1; i < argc; ++i)
    {
        if (std::string(argv[i]) == "--input_path" || std::string(argv[i]) == "-i")
        {
            raw_annotated_dir_path = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--output_path" || std::string(argv[i]) == "-o")
        {
            pre_res_path = argv[i + 1];
        }
        else if (std::string(argv[i]) == "--help" || std::string(argv[i]) == "-h")
        {
            std::cout << "Thanks for using OrthoSLC! (version: " << __version__ << ")\n\n";
            std::cout << "Usage: Step1_preparation -i input/ -o output.txt\n\n";
            std::cout << "  -i or --input_path ---> <dir> path/to/input/directory\n";
            std::cout << "  -o or --output_path --> <txt> path/to/output_file\n";
            std::cout << "  -h or --help ---------> display this information\n";
            exit(0);
        }
    }

    // if file path exist
    if (!(fs::exists(raw_annotated_dir_path))) {
        std::cerr << "Error: path provided to '-i or --input_path' does not exist. 路径不存在\n";
        exit(0);
    }

    fs::path parent_pre_res_path = fs::path(pre_res_path).parent_path();
    if (!(fs::exists(parent_pre_res_path))) {
        std::cerr << "Error: parent path provided to '-o or --output_path' does not exist. 路径不存在\n";
        exit(0);
    }

    std::ofstream pre_res_output(pre_res_path, std::ios::trunc); // m = overwrite

    int start_c = 0;
    std::string strain_naam;

    for (const auto& entry : fs::directory_iterator(raw_annotated_dir_path)) { // range iteration to get abs path
        fs::path current_path = entry.path();
        std::string abs_path = fs::absolute(current_path);
        
        strain_naam = current_path.filename();
        size_t lastindex = strain_naam.find_last_of(".");
        strain_naam = strain_naam.substr(0, lastindex);

        std::string str_start_c = std::to_string(start_c);
        if (str_start_c.length() == 4) {
            pre_res_output << '0' << start_c << "\t";
        } else {
            pre_res_output << start_c << "\t";
        }

        pre_res_output << strain_naam << "\t";

        pre_res_output << abs_path << "\n";

        start_c++;
    }
    
    pre_res_output.close();

    return 0;
}