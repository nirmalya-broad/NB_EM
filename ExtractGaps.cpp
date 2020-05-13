#include <iostream>
#include <fstream>
#include <string>

#include <boost/program_options.hpp>

#include "misc.hpp"
namespace po = boost::program_options;


class ExtractGapsC {

    private:
        po::options_description desc;
        std::string infile_str;
        std::string outfile_str;
        std::string strand_str;
        long start_pos;
        long end_pos;

    public:
        void print_help();
        void process_gap_file();
        bool parse_args(int argc, char* argv[]);

};

void ExtractGapsC::print_help() {
    std::cout << desc << "\n";
    std::cout << "Usage: ExtractGaps --infile <infile> --outfile <outfile>"
        " --gene_strand <gene_strand>  --start_pos <start_pos>"
        " --end_pos <end_pos>"
        "\n\n";
}

bool ExtractGapsC::parse_args(int argc, char* argv[]) {
    bool all_set = true;

    desc.add_options()
        ("help,h", "produce help message")
        ("infile,i", po::value<std::string>(&infile_str), "Infile for the gap data.")
        ("outfile,o", po::value<std::string>(&outfile_str), "Outfile for the gap data.")
        ("gene_strand,g", po::value<std::string>(&strand_str), "Strand of the genome.")
        ("start_pos,s", po::value(&start_pos), "Starting position")
        ("end_pos,e", po::value(&end_pos), "Ending position")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        return 0;
    } else {
    }

    if (vm.count("infile")) {
        std::cout << "Infile is set to: " << infile_str << "\n";
    } else {
        all_set = false;
        std::cout << "Error: infile is not set.\n";
    }

    if (vm.count("outfile")) {
        std::cout << "Outfile is set to: " << outfile_str << "\n";
    } else {
        all_set = false;
        std::cout << "Error: outfile is not set.\n";
    }

    if (vm.count("gene_strand")) {
        std::cout << "Gene_strand is set to: " << strand_str << "\n";
    } else {
        all_set = false;
        std::cout << "Error: Gene_strand is not set.\n";
    }

    if (vm.count("start_pos")) {
        std::cout << "Start_pos is set to: " << start_pos << "\n";
    } else {
        all_set = false;
        std::cout << "Error: Start is not set.\n";
    }

    if (vm.count("end_pos")) {
        std::cout << "End_pos is set to: " << end_pos << "\n";
    } else {
        all_set = false;
        std::cout << "Error: End_pos is not set.\n";
    }

    return all_set;
}


void ExtractGapsC::process_gap_file() {

    std::ifstream inf(infile_str);
    std::ofstream outf(outfile_str);

    std::string line1;

    long lcount = 0;
     while (std::getline(inf, line1)) {
    
        char ldelim;
        unsigned long lstart;
        unsigned long lend;
        extract_nums(line1, ldelim, lstart, lend);

        std::string lstrand_str(1, ldelim);

        if (lstrand_str == strand_str) {
            
            if ((start_pos <= lstart && lstart <= end_pos) ||
                (start_pos <= lend && lend <= end_pos)) {
                outf << line1 << "\n";
            }
        }

        lcount++;
    }
}

int main(int argc, char** argv) {
    ExtractGapsC ego;
    bool all_set = true;

    try {
        all_set = ego.parse_args(argc, argv);
    } catch(std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    } catch(...) {
        return 0;
    }

    if (!all_set) {
        ego.print_help();
        return 0;
    }

    try {
        ego.process_gap_file();

    } catch(const std::runtime_error& e) {

        std::cerr << "error: "  << e.what() << "\n";
        return 1;
    }

    return 0;
}

