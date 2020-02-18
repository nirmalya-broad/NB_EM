#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>

#include <boost/program_options.hpp>

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>

#include "misc.hpp"


namespace po = boost::program_options;

// This would optimize a mixture of two three densities; We shall derive the 
// parameters in the p,r format. Later we can convert them to the m,r format
// as required.

// A. one zero-point density and two negative binomials
// B. One zero-point density and one negative binomials
// C. A single negative binomial
// D. Hurdle distributions ?

class NBinomEMC {

    private:

        po::options_description desc;
        std::string infile_str;
        std::string outdir_str;
        std::string prefix_str;
        int density_count;
        long sample_count;
        std::vector<long> sample_vec;
        double* mem_prob_mat; 

    public:
       
        NBinomEMC();
        void print_help();
        bool parse_args(int argc, char* argv[]);
        void read_infile();
        void allocate_resources();
        void initialize();
        void main_func();
        void free_vars();
        double get_mem_prob(int density_ind, int sample_ind);
        void set_mem_prob(double mem_prob, int density_ind, int sample_ind);
        double M_step_density_mean (int mix_ind);
        double M_step_density_prob (int mix_ind);
        
};

NBinomEMC::NBinomEMC() { 
    
}

void NBinomEMC::free_vars() {
    // Release the resources allocated by new
    delete[]  mem_prob_mat;
}

void NBinomEMC::allocate_resources() {

    if (density_count == 0 || sample_count == 0) {
        std::string err_msg = "Invalid value(s); density_count: " + 
            std::to_string(density_count) + ", sample_count: " + 
            std::to_string(sample_count);
        throw std::runtime_error(err_msg);
    }

    mem_prob_mat = new double[density_count * sample_count];
}

void NBinomEMC::initialize() {
    read_infile();
    allocate_resources();
}

void NBinomEMC::read_infile() {

    std::ifstream inf(infile_str);
    std::string line1;     
    long temp_num1;   
    while (std::getline(inf, line1)) {

        std::vector<std::string> lparts = split(line1, '\t');
        std::string gap_str = lparts[2];
        temp_num1 = std::stol(gap_str);
        sample_vec.push_back(temp_num1);
    } 
    
    sample_count = sample_vec.size();
    std::cout << "sample_count: " << std::to_string(sample_count) << "\n";
}

double NBinomEMC::get_mem_prob(int density_ind, int sample_ind) {
    size_t actual_pos = density_ind * sample_count + sample_ind;
    double lmem_prob = mem_prob_mat[actual_pos];
}

void NBinomEMC::set_mem_prob(double mem_prob, int density_ind, 
    int sample_ind) {
    size_t actual_pos = density_ind * sample_count + sample_ind;
    mem_prob_mat[actual_pos] = mem_prob;
}

void NBinomEMC::print_help() {
    std::cout << desc << "\n";
    std::cout << "Usage: NBinomEM --infile <txt> --outdir <outdir>"
        " --prefix <prefix>  --density_c <density count>"
        "\n\n";
}

bool NBinomEMC::parse_args(int argc, char* argv[]) {

    bool all_set = true;

    desc.add_options()
        ("help,h", "produce help message")
        ("infile,i", po::value<std::string>(&infile_str), "Infile for the gap data.")
        ("outdir,o", po::value<std::string>(&outdir_str), "Output dir.")
        ("prefix,p", po::value<std::string>(&prefix_str), "Prefix str.")
        ("density_c,d", po::value(&density_count), "Density count.")
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

    if (vm.count("outdir")) {
        std::cout << "Outdir is set to " << outdir_str << "\n";
    } else {
        all_set = false;
        std::cout << "Error: outdir is not set.\n";
    }

    if (vm.count("prefix")) {
        std::cout << "Prefix is set to: " << prefix_str << "\n";
    } else {
        all_set = false;
        std::cout << "Error: Prefix is not set.\n";
    }

    if (vm.count("density_c")) {
        std::cout << "Density_c is set to: " << density_count << "\n";
    } else {
        all_set = false;
        std::cout << "Error: Density_c is not set.\n";
    }

    return all_set;
}

void NBinomEMC::main_func() {

}

int main(int argc, char** argv) {
    NBinomEMC nbemo;
    bool all_set = true;

    try {
        all_set = nbemo.parse_args(argc, argv);
    } catch(std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    } catch(...) {
        return 0;
    } 

    if (!all_set) {
        nbemo.print_help();
        return 0;
    }

    try {
        nbemo.initialize();
        nbemo.main_func();
        nbemo.free_vars();

    } catch(const std::runtime_error& e) {

        std::cerr << "error: "  << e.what() << "\n";
        return 1;
    }

    return 0;
}
