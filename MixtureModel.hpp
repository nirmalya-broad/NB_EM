#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>

#include <boost/program_options.hpp>

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>

#include <dlib/optimization.h>
#include <dlib/global_optimization.h>

#include "misc.hpp"


namespace po = boost::program_options;

// This would optimize a mixture of two three densities; We shall derive the 
// parameters in the p,r format. Later we can convert them to the m,r format
// as required.

// A. one zero-point density and two negative binomials
// B. One zero-point density and one negative binomials
// C. A single negative binomial
// D. Hurdle distributions ?

typedef dlib::matrix<double,0,1> column_vector;
using namespace std::placeholders;

class MixtureModelC {

    private:

        po::options_description desc;
        std::string infile_str;
        std::string outdir_str;
        std::string prefix_str;
        int density_count;
        long sample_count;
        std::vector<long> sample_vec;
        double* mem_prob_mat; 

        // For our model with one zero density + two NB density; the r_val[0]
        // would be undefined. 
        std::vector<double> prior_vec;
        std::vector<double> p_vec;
        std::vector<double> r_vec;

        double CD_LL_value_k(int density_ind, double r_val);
        double CD_LL_gradient_k(int density_ind, double r_val);

    public:
       
        MixtureModelC();
        void print_help();
        bool parse_args(int argc, char* argv[]);
        void read_infile();
        void allocate_resources();
        void initialize();
        void main_func();
        void free_vars();
        double get_mem_prob(int density_ind, int sample_ind);
        void set_mem_prob(double mem_prob, int density_ind, int sample_ind);

        double CD_LL_value_comp1(const column_vector& m);
        double CD_LL_value_comp2(const column_vector& m);

        const column_vector CD_LL_gradient_comp1(const column_vector& m);
        const column_vector CD_LL_gradient_comp2(const column_vector& m);

        double get_NB_density(long xi_val, double p_val, double r_val);

        void do_E_step();
        double get_LL();

        auto value_comp1;
        auto value_comp2;
        auto gradient_comp1;
        auto gradient_compe2;
        
};

double MixtureModelC::get_LL() {
    double ll = 0;

    for (int n = 0; n < sample_count; n++) {

        long xn_val = sample_vec[n];
        double comp_0_L = 0;

        if (xn_val == 0) {
            comp_0_L = 1;
        } else {
            comp_0_L = 0;
        }

        double comp_0_u = prior_vec[0] * comp_0_L;

        double comp_1_L = get_NB_density(xn_val, p_vec[1], r_vec[1]); 
        double comp_1_u = prior_vec[1] * comp_1_L;

        double comp_2_L = get_NB_density(xn_val, p_vec[2], r_vec[2]);
        double comp_2_u = prior_vec[2] * comp_2_L; 
       
        double l_n = comp_0_u + comp_1_u + comp_2_u;
        double ll_n = log(l_n)
        ll += ll_n;
    }

    return (ll);
}

void MixtureModelC::do_E_step() {

    for (int n = 0; n < sample_count; n++) {

        long xn_val = sample_vec[n];
        double comp_0_L = 0;

        if (xn_val == 0) {
            comp_0_L = 1;
        } else {
            comp_0_L = 0;
        }

        double comp_0_u = prior_vec[0] * comp_0_L;

        double comp_1_L = get_NB_density(xn_val, p_vec[1], r_vec[1]); 
        double comp_1_u = prior_vec[1] * comp_1_L;

        double comp_2_L = get_NB_density(xn_val, p_vec[2], r_vec[2]);
        double comp_2_u = prior_vec[2] * comp_2_L; 
       
        double total_d = comp_0_u + comp_1_u + comp_2_u;

        double comp_0_n = comp_0_u / total_d;
        double comp_1_n = comp_1_u / total_d;
        double comp_2_n = comp_2_u / total_d; 
    
        set_mem_prob(comp_0_n, 0, n);    
        set_mem_prob(comp_1_n, 1, n);    
        set_mem_prob(comp_2_n, 2, n);    
    } 
}

double MixtureModelC::get_NB_density(long xi_val, double p_val, double r_val) {

    double lval1 = lgamma(xi_val + r_val) -lgamma(xi_val + 1) -lgamma(r_val) + 
            xi_val * log (lpval) + r_val * log(1-lpval);
    double lval2 = exp(lval1);

    return (lval2);

}

double MixtureModelC::CD_LL_value_k(int density_ind, double r_val) {
    
    if (!(density_ind == 1 || density_ind == 2)) {
        std::string err_str = "Density index and function mismathch.";
        throw std::runtime_error(err_str);
    }

    double lpval = p_vec[density_ind];
    double lprior_val = prior_vec[density_ind];
    double lprior_ln = math::log(lprior_val);

    double total_ll = 0;

    for (int j_ind = 0; j_ind < sample_vec.size(); j_ind) {
        long xi_val = sample_vec[j_ind];

        lval1 = lgamma(xi_val + r_val) -lgamma(xi_val + 1) -lgamma(r_val) + 
            xi_val * log (lpval) + r_val * log(1-lpval);
        lval2 = lprior_ln + lval1;
        double lmem_prob = get_mem_prob(density_ind, j_ind);
        double LL_xi = lmem_prob * lval2;
        total_ll += LL_xi;
    }         

    return total_ll;
}

double MixtureModelC::CD_LL_value_comp1(const column_vector& m) {
    const r_val = m(0);
    double lval1 = CD_LL_value_k(1, r_val); 
    return (lval1);
}

double MixtureModelC::CD_LL_value_comp2(const column_vector& m) {
    const r_val = m(0);
    double lval1 = CD_LL_value_k(2, r_val); 
    return (lval1);
}

double MixtureModelC::CD_LL_gradient_k(int density_ind, double r_val) {
    
    if (!(density_ind == 1 || density_ind == 2)) {
        std::string err_str = "Density index and function mismathch.";
        throw std::runtime_error(err_str);
    }

    double lpval = p_vec[density_ind];

    double total_gradient = 0;

    for (int j_ind = 0; j_ind < sample_vec.size(); j_ind) {
        long xi_val = sample_vec[j_ind];

        lval1 = boost::math::digamma(xi_val + r_val) - 
            boost::math::digamma(r_val) + log(1-lpval);
        double lmem_prob = get_mem_prob(density_ind, j_ind);
        double gradient_xi = lmem_prob * lval1;
        total_gradient += gradient_xi;
    }         

    return total_gradient;
}

const column_vector MixtureModelC::CD_LL_gradient_comp1(const column_vector& m) {
    const r_val = m(0);
    double lval1 = CD_LL_gradient_k(1, r_val);
    column_vector lval2 = {lval1};
    return (lval2);
}

const column_vector MixtureModelC::CD_LL_gradient_comp2(const column_vector& m) {
    const r_val = m(0);
    double lval1 = CD_LL_gradient_k(2, r_val);
    column_vector = {lval1};
    return (lval2);
}

MixtureModelC::MixtureModelC() { 
        
}

void MixtureModelC::free_vars() {
    // Release the resources allocated by new
    delete[]  mem_prob_mat;
}

void MixtureModelC::allocate_resources() {

    if (density_count == 0 || sample_count == 0) {
        std::string err_msg = "Invalid value(s); density_count: " + 
            std::to_string(density_count) + ", sample_count: " + 
            std::to_string(sample_count);
        throw std::runtime_error(err_msg);
    }

    mem_prob_mat = new double[density_count * sample_count];
}

void MixtureModelC::initialize() {
    read_infile();
    allocate_resources();
}

void MixtureModelC::read_infile() {

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

double MixtureModelC::get_mem_prob(int density_ind, int sample_ind) {
    size_t actual_pos = density_ind * sample_count + sample_ind;
    double lmem_prob = mem_prob_mat[actual_pos];
}

void MixtureModelC::set_mem_prob(double mem_prob, int density_ind, 
    int sample_ind) {
    size_t actual_pos = density_ind * sample_count + sample_ind;
    mem_prob_mat[actual_pos] = mem_prob;
}

void MixtureModelC::print_help() {
    std::cout << desc << "\n";
    std::cout << "Usage: NBinomEM --infile <txt> --outdir <outdir>"
        " --prefix <prefix>  --density_c <density count>"
        "\n\n";
}

bool MixtureModelC::parse_args(int argc, char* argv[]) {

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

void MixtureModelC::main_func() {
    auto value_comp1 = std::bind(&MixtureModelC::CD_LL_value_comp1, this, _1);
    auto value_comp2 = std::bind(&MixtureModelC::CD_LL_value_comp2, this, _1);

    auto gradient_comp1 = std::bind(&MixtureModelC::CD_LL_gradient_comp1, this, _1);
    auto gradient_comp1 = std::bind(&MixtureModelC::CD_LL_gradient_comp2, this, _1);


    // Here would be the four main steps.


}

int main(int argc, char** argv) {
    MixtureModelC nbemo;
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
