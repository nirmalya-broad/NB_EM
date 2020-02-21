#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>

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

typedef matrix<double,0,1> column_vector;

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

        // For our model with one zero density + two NB density; the r_val[0]
        // would be undefined. 
        std::vector<double> prior_vec;
        std::vector<double> p_vec;
        std::vector<double> r_vec;

        double get_nb_ll_by_k(int density_ind, double r_val);
        double get_nb_ll_derivative_by_k(int density_ind, double r_val);

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

        double get_nb_ll(const column_vector &m);
        double get_nb_ll_derivative(const column_vector &m);

        double get_mix_ll_by_comp1(const column_vector& m);
        double get_mix_ll_by_comp2(const column_vector& m);

        double get_mix_ll_derivative_by_comp1(const column_vector& m);
        double get_mix_ll_derivative_by_comp2(const column_vector& m);

        double get_NB_density(long xi_val, double p_val, double r_val);

        void do_E_step();
        double get_LL();
        
};

double NBinomEMC::get_LL() {
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

std::vector<long> NBinomEMC::get_sub_vec(std::vector<long> ori_vec, 
        int start_pos, int end_pos) {

    std::vector<long> dest_vec;

    for (for j = start_pos; j <= end_pos; j++) {
        dest_vec.push_back(ori_vec[j]);
    return (dest_vec);
}

double NBinomEMC::get_nb_ll(const column_vector &m) {

    double lr_val = m(0);

    auto n = sample_vec.size();
    double lmean = 0.0;
    if ( n != 0) {
        lmean = accumulate(sample_vec.begin(), sample_vec.end(), 0.0) / n;
    }

    double total_ll = 0;

    double lpval = lmean / (lmean + lr_val);

    for (int j_ind = 0; j_ind < sample_vec.size(); j_ind) {
        long xi_val = sample_vec[j_ind];

        lval1 = lgamma(xi_val + r_val) -lgamma(xi_val + 1) -lgamma(r_val) +
            xi_val * log (lpval) + r_val * log(1 - lpval);
        total_ll += lval1;
    }

}


double NBinomEMC::get_nb_ll_derivative(const column_vector &m) {

}



std::map<std::string, double> NBinomEMC::get_NB_params(std::vector<long> 
    lsample) {

    // Get initial set of values by method of moment as done by MASS R package

    auto n = lsample.size(); 
    double lmean = 0.0;
    if ( n != 0) {
        lmean = accumulate( lsample.begin(), lsample.end(), 0.0) / n; 
    }

    double lsum = 0;
    for (int j = 0; j < n; j++) {
        lsum += pow((lsample[j] - lmean), 2.0);
    }
    double lvar = lsum / (n-1);

    double r_val_0 = 0;
    if (lvar > lmean) {
        r_val_0 = pow(lmean, 2.0) / (lvar - lmean);
    } else {
        r_val_0 = 100;
    }

    
    // We now shall get the r_val using optimization.
    column_vector starting_point = {r_val_0};

    double r_val = find_max_box_constrained(bfgs_search_strategy(),  
                             objective_delta_stop_strategy(1e-9),  
                             rosen, rosen_derivative, starting_point, 0.01, 1000); 

}

void NBinomEMC::init_EM_params() {

    // Take the data and make a rough estimate from sample_vec;

    // 1. Get the count of 0
    // 2. Get the count of less than 500

    int zero_count = 0;
    int less_500_count = 0;
    for (long lval1 : sample_vec) {
        if (0 == lval1) {
            zero_count++;
        }

        if (lval1 < 500) {
            less_500_count++;
        }
    }

    double zero_frac_model_0 = 0.9;
    int model_0_count  = (int) floor(zero_count * zero_frac_model_0);
    int model_1_count = less_500_count - model_0_count;
    int model_2_count = sample_vec.size() - (model_0_count + model_1_count);

    std::vector<long> sample_vec_s = sample_vec;
    std::sort(sample_vec_s.begin(), sample_vec_s.end());

    std::vector<long> sample_model_0 =  get_sub_vec(sample_vec_s, 0, 
        (model_0_count -1));
    std::vector<long> sample_model_1 = get_sub_vec(sample_vec_s, model_0_count, 
        (model_0_count + model_1_count -1));
    std::vector<long> sample_model_2 = get_sub_vec(sample_vec_s, 
        (model_0_count + model_1_count), (sample_vec_s.size() -1);


    p_vec[0] = (double) model_0_count / sample_vec_s.size();

    // We shall get the p_val and r_val for model_1 and model_2 by doing
    // another EM?

    std::map<std::string, double> model_1_vals = get_NB_params(sample_model_1);
    p_vec[1] = model_1_vals["p_val"];
    r_vec[1] = model_1_vals["r_val"];
    
    std::map<std::string, double> model_2_vals = get_NB_params(sample_model_2);
    p_vec[2] = model_2_vals["p_val"];
    r_vec[2] = model_2_vals["r_val"];
    

}


void NBinomEMC::do_E_step() {

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

double NBinomEMC::get_NB_density(long xi_val, double p_val, double r_val) {

    double lval1 = lgamma(xi_val + r_val) -lgamma(xi_val + 1) -lgamma(r_val) + 
            xi_val * log (lpval) + r_val * log(1-lpval);
    double lval2 = exp(lval1);

    return (lval2);

}

double NBinomEMC::get_nb_ll_by_k(int density_ind, double r_val) {
    
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

double NBinomEMC::get_mix_ll_by_comp1(const column_vector& m) {
    const r_val = m(0);
    double lval1 = get_nb_ll_by_k(1, r_val); 
    return (lval1);
}

double NBinomEMC::get_mix_ll_by_comp2(const column_vector& m) {
    const r_val = m(0);
    double lval1 = get_nb_ll_by_k(2, r_val); 
    return (lval1);
}

double NBinomEMC::get_nb_ll_derivative_by_k(int density_ind, double r_val) {
    
    if (!(density_ind == 1 || density_ind == 2)) {
        std::string err_str = "Density index and function mismathch.";
        throw std::runtime_error(err_str);
    }

    double lpval = p_vec[density_ind];

    double total_derivative = 0;

    for (int j_ind = 0; j_ind < sample_vec.size(); j_ind) {
        long xi_val = sample_vec[j_ind];

        lval1 = boost::math::digamma(xi_val + r_val) - 
            boost::math::digamma(r_val) + log(1-lpval);
        double lmem_prob = get_mem_prob(density_ind, j_ind);
        double derivative_xi = lmem_prob * lval1;
        total_derivative += derivative_xi;
    }         

    return total_derivative;
}

double NBinomEMC::get_mix_ll_derivative_by_comp1(const column_vector& m) {
    const r_val = m(0);
    double lval1 = get_nb_ll_derivative_by_k(1, r_val);
    return (lval1);
}

double NBinomEMC::get_mix_ll_derivative_by_comp2(const column_vector& m) {
    const r_val = m(0);
    double lval1 = get_nb_ll_derivative_by_k(2, r_val);
    return (lval1);
}

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
