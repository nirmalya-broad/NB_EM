#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>

#include <boost/program_options.hpp>

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>

#include <dlib/optimization.h>

#include "misc.hpp"
#include "NBModel.hpp"


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
        //std::vector<double> prior_vec;
        double prior_vec[3];
        double p_vec[3];
        double r_vec[3];
        double mean_vec[3];

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

        std::vector<long> get_sub_vec(std::vector<long>& lvec, 
            int start_p, int end_p);
        double CD_LL_value_comp1(const column_vector& m);
        double CD_LL_value_comp2(const column_vector& m);

        const column_vector CD_LL_gradient_comp1(const column_vector& m);
        const column_vector CD_LL_gradient_comp2(const column_vector& m);

        double get_NB_density(long xi_val, double p_val, double r_val);

        void init_EM_params();
        void do_E_step();
        void do_M_step();
        void update_p_val(int k);
        void update_prior_val(int k);
        double get_LL();
        
};

double MixtureModelC::get_NB_density(long xi_val, double p_val, double r_val) {

    double lval1 = lgamma(xi_val + r_val) -lgamma(xi_val + 1) -lgamma(r_val) + 
            xi_val * log (p_val) + r_val * log(1 - p_val);
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
    double lprior_ln = log(lprior_val);

    double total_ll = 0;


    for (int j_ind = 0; j_ind < sample_vec.size(); j_ind++) {
        long xi_val = sample_vec[j_ind];

        double lval1 = lgamma(xi_val + r_val) -lgamma(xi_val + 1) 
            - lgamma(r_val) + xi_val * log (lpval) + r_val * log(1-lpval);
        double lval2 = lprior_ln + lval1;
        double lmem_prob = get_mem_prob(density_ind, j_ind);
        double LL_xi = lmem_prob * lval2;
        total_ll += LL_xi;
    }         

    return total_ll;
}

double MixtureModelC::CD_LL_value_comp1(const column_vector& m) {
    const double r_val = m(0);
    double lval1 = CD_LL_value_k(1, r_val); 
    return (lval1);
}

double MixtureModelC::CD_LL_value_comp2(const column_vector& m) {
    const double r_val = m(0);
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

    for (int j_ind = 0; j_ind < sample_vec.size(); j_ind++) {
        long xi_val = sample_vec[j_ind];

        double lval1 = boost::math::digamma(xi_val + r_val) - 
            boost::math::digamma(r_val) + log(1-lpval);
        double lmem_prob = get_mem_prob(density_ind, j_ind);
        double gradient_xi = lmem_prob * lval1;
        total_gradient += gradient_xi;
    }         


    //std::cout << "total_gradient: " << total_gradient << "\n";

    return total_gradient;
}

const column_vector MixtureModelC::CD_LL_gradient_comp1(const column_vector& m) {
    const double r_val = m(0);
    double lval1 = CD_LL_gradient_k(1, r_val);
    column_vector lval2 = {lval1};
    return (lval2);
}

const column_vector MixtureModelC::CD_LL_gradient_comp2(const column_vector& m) {
    const double r_val = m(0);
    double lval1 = CD_LL_gradient_k(2, r_val);
    column_vector lval2 = {lval1};
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
    return (lmem_prob);
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

std::vector<long> MixtureModelC::get_sub_vec(std::vector<long>& lvec, 
        int start_p, int end_p) {

    auto first = lvec.cbegin() + start_p;
    auto last = lvec.cbegin() + end_p + 1;
    std::vector<long> new_v(first, last);
    return new_v;

}
 
void MixtureModelC::init_EM_params() {

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

    double zero_frac_m0 = 0.9;
    int m0_count  = (int) floor(zero_count * zero_frac_m0);
    int m1_count = less_500_count - m0_count;
    int m2_count = sample_vec.size() - (m0_count + m1_count);

    std::cout << "model0_count: " << m0_count << "\n";
    std::cout << "model1_count: " << m1_count << "\n";
    std::cout << "model2_count: " << m2_count << "\n";
    std::cout << ".................\n";

    std::vector<long> sample_vec_s = sample_vec;
    std::sort(sample_vec_s.begin(), sample_vec_s.end());

    int m0_s = 0; 
    int m0_e = m0_count -1;
    std::vector<long> sample_model0 =  get_sub_vec(sample_vec_s, m0_s, m0_e);

    
    int m1_s = m0_count;
    int m1_e = m0_count + m1_count -1;
    std::vector<long> sample_model1 =  get_sub_vec(sample_vec_s, m1_s, m1_e);

    int m2_s = m0_count + m1_count;
    int m2_e = sample_vec_s.size() -1;
    std::vector<long> sample_model2 =  get_sub_vec(sample_vec_s, m2_s, m2_e);

    int sample_vec_s_size = sample_vec_s.size();
    std::cout << "sample_vec_s_size: " << sample_vec_s_size << "\n";

    prior_vec[0] = (double) m0_count / sample_vec_s_size;
    prior_vec[1] = (double) m1_count / sample_vec_s_size;
    prior_vec[2] = (double) m2_count / sample_vec_s_size;

    std::cout << "Reached here." << "\n";
    std::cout << "model0_prior: " << prior_vec[0] << "\n";
    std::cout << "model1_prior: " << prior_vec[1] << "\n";
    std::cout << "model2_prior: " << prior_vec[2] << "\n";
    std::cout << ".................\n";

    // We shall get the p_val and r_val for model_1 and model_2 by doing
    // another EM?
    
    std::cout << "sample_model1: " << sample_model1.size() << "\n";

    NBModelC model1(sample_model1);
    r_vec[1] = model1.find_r_val();
    p_vec[1] = model1.find_p_val();
    std::cout << "model1_r_val: " << r_vec[1] << "\n";
    std::cout << "model1_p_val: " << p_vec[1] << "\n";
    std::cout << "model1_mean_val: " << model1.find_mean_val() << "\n"; 
    std::cout << ".................\n";

    NBModelC model2(sample_model2);
    r_vec[2] = model2.find_r_val();
    p_vec[2] = model2.find_p_val(); 

    std::cout << "model2_r_val: " << r_vec[2] << "\n";
    std::cout << "model2_p_val: " << p_vec[2] << "\n";
    std::cout << "model2_mean_val: " << model2.find_mean_val() << "\n"; 
    std::cout << ".................\n";
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


void MixtureModelC::update_prior_val(int k) {

    double N_k = 0;
    for (int n = 0; n < sample_vec.size(); n++) {
        N_k += get_mem_prob(k, n);
    }

    int N = sample_vec.size();
    prior_vec[k] = N_k / N;
}

void MixtureModelC::update_p_val(int k) {
    
    double N_k = 0;
    for (int n = 0; n < sample_vec.size(); n++) {
        N_k += get_mem_prob(k, n);
    }

    double lval1 = 0;

    for (int n = 0; n < sample_vec.size(); n++) {
        double xn_val = sample_vec[n];
        double x_mem_prob = get_mem_prob(k, n);
        lval1 += ( xn_val * x_mem_prob);
    }

    double mean_k = lval1 / N_k;
    double r_val_k = r_vec[k];
    
    p_vec[k] = mean_k / (mean_k + r_val_k);
}

void MixtureModelC::do_M_step() {


    // Estimate the parameters using the BFGS in dlib
    auto value_comp1 = std::bind(&MixtureModelC::CD_LL_value_comp1, this, _1);
    auto gradient_comp1 = std::bind(&MixtureModelC::CD_LL_gradient_comp1, this, _1);

    // Get updated r_val for component 1
    
     const column_vector x_lower1 = {0.0000001};
     const column_vector x_upper1 = {1000000};

     column_vector starting_point1 = {r_vec[1]};
     find_max_box_constrained(dlib::bfgs_search_strategy(),
                                     dlib::objective_delta_stop_strategy(1e-9),
                                     value_comp1, gradient_comp1,
                                    starting_point1,
                                    x_lower1, x_upper1);
    r_vec[1] = starting_point1(0);
    // update p_val for component 1
    update_p_val(1);
    

    // Similarly update r_val for component 2
    
    auto value_comp2 = std::bind(&MixtureModelC::CD_LL_value_comp2, this, _1);
    auto gradient_comp2 = std::bind(&MixtureModelC::CD_LL_gradient_comp2, this, _1);
     const column_vector x_lower2 = {0.0000001};
     const column_vector x_upper2 = {1000000};
     column_vector starting_point2 = {r_vec[2]};
     find_max_box_constrained(dlib::bfgs_search_strategy(),
                                     dlib::objective_delta_stop_strategy(1e-9),
                                     value_comp2, gradient_comp2,
                                    starting_point2,
                                    x_lower1, x_upper1);
    r_vec[2] = starting_point2(0);
    // update p_val for component2
    update_p_val(2);


    // calculate the prior_val for all three component
    update_prior_val(0);
    update_prior_val(1);
    update_prior_val(2);

}

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
        double ll_n = log(l_n);
        ll += ll_n;
    }

    return (ll);
}


void MixtureModelC::main_func() {

    // Here would be the four main steps.

    double target_ll_change = 1e-9;
    double ll_change = 0;
    // 1. initialize the parameters
    init_EM_params();

    double old_ll = get_LL();

    std::cout << "It 0, LL val: " << old_ll << "\n";

    // 2. do ... while not converged in terms of likelihood
    int it_count = 0;
    do {
        it_count++;

        // 3a. E step
        do_E_step(); 

        // 3b. M step
        do_M_step();

        //4. Calculate the log linklihood
    double new_ll = get_LL();

    double mean_1 = (p_vec[1] * r_vec[1])/(1 - p_vec[1]);
    double mean_2 = (p_vec[2] * r_vec[2])/(1 - p_vec[2]);
    std::cout << "It " << it_count << ", LL val " << new_ll << "\n";
    std::cout << "prior_0: " << prior_vec[0] << ", prior_1: " << 
        prior_vec[1] << ", prior_2: " << prior_vec[2] << "\n";
    std::cout << "p_val_1: " << p_vec[1] << ", p_val_2: " << p_vec[2] << 
        ", r_val_1: " << r_vec[1]  << ", r_val_2: " << r_vec[2] << 
        ", mean_1: " << mean_1 << ", mean_2: " << mean_2 << "\n";
    
    std::cout << "................\n";

    ll_change = abs(new_ll - old_ll);
    old_ll = new_ll;

    } while(ll_change > target_ll_change);

    for (int n = 0; n < sample_vec.size(); n++) {
        std::cout << sample_vec[n] << ", g_0: " << get_mem_prob(0, n) <<
            ", g_1: " << get_mem_prob(1, n) << ", g_2: " << 
            get_mem_prob(2, n) << "\n";
    }

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
