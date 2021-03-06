#ifndef _NB_SAMPLE_HPP
#define _NB_SAMPLE_HPP

#include <vector>
#include <cmath>
#include <numeric>

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>

#include <dlib/optimization.h>


#include "misc.hpp" 

typedef dlib::matrix<double,0,1> column_vector;
using namespace std::placeholders;

class NBModelC {

    public:
        NBModelC (std::vector<long> lvec) : sample_vec(lvec){

        }

        double value(const column_vector &m) {

            double r_val = m(0);
            double lmean = do_mean(sample_vec);
            double lpval = lmean / (lmean + r_val);
            double total_ll = 0;

            for (int j_ind = 0; j_ind < sample_vec.size(); j_ind++) {
                double xi_val = (double)sample_vec[j_ind];

                double lval1 = lgamma(xi_val + r_val) -lgamma(xi_val + 1) -
                    lgamma(r_val) + xi_val * log (lpval) + 
                    r_val * log(1 - lpval);
                total_ll += lval1;
            }

            return (total_ll);
        }


        const column_vector gradient(const column_vector &m) {
           
            double r_val = m(0);
            double lmean = do_mean(sample_vec);
            double lpval = lmean / (lmean + r_val);
         
            double total_derivative = 0;

            for (int j_ind = 0; j_ind < sample_vec.size(); j_ind++) {
                double xi_val = sample_vec[j_ind];

                double lval1 = boost::math::digamma(xi_val + r_val) -
                    boost::math::digamma(r_val) + log(1-lpval);
                
                total_derivative += lval1;
            }
            const column_vector final_grad = {total_derivative};
            return final_grad;

        }


        double get_r_val_0 () {
            int n = sample_vec.size();
            double lmean = do_mean(sample_vec);

            double lsum = 0;
            for (int j = 0; j < n; j++) {
                lsum += pow((sample_vec[j] - lmean), 2.0);
            }
            double lvar = lsum / (n-1);

            double r_val_0 = 0;
            if (lvar > lmean) {
                r_val_0 = pow(lmean, 2.0) / (lvar - lmean);
            } else {
                r_val_0 = 100;
            }
            return r_val_0;
        }

        double find_r_val() {
        
            if (r_val_ready) {
                return (final_r_val);
            }

            auto value_funct = std::bind(&NBModelC::value, this, _1);
            auto gradient_funct = std::bind(&NBModelC::gradient, this, _1);
            double r_val_0 = get_r_val_0();
            column_vector starting_point = {r_val_0};
            const column_vector x_lower1 = {0.0000001};
            const column_vector x_upper1 = {1000000};


            find_max_box_constrained(dlib::bfgs_search_strategy(),
                                     dlib::objective_delta_stop_strategy(1e-9),
                                     value_funct, gradient_funct, 
                                    starting_point,
                                    x_lower1, x_upper1);

            final_r_val = starting_point(0);
            r_val_ready = true;
            return final_r_val;
        }

        double find_p_val() {

            if (p_val_ready) {
                return (final_p_val);
            }
            double lrval = find_r_val();
 
            double lmean = do_mean(sample_vec);
            final_p_val = lmean / (lmean + lrval);
            p_val_ready = true;
            return (final_p_val);
            
        }

        double find_mean_val() {

            double lmean = do_mean(sample_vec);
            return (lmean);
        }


    private:
        std::vector<long> sample_vec;
        double final_p_val;
        double final_r_val;
        bool p_val_ready = false;
        bool r_val_ready = false;

};

#endif
