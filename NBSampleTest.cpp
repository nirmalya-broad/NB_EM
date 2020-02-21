#include <iostream>
#include <vector>
#include <functional>

#include <dlib/optimization.h>
#include <dlib/global_optimization.h>


#include "NBSample.hpp"

int main() {

    std::vector<long> lsample{16, 2, 77, 29,  101, 525, 673, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 16, 2, 77, 29, 101, 525};

    NBSampleC nbso(lsample);
    double r_val = nbso.find_r_val();
    std::cout << "r_val: " << r_val << "\n";


    
}
