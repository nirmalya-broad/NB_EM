#include <iostream>
#include <vector>
#include <functional>


#include "NBModel.hpp"

int main() {

    std::vector<long> lsample{16, 2, 77, 29,  101, 525, 673, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 16, 2, 77, 29, 101, 525};

    NBModelC nbmo(lsample);
    double r_val = nbmo.find_r_val();
    std::cout << "r_val: " << r_val << "\n";


    
}
