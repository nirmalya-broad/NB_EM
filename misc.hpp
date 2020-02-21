#ifndef _MISC_HPP
#define _MISC_HPP

#include <vector>

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}


std::vector<std::string> split(const std::string &s,
        char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

double do_mean(std::vector<long> lvec) {
    double ltotal = 0;
    for (long lval1 : lvec) {
        ltotal += (double)lval1;
    }

    return ((double)ltotal/lvec.size());
}

#endif
