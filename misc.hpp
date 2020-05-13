#ifndef _MISC_HPP
#define _MISC_HPP

#include <string>
#include <cstring>
#include <vector>

void inplace_reverse(char * str)
{
  if (str)
  {
    char * end = str + strlen(str) - 1;

    // swap the values in the two given variables
    // XXX: fails when a and b refer to same memory location
#   define XOR_SWAP(a,b) do\
    {\
      a ^= b;\
      b ^= a;\
      a ^= b;\
    } while (0)

    // walk inwards from both ends of the string, 
    // swapping until we get to the middle
    while (str < end)
    {
      XOR_SWAP(*str, *end);
      str++;
      end--;
    }
#   undef XOR_SWAP
  }
}

void extract_nums (std::string ori_str1, char& delim, unsigned long& start_pos, 
        unsigned long& end_pos) {

    const char * ori_str = ori_str1.c_str();

    int ori_len = strlen(ori_str);
    char num_rev_arr[6][25];
    int tab_num = 0;
    int i = 0;
    int j = strlen(ori_str) -1;
    for (; j >=0, tab_num < 5; j--) {
        if (ori_str[j] == '\t') {
            num_rev_arr[tab_num][i] = '\0';
            i = 0;
            tab_num++;
        } else {
            num_rev_arr[tab_num][i] = ori_str[j];
            i++;
        }
        
    }

    inplace_reverse(num_rev_arr[0]);
    inplace_reverse(num_rev_arr[1]);
    
    
    start_pos = strtoul(num_rev_arr[1],NULL,10);
    end_pos = strtoul(num_rev_arr[0],NULL,10);

    delim = num_rev_arr[4][0];
    
}


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

