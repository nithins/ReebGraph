#include <iostream>
#include <sstream>

#include "../DisjointSets.hpp"

using namespace contourtree;

inline void __dump_error__(const std::string & s)
{std::cerr << s << std::endl;throw std::runtime_error(s);}

#define ENSURES(cond) \
  if(!(cond)) \
  for(std::stringstream ss ; true ; __dump_error__(ss.str())) \
  ss<<"Failed to ensure condition " << #cond <<"\n" \
    <<"at ("<<__FILE__<<","<<__func__<<","<<__LINE__<<") \n "



void testDisjointSets() {
    int numElements = 128;
    int numInSameSet = 16;
    DisjointSets<int64_t> ds(numElements);
    int set1, set2;

    for (int k = 1; k < numInSameSet; k *= 2) {
        for (int j = 0; j + k < numElements; j += 2 * k) {
            set1 = ds.find(j);
            set2 = ds.find(j + k);
            ds.merge(set1, set2);
        }
    }

    for (int i = 0; i < numElements; i++) {
//        std::cout << ds.find(i) << "*";
//        if (i % numInSameSet == numInSameSet - 1)
//            std::cout << "\n";

        ENSURES(ds.find(i) == numInSameSet*int(i/numInSameSet))
                << " ds.find(i)=" << ds.find(i)
                << " int(i/numInSameSet)" << int(i/numInSameSet);
    }
//    std::cout << "\n";
}

int main() {
    testDisjointSets();
}

