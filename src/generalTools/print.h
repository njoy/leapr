
#ifndef LEAPR_GENERALTOOLS_PRINTRANGE
#define LEAPR_GENERALTOOLS_PRINTRANGE
#include <iostream>
#include <range/v3/all.hpp>

template <typename Range>
void printRange( Range range){
  ranges::for_each(range, [](auto x){
        std::cout << x << ' ';
    });
  std::cout << std::endl;
}
#endif
