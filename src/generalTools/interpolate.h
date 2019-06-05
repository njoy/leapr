#include <range/v3/all.hpp>
#include <iostream>

template <typename RangeZip, typename Float>
Float interpolate( RangeZip xyRange, Float x ){
  auto xVec = xyRange | ranges::view::keys;
  auto yVec = xyRange | ranges::view::values;
  std::cout << x << std::endl;

}





