#include <range/v3/all.hpp>

template <typename Range, typename Callable >
auto trapezoidIntegral( Range inputXY, Callable callable ){
  auto xVec = inputXY | ranges::view::keys;
  auto yVec = inputXY | ranges::view::values;
  auto binWidths = xVec | ranges::view::sliding(2) 
                        | ranges::view::transform([](auto pair){ 
                            return pair[1]-pair[0]; } );
  auto argument = inputXY | ranges::view::transform(callable);
  auto outputWindows = argument | ranges::view::sliding(2);
  auto trapezoid = [](auto binWidth, auto leftRightPair){ 
    return (leftRightPair[0]+leftRightPair[1])*0.5*binWidth;
  };
  auto integral = ranges::view::zip_with(trapezoid,binWidths,outputWindows);
  return ranges::accumulate(integral,0.0);
}
