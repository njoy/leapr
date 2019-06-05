#ifndef LEAPR_CONTIN_START_FSUM_HH
#define LEAPR_CONTIN_START_FSUM_HH
#include <range/v3/all.hpp>
#include <tuple>
#include <iostream>

template <typename F, typename A> 
double fsum( const int& n, const A& p, const F& tau, const A& betaGrid ){
  /* Inputs
   * ------------------------------------------------------------------------
   * n       : appears in equation being evaluated
   * p       : P(beta), which is defined in Eq. 507
   * tau     : appears in equation being evaluated
   * delta_b : increment size used for integral. Doing a Riemann sum with
   *           rectangle sum delta_b wide
   *
   * Operations
   * ------------------------------------------------------------------------
   * * Compute an integral of the form
   *
   *   int 0->infty  2 * P(beta) * beta^n * sinh( tau*beta ) dbeta
   *                                or
   *   int 0->infty  2 * P(beta) * beta^n * cosh( tau*beta ) dbeta
   *
   *  depending on whether n is odd of even, respectively.
   *
   * Outputs
   * ------------------------------------------------------------------------
   * * The computed integral is returned
   *
   */

  double func_sum = 0, func_val = 0, dx_left, dx_right;
  bool even = ( 1 - 2*(n%2) == 1 ); // +1 if even, -1 if odd. This is to
                                    // help differentiate betwene sinh and
                                    // cosh while evaluating the integrand
  for( size_t i = 0; i < p.size(); ++i ){
      func_val = even ? 2*p[i]*cosh(betaGrid[i]*tau)*std::pow(betaGrid[i],n) :
                        2*p[i]*sinh(betaGrid[i]*tau)*std::pow(betaGrid[i],n) ;
      dx_left  = 0.0; // We look a half space left, if there's anything <---
      dx_right = 0.0; // We look a half space right, if there's anything --->
      if ( i != 0 )         { dx_left  = (betaGrid[i]-betaGrid[i-1])*0.5; }
      if ( i != p.size()-1 ){ dx_right = (betaGrid[i+1]-betaGrid[i])*0.5; }
      func_val = func_val * ( dx_left + dx_right ); 
    func_sum += func_val;
  } // for i in p
  return func_sum;
}


template <typename Float>
auto even(Float tau, int n, Float beta){
  using std::cosh; using std::pow;
  return 2.0*cosh(beta*tau)*pow(beta,n);
}

template <typename Float>
auto odd(Float tau, int n, Float beta){
  using std::sinh; using std::pow;
  return 2.0*sinh(beta*tau)*pow(beta,n);
}

  //ranges::for_each(xVec, [](auto x){
  //      std::cout << x << ' ';
  //  });


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


template <typename Float, typename Tuple>
auto even2(Float tau, int n, Tuple xyPair){
  using std::cosh; using std::pow;
  Float beta = std::get<0>(xyPair);
  Float pVal = std::get<1>(xyPair);
  return 2.0*pVal*cosh(beta*tau)*pow(beta,n);
}

template <typename Float, typename Tuple>
auto odd2(Float tau, int n, Tuple xyPair){
  using std::sinh; using std::pow;
  Float beta = std::get<0>(xyPair);
  Float pVal = std::get<1>(xyPair);
  return 2.0*pVal*sinh(beta*tau)*pow(beta,n);
}


template <typename Range, typename Float>// = Range::value_type >
auto fsum3( int n, Range p, Float tau, Range betas ){
  auto inputXY = ranges::view::zip(betas,p);
  if (n%2 == 0){ 
    auto evenLambda = [tau, n](auto xyPair){ return even2(tau,n,xyPair); };
    return trapezoidIntegral( inputXY, evenLambda );
  }
  else { 
    auto oddLambda  = [tau, n](auto xyPair){ return odd2(tau,n,xyPair);  };
    return trapezoidIntegral( inputXY, oddLambda );
  }

 
}





template <typename Range, typename Float>// = Range::value_type >
auto fsum2( int n, Range p, Float tau, Range betas ){

  auto do_the_thing = [&](auto lambda){  
    auto binWidths = betas | ranges::view::sliding(2) 
                           | ranges::view::transform([](auto pair){ 
                               return pair[1]-pair[0]; } );
    auto argument = ranges::view::zip_with([](auto a, auto b){return a*b;},
                      betas | ranges::view::transform(lambda), p );
    auto outputWindows = argument | ranges::view::sliding(2);
    auto trapezoid = [](auto binWidth, auto leftRightPair){ 
      return (leftRightPair[0]+leftRightPair[1])*0.5*binWidth;
    };
    auto integral = ranges::view::zip_with(trapezoid,binWidths,outputWindows);
    auto XY = ranges::view::zip(betas,argument);
    return ranges::accumulate(integral,0.0);
  };

  if (n%2 == 0){ 
    auto evenLambda = [tau, n](auto beta){ return even(tau,n,beta); };
    return do_the_thing(evenLambda);
  }
  else { 
    auto oddLambda  = [tau, n](auto beta){ return odd(tau,n,beta);  };
    return do_the_thing(oddLambda);
  }

}










/*
template <typename Callable, typename Range,
  typename Inputs = Range::value_type, 
  typename Result = std::invoke_result_t<Callable,Inputs> >
auto fsum_new( Callable callable, Range grid ){
  std::vector<Result> result;
  result.reserve(grid.size());
  for (auto&& gridVals : grid){
    result.push_back(callable(gridVals));
  }
  return result;
}
*/











#endif
