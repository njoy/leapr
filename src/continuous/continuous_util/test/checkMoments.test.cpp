#include "catch.hpp"
#include "continuous/continuous_util/checkMoments.h"
#include <range/v3/all.hpp>

auto equal = [](auto x, auto y, double tol = 1e-6){ return x == Approx(y).epsilon(tol); };

TEST_CASE( "check moments" ){

  int ntempr = 1;
  std::vector<double> alpha { 1.008e-2, 1.5e-2, 2.52e-2, 3.3e-2, 5.0406e-2 },
  beta { 0.0, 6.375e-3, 1.275e-2, 2.55e-2, 3.825e-2, 5.1e-2, 6.575e-2 };

  std::vector<double> ssm(alpha.size()*beta.size(),0.0);

  GIVEN( "inputs" ){

    std::vector<int> maxt ( 1000, 0.0 );
    maxt[0] = 6;
    double sc = 0.9918667, f0 = 0.23520651, tbeta = 0.444444, arat = 1, 
           tbar = 1.93448465, explim = -250;

    for ( auto& a : alpha ){ a *= sc/arat; }
    for ( auto& b : beta  ){ b *= sc; }

    checkMoments( alpha, beta, maxt, f0, tbeta, tbar, ssm );

    std::vector<double> correctSab {0.00000000, 3.04230221, 3.03666664, 
      3.00439217, 2.94493755, 2.85993079, 2.73274581, 0.00000000, 2.49419786, 
      2.49242782, 2.47724974, 2.44682038, 2.40170401, 2.33228656, 0.00000000, 
      1.92380592, 1.92426566, 1.91982023, 1.90827115, 1.88974668, 1.85988639, 
      0.00000000, 1.68058049, 1.68153661, 1.67986667, 1.67343655, 1.66230091, 
      1.64367368, 0.00000000, 1.35861931, 1.35989255, 1.36054304, 1.35866398, 
      1.35426586, 1.34606790 };

    ranges::equal(correctSab,ssm,equal);
    
  } // GIVEN

  GIVEN( "other inputs" ){

    std::vector<int> maxt ( 1000, 0.0 );
    maxt[0] = 6; maxt[1] = 1; maxt[2] = 2; maxt[3] = 1; maxt[4] = 2;
    maxt[5] = 1; maxt[6] = 2;

    double sc = 0.9918667, f0 = 0.23520651, tbeta = 0.444444, arat = 1, 
           tbar = 1.93448465, explim = -250;

    for ( auto& a : alpha ){ a *= sc/arat; }
    for ( auto& b : beta  ){ b *= sc; }


    checkMoments( alpha, beta, maxt, f0, tbeta, tbar, ssm );

    std::vector<double> correctSab { 0.00000000, 3.04230221, 0.00000000, 
      3.00439217, 0.00000000, 2.85993079, 0.00000000, 0.00000000, 2.49419786, 
      2.49242782, 2.47724974, 2.44682038, 2.40170401, 2.33228656, 0.00000000, 
      1.92380592, 1.92426566, 1.91982023, 1.90827115, 1.88974668, 1.85988639, 
      0.00000000, 1.68058049, 1.68153661, 1.67986667, 1.67343655, 1.66230091, 
      1.64367368, 0.00000000, 1.35861931, 1.35989255, 1.36054304, 1.35866398, 
      1.35426586, 1.34606790 };

    ranges::equal(correctSab,ssm,equal);
    
  } // GIVEN

} // TEST_CASE
