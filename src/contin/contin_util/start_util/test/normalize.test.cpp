#include "catch.hpp"
#include "contin/contin_util/start_util/normalize.h"
#include "generalTools/print.h"

auto equal = [](auto x, auto y, double tol = 1e-6){return x == Approx(y).epsilon(tol);};

TEST_CASE( "normalize" ){
  std::vector<double> p, correct;
  double continWgt, delta;
  GIVEN( "some vector p, spacing delta, normalizing value continWgt" ){
    p = {0.02, 0.04, 0.06, 0.08, 0.10, 0.12};
    continWgt = 0.8, delta = 0.5;
    std::vector<double> betaGrid(p.size());
    for (size_t i = 0; i < betaGrid.size(); ++i){ betaGrid[i] = i*delta; }

    auto beta_P_zipped = ranges::view::zip(betaGrid,p);
    auto normalizedP = normalize( beta_P_zipped, continWgt );
    correct= {2.62155807E-2, 5.24311614E-2, 7.86467421E-2, 0.104862322, 
              0.1310779, 0.15729348};

    THEN( "each value in vector has as most a 1e-6 percent error" ){
      REQUIRE(ranges::equal(normalizedP, correct, equal));
    } // THEN

    {
      p = {7.802837847e-3, 1.560567569e-2, 2.340851354e-2, 3.121135138e-2, 
        3.901418923e-2, 4.681702708e-2, 5.461986492e-2, 6.242270277e-2, 
        7.022554062e-2, 7.802837847e-2, 8.583121631e-2, 9.363405416e-2};

      delta = 0.2; continWgt = 2.0;
      betaGrid.resize(p.size());
      for (size_t i = 0; i < betaGrid.size(); ++i){ betaGrid[i] = i*delta; }
    
      auto beta_P_zipped = ranges::view::zip(betaGrid,p);
      auto normalizedP = normalize( beta_P_zipped, continWgt );
      correct = {5.2978935E-2, 0.10595787, 0.15893680, 0.21191574, 
        0.26489467, 0.31787361, 0.37085254, 0.42383148, 
        0.47681042, 0.52978935, 0.58276829, 0.63574722};


      THEN( "each value in vector has as most a 1e-6 percent error" ){
        REQUIRE(ranges::equal(normalizedP, correct, equal));
      } // THEN
    }
  } // GIVEN


  GIVEN( "two vectors where one is the first scaled by some constant" ){
    std::vector<double> p1 = {0.02, 0.04, 0.06, 0.08, 0.10, 0.12}, 
                        p2 = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2};
    continWgt = 5.0, delta = 0.8;
    auto betaGrid = ranges::view::iota(0,int(p1.size())) 
             | ranges::view::transform([delta](auto x){ return x * delta;});

    auto normalizedP1 = normalize( ranges::view::zip(betaGrid,p1), continWgt );
    auto normalizedP2 = normalize( ranges::view::zip(betaGrid,p2), continWgt );

    std::vector<double> correct= { 3.095826174E-2, 6.191652348E-2, 
      9.287478752E-2, 0.1238330469, 0.1547913063, 0.1857495750};

    THEN( "outputs should be same, and both within 1e-6 to true value" ){
      REQUIRE(ranges::equal(normalizedP1, correct, equal));
      REQUIRE(ranges::equal(normalizedP2, correct, equal));
    } // THEN
  } // GIVEN
} // TEST CASE

