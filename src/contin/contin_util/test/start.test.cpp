#include "catch.hpp"
#include "contin/contin_util/start.h"
#include <vector>
#include <tuple>
#include <memory>
#include <iostream>

template <typename A, typename B, typename C> 
void check( const A& correctT1,     const B& outputTuple, 
            const C& correctLambda, const C& correctEffTemp ){
  REQUIRE( (std::get<2>(outputTuple).size() == correctT1.size()) );

  RANGES_FOR( auto t, ranges::view::zip(std::get<2>(outputTuple),correctT1) ){
    REQUIRE( (std::get<1>(t) == Approx(std::get<0>(t)).epsilon(1e-6)) );
  }
  REQUIRE( (std::get<0>(outputTuple) == Approx(correctLambda).epsilon(1e-6))  );
  REQUIRE( (std::get<1>(outputTuple) == Approx(correctEffTemp).epsilon(1e-6)) );

}


TEST_CASE( "start function" ){

  std::vector<double> p, correct;
  double delta, tev, tbeta, correctLambda, correctEffectiveTemp;

  GIVEN( "input phonon distribution given as vector" ){
    WHEN( "small phonon distribution is provided" ){

      p = {0.1, 0.2, 0.3, 0.5, 0.8, 1.3};
      delta = 0.0001; tev = 0.001; tbeta = 1.0;
      auto output = start(p, delta, tev, tbeta);

      THEN( "T1, debye waller, and effective temp are correctly computed" ){
        correct = {1.9662112, 2.06616, 0.8135182, 0.632185, 0.5964, 0.649624};
        correctLambda = 41.517752; correctEffectiveTemp = 1.0118507;
        check( correct, output, correctLambda, correctEffectiveTemp );
      } // THEN
    } // WHEN

    WHEN( "slightly larger phonon distribution is provided" ){

      p = {0.001, 0.003, 0.006, 0.008, 0.01, 0.03, 0.06, 0.08, 1.0, 3.0, 6.0, 8.0};
      tbeta = 2.0;

      AND_WHEN( "spacing and temperature are very medium" ){
        delta = 0.015, tev = 1.723477E-2;
        auto output = start(p, delta, tev, tbeta);

        THEN( "T1, debye waller, and effective temp are correctly computed" ){
          correct = {2.691905E-3, 4.031148E-3, 2.841214E-3, 2.247653E-3, 
                     2.014354E-3, 4.746877E-3, 7.851893E-3, 8.945384E-3, 
                     9.771149E-2, 0.2604205, 0.4686490, 0.5680046};
          correctLambda = 0.238142827, correctEffectiveTemp = 4.298571365;
          check( correct, output, correctLambda, correctEffectiveTemp );
        } // THEN
      } // AND WHEN

      AND_WHEN( "temperature is changed" ){
        delta = 0.015; tev = 4.3086925E-2;
        auto output = start(p, delta, tev, tbeta);
 
        THEN( "T1, debye waller, and effective temp are correctly computed" ){
          correct = {1.54443108E-2, 1.82883208E-2, 1.07199712E-2, 7.37431363E-3, 
                     5.96172958E-3, 1.30407541E-2, 2.04553735E-2, 2.24449206E-2, 
                     0.238766209, 0.62462916, 1.10947058, 1.332378};
          correctLambda = 0.6485580; correctEffectiveTemp = 1.836531025;
          check( correct, output, correctLambda, correctEffectiveTemp );

        } // THEN
      } // AND WHEN
  
      AND_WHEN( "spacing is very large" ){
        tev = 4.3086925E-2; delta = 4.015;
        auto output = start(p, delta, tev, tbeta);

        THEN( "T1, debye waller, and effective temp are correctly computed" ){
          correct = {2.36543334e-7, 2.20419880e-5, 2.20419880e-5, 1.95928782e-5, 
                     1.83683233e-5, 4.40839761e-5, 7.34732935e-5, 8.39694783e-5, 
                     9.18416169e-4, 2.44910978e-3, 4.40839761e-3, 5.34351226e-3};
          correctLambda = 2.2081257E-3; correctEffectiveTemp = 459.94244303;
          check( correct, output, correctLambda, correctEffectiveTemp );
        } // THEN
      } // AND WHEN
    } // WHEN
  } // GIVEN

  GIVEN( "input phonon distribution given as ranges view" ){
    WHEN( "phonon distribution, spacing, and temp in ev" ){
      p = {1, 2, 3, 5, 8, 13};
      auto p1 = p | ranges::view::transform([](auto x){ return x * 0.1; });

      delta = 0.0001; tev = 0.001; tbeta = 1.0;
      auto output = start(p1, delta, tev, tbeta);

      THEN( "T1, debye waller, and effective temp are correctly computed" ){
        correct = {1.9662112, 2.06616, 0.8135182, 0.632185, 0.5964, 0.649624};
        correctLambda = 41.517752; correctEffectiveTemp = 1.0118507;
        check( correct, output, correctLambda, correctEffectiveTemp );

      } // THEN
    } // WHEN
  } // GIVEN
} // TEST_CASE
