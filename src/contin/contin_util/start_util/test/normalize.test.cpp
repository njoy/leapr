#include "catch.hpp"
#include "contin/contin_util/start_util/normalize.h"

TEST_CASE( "normalize" ){
  double tbeta, delta;
  GIVEN( "a ranges view as input" ){
    auto p1 = ranges::view::iota(1,7) 
            | ranges::view::transform([](auto x){ return 0.02*x; }); 

    WHEN( "vector is being scaled" ){

      tbeta = 0.8, delta = 0.5;

      std::vector<double> correct= {2.62155807E-2, 5.24311614E-2, 7.86467421E-2,
        0.104862322, 0.1310779, 0.15729348};
      auto p3 = ranges::view::zip( normalize( p1, delta, tbeta ), correct );

      THEN( "each value in vector has as most a 1e-6 percent error" ){
        //REQUIRE(ranges::equal(p2,correct));
        RANGES_FOR( const auto& t, p3 ){ 
          REQUIRE( std::get<0>(t) == Approx(std::get<1>(t)).epsilon(1e-6) );
        }
      } // THEN
    } // WHEN

    WHEN( "two vectors are scaled versions of each other" ){
      auto p2 = p1 | ranges::view::transform([](auto x){ return 10*x; });

      tbeta = 5.0, delta = 0.8;

      std::vector<double> correct= {3.095826174E-2, 6.191652348E-2, 
        9.287478752E-2, 0.1238330469, 0.1547913063, 0.1857495750};

      auto p3 = ranges::view::zip(normalize( p1, delta, tbeta ),
                                  normalize( p2, delta, tbeta ),correct);

      THEN( "outputs should be same, and both within 1e-6 to true value" ){
        //REQUIRE(ranges::equal(p2,correct));
        RANGES_FOR( const auto& t, p3 ){ 
          REQUIRE( std::get<0>(t) == Approx(std::get<2>(t)).epsilon(1e-6) );
          REQUIRE( std::get<1>(t) == Approx(std::get<2>(t)).epsilon(1e-6) );
        }
      } // THEN
    } // WHEN
  } // GIVEN

  GIVEN( "a vector as input" ){
    std::vector<double> p;

    WHEN( "vector is being scaled" ){

      p = {0.02, 0.04, 0.06, 0.08, 0.10, 0.12};
      tbeta = 0.8, delta = 0.5;

      p = normalize( p, delta, tbeta );
      std::vector<double> correct= {2.62155807E-2, 5.24311614E-2, 
        7.86467421E-2, 0.104862322, 0.1310779, 0.15729348};

      THEN( "each value in vector has as most a 1e-6 percent error" ){
        REQUIRE( p.size() == correct.size() );
        for ( size_t i = 0; i < p.size(); ++i ){
          REQUIRE( p[i] == Approx( correct[i] ).epsilon(1e-6 ) );  
        }
      } // THEN

      p = {7.802837847e-3, 1.560567569e-2, 2.340851354e-2, 3.121135138e-2, 
           3.901418923e-2, 4.681702708e-2, 5.461986492e-2, 6.242270277e-2, 
           7.022554062e-2, 7.802837847e-2, 8.583121631e-2, 9.363405416e-2};

      delta = 0.2; tbeta = 2.0;
      p = normalize( p, delta, tbeta );
      correct = {5.2978935E-2, 0.10595787, 0.15893680, 0.21191574, 
        0.26489467, 0.31787361, 0.37085254, 0.42383148, 
        0.47681042, 0.52978935, 0.58276829, 0.63574722};

      THEN( "each value in vector has as most a 1e-6 percent error" ){
        REQUIRE( p.size() == correct.size() );
        for ( size_t i = 0; i < p.size(); ++i ){
          REQUIRE( p[i] == Approx( correct[i] ).epsilon(1e-6 ) );  
        }
      } // THEN
    } // WHEn

    WHEN( "two vectors are scaled versions of each other" ){
      std::vector<double> p1 = {0.02, 0.04, 0.06, 0.08, 0.10, 0.12}, 
                          p2 = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2};
      tbeta = 5.0, delta = 0.8;

      p1 = normalize( p1, delta, tbeta );
      p2 = normalize( p2, delta, tbeta );

      std::vector<double> correct= {3.095826174E-2, 6.191652348E-2, 
        9.287478752E-2, 0.1238330469, 0.1547913063, 0.1857495750};

      THEN( "outputs should be same, and both within 1e-6 to true value" ){
        REQUIRE( p1.size() == correct.size() );
        REQUIRE( p2.size() == correct.size() );
        for ( size_t i = 0; i < p.size(); ++i ){
          REQUIRE( p1[i] == Approx( correct[i] ).epsilon(1e-6 ) );  
          REQUIRE( p2[i] == Approx( correct[i] ).epsilon(1e-6 ) );  
        }
      } // THEN
    } // WHEN
  } // GIVEN

} // TEST CASE

