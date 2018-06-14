#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "coh/coh_util/sigcoh_util/form.h"


TEST_CASE( "form" ){
  int lat, l1, l2, l3;
  GIVEN( "the material of interest is graphite" ){
    lat = 1;
    THEN( "computed form factors are correct" ){
      std::vector<std::tuple<int,int,int>> inputs { { 0, 0, 0 }, { 1, 0, 0 }, 
        { 0, 1, 0 }, { 0, 0, 1 }, { 1, 1, 0 }, { 1, 0, 1 }, { 0, 1, 1 }, 
        { 1, 1, 1 }, {-1, 0, 0 }, { 0,-1, 0 }, { 0, 0,-1 }, {-1, 1, 0 },
        { 1, 0,-1 }, { 0,-1, 1 }, {-1, 1,-1 }, { 2, 3,-4 }, { 3, 1, 0 }, 
        { 0,-2,-3 }, { 1,-2, 1 } };
      std::vector<double> correctVals { 4, 0.25, 0.25, 0, 4, 0.75, 0.75, 0, 
        0.25, 0.25, 0, 0.25, 0.75, 0.75, 0.75, 0.25, 0.25, 0.75, 0};
      for ( size_t i = 0; i < inputs.size(); ++i ){
        REQUIRE( correctVals[i] == Approx( form( lat, std::get<0>(inputs[i]), 
          std::get<1>(inputs[i]), std::get<2>(inputs[i]) ) ).epsilon(1e-6) );
      }
      
    } // THEN
  } // GIVEN
  GIVEN( "the material of interest is beryllium" ){
    lat = 2;
    THEN( "computed form factors are correct" ){
      std::vector<std::tuple<int,int,int>> inputs { { 0, 0, 0 }, { 1, 0, 0 }, 
        { 0, 1, 0 }, { 0, 0, 1 }, { 1, 1, 0 }, { 1, 0, 1 }, { 0, 1, 1 }, 
        { 1, 1, 1 }, {-1, 0, 0 }, { 0,-1, 0 }, { 0, 0,-1 }, {-1, 1, 0 },
        { 1, 0,-1 }, { 0,-1, 1 }, {-1, 1,-1 }, { 2, 3,-4 }, { 3, 1, 0 }, 
        { 0,-2,-3 }, { 1,-2, 1 } };
      std::vector<double> correctVals { 
         2, 0.5, 0.5, 0, 2, 1.5, 1.5, 0, 0.5, 0.5, 0, 0.5, 1.5, 1.5, 1.5, 0.5, 0.5, 1.5, 0 };

      for ( size_t i = 0; i < inputs.size(); ++i ){
        REQUIRE( correctVals[i] == Approx( form( lat, std::get<0>(inputs[i]), 
          std::get<1>(inputs[i]), std::get<2>(inputs[i]) ) ).epsilon(1e-6) );
      }
    } // THEN
  } // GIVEN
  GIVEN( "the material of interest is beryllium oxide" ){
    lat = 3;
    THEN( "computed form factors are correct" ){
      std::vector<std::tuple<int,int,int>> inputs { { 0, 0, 0 }, { 1, 0, 0 }, 
        { 0, 1, 0 }, { 0, 0, 1 }, { 1, 1, 0 }, { 1, 0, 1 }, { 0, 1, 1 }, 
        { 1, 1, 1 }, {-1, 0, 0 }, { 0,-1, 0 }, { 0, 0,-1 }, {-1, 1, 0 },
        { 1, 0,-1 }, { 0,-1, 1 }, {-1, 1,-1 }, { 2, 3,-4 }, { 3, 1, 0 }, 
        { 0,-2,-3 }, { 1,-2, 1 } };
      std::vector<double> correctVals { 46.18, 11.545, 11.545, 0, 46.18, 5.67393351, 5.67393351, 0, 11.545, 11.545, 0, 11.545, 5.67393351, 5.67393351, 5.67393351, 0.235, 11.545, 29.6660671, 0 };

      for ( size_t i = 0; i < inputs.size(); ++i ){
        REQUIRE( correctVals[i] == Approx( form( lat, std::get<0>(inputs[i]), 
          std::get<1>(inputs[i]), std::get<2>(inputs[i]) ) ).epsilon(1e-6) );
      }
    } // THEN
  } // GIVEN
  GIVEN( "the material of interest is not one of the three preset materials" ){
    std::vector<int> latVals = { 4, 5, 7, 10, 100, -1, -2, -3, 0 };
    THEN( "computed form factors are correct" ){
      std::vector<std::tuple<int,int,int>> inputs { { 0, 0, 0 }, { 1, 0, 0 }, 
        { 0, 1, 0 }, { 0, 0, 1 }, { 1, 1, 0 }, { 1, 0, 1 }, { 0, 1, 1 }, 
        { 1, 1, 1 }, {-1, 0, 0 }, { 0,-1, 0 }, { 0, 0,-1 }, {-1, 1, 0 },
        { 1, 0,-1 }, { 0,-1, 1 }, {-1, 1,-1 }, { 2, 3,-4 }, { 3, 1, 0 }, 
        { 0,-2,-3 }, { 1,-2, 1 } };

      for ( size_t i = 0; i < inputs.size(); ++i ){
        for ( auto lat : latVals ){
          REQUIRE( 0 == Approx( form( lat, std::get<0>(inputs[i]), 
            std::get<1>(inputs[i]), std::get<2>(inputs[i]) ) ).epsilon(1e-6) );
        }
      }
    } // THEN
  } // GIVEN



} // TEST CASE
