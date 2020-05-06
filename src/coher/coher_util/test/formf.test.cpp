#include <iostream>
#include "catch.hpp"
#include "coher/coher_util/formf.h"


TEST_CASE( "Form factor |F|^2 generation for coherent elastic scattering" ){
  using std::get;
  std::vector<double> outputs;
  std::vector<std::tuple<int,int,int,int>> inputs;
  std::tuple<int,int,int,int> in;

  GIVEN( "graphite input" ){
    THEN( "for various l1,l2,l3 combinations, correct form factor is returned" ){
      inputs = {{1,0,0,0}, {1,1,0,0}, {1,0,1,0}, {1,0,0,1}, {1,1,1,0}, {1,1,0,1}, 
                {1,0,1,1}, {1,1,1,1}, {1,1,2,3}, {1,1,0,3}, {1,1,4,0}, {1,0,4,8}};
      outputs = {4, 0.25, 0.25, 0, 4, 0.75, 0.75, 0, 0.75, 0.75, 4, 0.25};
      for (size_t i = 0; i < inputs.size(); ++i ){
        in = inputs[i];
        REQUIRE(outputs[i] == Approx(formf(get<0>(in),get<1>(in),
                                           get<2>(in),get<3>(in))).epsilon(1e-6));
      }
    } // THEN
  } // GIVEN
  GIVEN( "beryllium metal input" ){
    THEN( "for various l1,l2,l3 combinations, correct form factor is returned" ){
      inputs = {{2,0,0,0}, {2,1,0,0}, {2,0,1,0}, {2,0,0,1}, {2,1,1,0}, {2,1,0,1}, 
                {2,0,1,1}, {2,1,1,1}, {2,1,2,3}, {2,1,0,3}, {2,1,4,0}, {2,0,4,8}};
      outputs = {2, 0.5, 0.5, 0, 2, 1.5, 1.5, 0, 1.5, 1.5, 2, 0.5};
      for (size_t i = 0; i < inputs.size(); ++i ){
        in = inputs[i];
        REQUIRE(outputs[i] == Approx(formf(get<0>(in),get<1>(in),
                                           get<2>(in),get<3>(in))).epsilon(1e-6));
      }
    } // THEN
  } // GIVEN
  GIVEN( "beryllium oxide input" ){
    THEN( "for various l1,l2,l3 combinations, correct form factor is returned" ){
      inputs = {{3,0,0,0}, {3,1,0,0}, {3,0,1,0}, {3,0,0,1}, {3,1,1,0}, {3,1,0,1}, 
                {3,0,1,1}, {3,1,1,1}, {3,1,2,3}, {3,1,0,3}, {3,1,4,0}, {3,0,4,8}};
      outputs = {46.18, 11.545, 11.545, 0, 46.18, 5.67393345, 5.67393345, 0, 
                 29.6660665, 29.6660665, 46.18, 11.545 };
      for (size_t i = 0; i < inputs.size(); ++i ){
        in = inputs[i];
        REQUIRE(outputs[i] == Approx(formf(get<0>(in),get<1>(in),
                                           get<2>(in),get<3>(in))).epsilon(1e-6));
      }
    } // THEN 
  } // GIVEN
  GIVEN( "fcc (like aluminum or lead) lattice input" ){
    THEN( "for various l1,l2,l3 combinations, correct form factor is returned" ){
      inputs = {{4,0,0,0}, {4,1,0,0}, {4,0,1,0}, {4,0,0,1}, {4,1,1,0}, {4,1,0,1}, 
                {4,0,1,1}, {4,1,1,1}, {4,1,2,3}, {4,1,0,3}, {4,1,4,0}, {4,0,4,8}, 
                {5,0,0,0}, {5,1,0,0}, {5,0,1,0}, {5,0,0,1}, {5,1,1,0}, {5,1,0,1}, 
                {5,0,1,1}, {5,1,1,1}, {5,1,2,3}, {5,1,0,3}, {5,1,4,0}, {5,0,4,8}};
      for (auto in : inputs){
        REQUIRE( 16 == Approx(formf(get<0>(in),get<1>(in),
                                    get<2>(in),get<3>(in))).epsilon(1e-6) );
      }
    } // THEN
  } // GIVEN
  GIVEN( "bcc (iron) lattice input" ){
    THEN( "for various l1,l2,l3 combinations, correct form factor is returned" ){
      inputs = {{6,0,0,0}, {6,1,0,0}, {6,0,1,0}, {6,0,0,1}, {6,1,1,0}, {6,1,0,1},
                {6,0,1,1}, {6,1,1,1}, {6,1,2,3}, {6,1,0,3}, {6,1,4,0}, {6,0,4,8}};
      for (auto in : inputs){
        REQUIRE( 4 == Approx(formf(get<0>(in),get<1>(in),
                                   get<2>(in),get<3>(in))).epsilon(1e-6) );
     }
   } // THEN
  } // GIVEN
} // TEST CASE
