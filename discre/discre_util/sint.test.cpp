#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "sint.h"

void equal( double a, double b ){
  if (b == 0.0){ 
    REQUIRE( b-a < 1e-6 );
    return;
  }
  REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 );
}



TEST_CASE( "sint" ){
  GIVEN( "inputs" ){
    double x = -0.1, alpha = 0.1, wt = 2.3, tbart = 407.4545311;
    int b = 4, nbx = 10;
    std::vector<double> bex {-1.2, -0.6, -0.3, -0.15, -0.1, 0.1, 0.15, 0.3, 
      0.6, 1.2, 0.0};
    std::vector<double> rdbex  {1.66666667, 3.33333333, 6.666666667, 20.0, 5.0,
      20.0, 6.66666667, 3.33333333, 1.66666667, 0.0, 0.0};
    std::vector<double> betan {0.1, 0.15, 0.3, 0.6, 1.2};
    std::vector<double> sex  {5.0, 4.0, 3.0, 2.0, 1.0, 1.0, 1.7214159, 
      2.2224546, 2.1952465, 1.5059710, 0.0};
  
    double sintOut;
   // std::cout << "TEST SINT" << std::endl;
   // auto sintOut = sint( x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx );
   // equal( sintOut, 1.0 );


    WHEN( "OPTION A" ){
      THEN( "correct values are output" ){
        x = -1.201; alpha = -0.1;
    std::cout << "TEST A" << std::endl;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        equal( sintOut, 0.0 );
        x = 1.201; alpha = -1e-5;
    std::cout << "TEST A" << std::endl;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        equal( sintOut, 0.0 );
      } // THEN
    } // WHEN

    WHEN( "OPTION B with B2" ){
      THEN( "correct values are output" ){
        x = 1.201; alpha = 0.1;
    std::cout << "TEST B with B2" << std::endl;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        equal( sintOut, 2.548608E-4 );

      } // THEN
    } // WHEN


    WHEN( "OPTION B no B2" ){
      THEN( "correct values are output" ){
        x = -1.201; alpha = 1e-5;
    std::cout << "TEST B no B2" << std::endl;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        equal( sintOut, 1.654202E-16 );
      
      } // THEN
    } // WHEN

    WHEN( "OPTION C" ){
      THEN( "correct values are output" ){
        x = -0.1;
    std::cout << "TEST C" << std::endl;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        equal( sintOut, 1.0 );
        x = -0.15;
    std::cout << "TEST C" << std::endl;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        equal( sintOut, 2.0 );

        x = -0.3;
    std::cout << "TEST C" << std::endl;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        equal( sintOut, 3.0 );

      } // THEN
    } // WHEN
    WHEN( "alpha is changed in some non option A or B situationn" ){
      THEN( "no change in output" ){ 
        std::cout << "TEST CHANGE ALPHA" << std::endl;
        x = -0.55; alpha = 0.1;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        equal( sintOut, 3.812737215 );
        std::cout << "TEST CHANGE ALPHA" << std::endl;
        alpha = 1.0e-5;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        equal( sintOut, 3.812737215 );
      } // THEN
    } // WHEN


    WHEN( "OPTION E and G and I and J" ){
      THEN( "correct values are output" ){
        x = -0.35;
        std::cout << "TEST E and G and I and J" << std::endl;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
        equal( sintOut, 3.14734517 );

        x = -0.65; 
      
        std::cout << "TEST E and G and I and J" << std::endl;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
  //  std::cout << "TEST E" << std::endl;
      //  sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
      //  equal( sintOut, 2.0 );

     //   x = -0.3;
   // std::cout << "TEST E" << std::endl;
      //  sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);
      //  equal( sintOut, 3.0 );

      } // THEN
    } // WHEN
     
    WHEN( "OPTION D and G and I and J" ){
      THEN( "correct values are output" ){
        x = -0.11; alpha = 1.0e-5; 
    std::cout << "TEST D, G, I and J" << std::endl;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);

        equal( sintOut, 1.148698345 );
      } // THEN
    } // WHEN

    WHEN( "OPTION D and F and H and J" ){
      THEN( "correct values are output" ){
        sex =  {-5.0, -4.0, -3.0, -2.0, -1.0, -1.0, -1.7214159,  -2.2224546, 
          -2.1952465, -1.5059710, 0.0};
        x = -0.11; alpha = 1.0e-5; 
    std::cout << "TEST D, F, H and J" << std::endl;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);

        equal( sintOut, 1.921947E-98 );
      } // THEN
    } // WHEN
    WHEN( "OPTION D and F and H, but no J" ){
      THEN( "value of zero is returned because option j not taken" ){
        sex =  {-5.0, -4.0, -3.0, -2.0, -1.0, -1.0, -1.7214159,  -2.2224546, 
          -2.1952465, -1.5059710, 0.0};

        rdbex = {1.66666667, 3.33333333, 6.666666667, 25.0, 5.0, 25.0, 
          6.66666667, 3.33333333, 1.66666667, 0.0, 0.0};
        x = -0.11; alpha = 1.0e-5; 
    std::cout << "TEST D, F, H" << std::endl;
        sintOut = sint(x, bex, rdbex, sex, betan, b, alpha, wt, tbart, nbx);

        equal( sintOut, 0.0 );
      } // THEN
    } // WHEN


  } // GIVEN
} // TEST CASE
