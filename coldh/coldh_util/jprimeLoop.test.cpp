#define CATCH_CONFIG_MAIN
#include <iostream>
#include "../../catch.hpp"
#include "jprimeLoop.h"


TEST_CASE( "jprime loop for odd values" ){
  GIVEN( "inputs" ){
    double total = 0.02;
    double alp = 0.23;
    double wt = 2.3;
    double be = -1.2; 
    double al = 0.1;
    double x = 0.85, swo = 0.87, pj = 0.48;
    double tbart = 950.0;
    int j = 1, ifree = 0, jj = 0, maxbb = 11, nbeta = 5;
    std::vector<double> bex { -1.2, -0.6, -0.3, -0.15, -0.1, 0.1, 0.15, 0.3, 
      0.6, 1.2, 0.0 };
    std::vector<double> rdbex { 1.6666666, 3.3333333, 6.6666666, 20.0, 5.0, 
      20.0, 6.66666, 3.333333, 1.666666, 0.0, 0.0 };
    std::vector<double> sex { 5.0, 4.0, 3.0, 2.0, 1.0, 1.0, 1.721415952, 
      2.2224546, 2.1952465, 1.5059710, 5.0 };
    std::vector<double> betan { 0.1, 0.15, 0.3, 0.6, 1.2 };

    double snlk = jPrimeOdd( total, alp, j, ifree, be, x, swo, pj, jj, bex, 
        rdbex, sex, betan, al, wt, tbart, maxbb, nbeta );
    std::cout << snlk << std::endl;
    
    


  } // GIVEN
} // TEST CASE
