#include <iostream> 
#include <vector>
#include <cmath>

#include "contin/contin.h"
#include "trans/trans.h"

int main(){
  // Card1
    int nout = 20; 
  // Card2
    // Title goes here
  // Card3
    int ntempr = 1, iprint = 1, nphon = 6;
  // Card4
    int mat = 26; double za = 126.0;
  // Card5
    double awr = 8.93478, spr = 6.15; int npr = 1, iel = 2;
  // Card6
    int nss = 0; double aws = 0.0;
  // Card7
    int nalpha = 5, nbeta = 5, lat = 3;
  // Card8
    std::vector<double> alpha { 0.10, 0.20, 0.40, 0.80, 1.60 };
  // Card9
    std::vector<double> beta  { 0.00, 0.15, 0.30, 0.60, 1.20 };
  // Card10
   // double temp = 200.0;
    std::vector<double> temp_vec {200.0};
  // Card11
    double delta = 3.8; int ni = 6;
  // Card12
    std::vector<double> rho { 0.002, 0.004, 0.02, 0.04, 0.2, 0.4 };
  // Card13
//    double twt = 0.3, c = 0.1, tbeta = 2.0;
    double trans_weight = 0.03, diffusion_const = 1.0, tbeta = 2.0;
  // Card14
    int nd = 0;
  
  double bk = 8.617385e-5;
  double therm = 0.0253;

  // Loop over scatterers and temperatures
  int isecs = 0;
  double arat  = 1.0; // This is for scaling alpha and beta values later
  bool done = false;

  std::vector<std::vector<std::vector<double>>> sym_sab( alpha.size(),
      std::vector<std::vector<double>> (beta.size(), 
      std::vector<double> ( ntempr, 0.0 ) ) );



  while ( not done ){
      if ( isecs == 0 ){ std::cout << "Principal scatterer" << std::endl; }
      if ( isecs >  0 ){ 
        std::cout << "Secondary scatterer" << std::endl;
        arat = aws / awr; 
      }
      
      for ( int itemp = 0; itemp < temp_vec.size(); ++itemp ){ 
        double temp = temp_vec[itemp];
        double tev = bk * temp;
        double sc = 1.0;
        if ( lat == 1 ){ sc = therm/tev; }

        if ( itemp == 1 or temp >= 0 ){
          std::cout << "we want to read in tempdependent parameters" << std::endl;
        } // if 1st temp or some positive temp, we want to calculate
          // the temperature dependent parameters for this specifically 

        // Continuous part of the distribution
        double lambda_s = contin( sym_sab, ntempr, nphon, alpha, beta, lat, delta, rho, tbeta, arat, tev, sc, itemp );
       
        std::cout << sym_sab[2][2][0] << std::endl;
        // Translational part of distribution, if any
        trans( alpha, beta, lat, trans_weight, delta, diffusion_const, sc, arat, 
               tev, sym_sab, itemp, lambda_s );

        std::cout << sym_sab[2][2][0] << std::endl;
      }
      done = true;
  }


    return 0;
}
