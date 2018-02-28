//#include <Eigen/Dense>
#include <iostream> 
#include <vector>
#include <cmath>
#include <tuple>
#include <string>
#include "contin/contin.h"
#include "trans/trans.h"
#include "discre/discre.h"
#include "coldh/coldh.h"
#include <time.h>


int leapr( int nout, std::string title, int ntempr, int iprint, int nphon, 
  int mat, double za, double awr, double spr, int npr, int iel, int ncold,
  int nss, double aws, int nalpha, int nbeta, int lat, 
  std::vector<double> alpha, std::vector<double> beta, 
  std::vector<double> temp_vec, double delta, int ni, std::vector<double> rho,
  double trans_weight, double diffusion_const, double tbeta, int nd, 
  std::vector<double> oscEnergies, std::vector<double> oscWeights, int nka, 
  double dka, std::vector<double> kappaVals ){
  /*
  // Card1
    int nout = 20; 
  // Card2
    // Title goes here
  // Card3
    int ntempr = 1, iprint = 1, nphon = 3;
  // Card4
    int mat = 26; double za = 126.0;
  // Card5
    double awr = 8.93478, spr = 6.15; int npr = 1, iel = 2; int ncold = 0;
  // Card6
    int nss = 0; double aws = 0.0;
  // Card7
    int nalpha = 5, nbeta = 5, lat = 3;
  // Card8
    std::vector<double> alpha { 0.10, 0.20, 0.40, 0.80, 1.60 };
  // Card9
    std::vector<double> beta  { 0.10, 0.15, 0.30, 0.60, 1.20 };
  // Card10
   // double temp = 200.0;
    std::vector<double> temp_vec {200.0};
  // Card11
    double delta = 3.8; int ni = 6;
  // Card12
    std::vector<double> rho { 0.002, 0.004, 0.02, 0.04, 0.2, 0.4 };
  // Card13
//    double twt = 0.3, c = 0.1, tbeta = 2.0;
    double trans_weight = 0.3, diffusion_const = 1.0, tbeta = 2.0;
  // Card14
    int nd = 2;
  // Card15
    std::vector<double> oscEnergies {1.0, 2.0};
  // Card16
    std::vector<double> oscWeights {0.3, 0.4};
  // Card17 
    int nka = 4; double dka = 0.01;
    std::vector<double> kappaVals { 0.1, 0.2, 0.4, 0.7 };
    */
  clock_t t;
  t = clock();


    //Eigen::MatrixXd matrix1 = Eigen::MatrixXd::Random(2,3);
    //std::cout << matrix1 << std::endl;
  double bk = 8.617385e-5;
  double therm = 0.0253;

  // Loop over scatterers and temperatures
  int isecs = 0;
  double arat  = 1.0; // This is for scaling alpha and beta values later
  bool done = false;

  std::vector<std::vector<std::vector<double>>> sym_sab( alpha.size(),
    std::vector<std::vector<double>> (beta.size(), 
    std::vector<double> ( ntempr, 0.0 ) ) );

  // This is only going to be used by coldh
  std::vector<std::vector<std::vector<double>>> sym_sab_2( alpha.size(),
    std::vector<std::vector<double>> (beta.size(), 
    std::vector<double> ( ntempr, 0.0 ) ) );



  std::vector<double> t_eff_vec ( temp_vec.size(), 0.0 );

  while ( not done ){
    if ( isecs == 0 ){ std::cout << "Principal scatterer" << std::endl; }
    if ( isecs >  0 ){ 
      std::cout << "Secondary scatterer" << std::endl;
      arat = aws / awr; 
    }
      
    for ( size_t itemp = 0; itemp < temp_vec.size(); ++itemp ){ 
      double temp = temp_vec[itemp];
      double tev = bk * temp;
      double sc = 1.0;
      if ( lat == 1 ){ sc = therm/tev; }
      double scaling = sc/arat;
      if ( itemp == 1 or temp >= 0 ){
        std::cout << "we want to read in tempdependent parameters" << std::endl;
      } // if 1st temp or some positive temp, we want to calculate
        // the temperature dependent parameters for this specifically 

      /*
      std::cout << sym_sab[0][0][0] << std::endl;
      std::cout << sym_sab[1][1][0] << std::endl;
      std::cout << sym_sab[2][2][0] << std::endl;
      std::cout << sym_sab[3][3][0] << std::endl;
      std::cout << sym_sab[4][4][0] << std::endl;
      std::cout << "    " << std::endl;
      */


      // Continuous part of the distribution
      //std::cout << "\n-------- contin" << std::endl;

      auto lambda_s_t_eff = contin( itemp, nphon, delta, tbeta, scaling, tev,
        sc, rho, alpha, beta, sym_sab );
      double lambda_s = std::get<0>(lambda_s_t_eff);

      double t_contin = clock();

     // update the effective temperature list
      t_eff_vec[itemp] = std::get<1>(lambda_s_t_eff) * temp;

 
      // Translational part of distribution, if any
      //std::cout << "\n-------- trans" << std::endl;
      if ( trans_weight > 0.0 ){
        trans( alpha, beta, trans_weight, delta, diffusion_const, sc, scaling,
          itemp, lambda_s, tbeta, t_eff_vec, temp_vec, sym_sab );
      }
      double t_trans = clock();


      if ( oscEnergies.size() > 0 ){
      //std::cout << "\n-------- discre" << std::endl;
        discre( itemp, sc, scaling, tev, lambda_s, trans_weight, tbeta, alpha,
          beta, temp_vec, oscEnergies, oscWeights, t_eff_vec, sym_sab );
      }
      double t_discre = clock();

      return 0;
      std::cout << "contin: " << ((float)(t_contin-t))/CLOCKS_PER_SEC 
        << " seconds" << std::endl;
      std::cout << "trans:  " << ((float)(t_trans-t_contin))/CLOCKS_PER_SEC 
        << " seconds" << std::endl;
      std::cout << "discre: " << ((float)(t_discre-t_trans))/CLOCKS_PER_SEC 
        << " seconds" << std::endl;
      std::cout << "\n-------------------\nTotal: " << ((float)clock()-t)/CLOCKS_PER_SEC 
        << " seconds" << std::endl;

      if ( ncold != 0 ){
	bool free = false;
        coldh( itemp, temp, tev, ncold, trans_weight, tbeta, t_eff_vec, 
	  scaling, alpha, beta, dka, kappaVals, nbeta, lat, free, sym_sab, 
	  sym_sab_2 );
      }

    }
    done = true;
  }

    return 0;
    std::cout << "Card1: " << nout << std::endl;
    std::cout << "Card2: " <<  title<< std::endl;
    std::cout << "Card3: " << ntempr <<  "     " << iprint << "     " << nphon<< std::endl;
    std::cout << "Card4: " << mat    <<  "     " << za << std::endl;
    std::cout << "Card5: " << awr    <<  "     " << spr    << "     " << npr << "      " << iel << "     " << ncold << std::endl;
    std::cout << "Card6: " << nss    <<  "     " << aws << std::endl;
    std::cout << "Card7: " << nalpha <<  "     " << nbeta  << "     " << lat<<   std::endl;
    std::cout << "Card8: " << std::endl;
    for( auto entry : alpha   ){ std::cout << " -- " << entry << std::endl; }
    std::cout << "Card9: " << std::endl;
    for( auto entry : beta    ){ std::cout << " -- " << entry << std::endl; }
    std::cout << "Card10: " << std::endl;
    for( auto entry : temp_vec){ std::cout << " -- " << entry << std::endl; }
    std::cout << "Card11: " << delta << "     " << ni <<  std::endl;
    std::cout << "Card12: " << std::endl;
    for( auto entry : rho     ){ std::cout << " -- " << entry << std::endl; }
    std::cout << "Card13: " << trans_weight << "     " << diffusion_const << "      " << tbeta << std::endl;
    std::cout << "Card14: " << nd << std::endl;
    std::cout << "Card15: " << std::endl;
    for( auto entry : oscEnergies){ std::cout << " -- " << entry << std::endl; }
    std::cout << "Card16: " << std::endl;
    for( auto entry : oscWeights){ std::cout << " -- " << entry << std::endl; }
    std::cout << "Card17: " << nka << "     " << dka << std::endl;
    std::cout << "Card18: " << std::endl;
    for( auto entry : kappaVals){ std::cout << " -- " << entry << std::endl; }

    std::cout << "\n\n\n " << std::endl;


}

