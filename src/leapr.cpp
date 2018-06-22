#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
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
#include <range/v3/all.hpp>


auto leapr( int nout, std::string title, int ntempr, int iprint, int nphon, 
  int mat, double za, double awr, double spr, int npr, int iel, int ncold,
  int nss, double aws, int nalpha, int nbeta, int lat, 
  std::vector<double> alpha, std::vector<double> beta, 
  std::vector<double> temp_vec, double delta, int ni, std::vector<double> rho,
  double trans_weight, double diffusion_const, double tbeta, int nd, 
  std::vector<double> oscEnergies, std::vector<double> oscWeights, int nka, 
  double dka, std::vector<double> kappaVals ){

  //clock_t t;
  //t = clock();
  int b_size = 3, a_size = 5;
  auto sab1 = ranges::view::generate_n([](){ return 0.0; },
                a_size*b_size) 
	    | ranges::view::chunk(b_size);


  //          | ranges::view::chunk(10)
  //          | ranges::view::transform([](auto range){ 
  //              return ranges::accumulate(range,0);
  //            });
  std::cout << sab1 << std::endl;


  //Eigen::MatrixXd matrix1 = Eigen::MatrixXd::Random(2,3);
  //Eigen::MatrixXd matrix1 = Eigen::MatrixXd(2,3);
  //Eigen::MatrixXd matrix1(2,3);
  //matrix1(0,0) = 100;
  //std::cout << matrix1 << std::endl;

  double bk = 8.617385e-5;
  double therm = 0.0253;

  // Loop over scatterers and temperatures
  int isecs = 0;
  double arat  = 1.0; // This is for scaling alpha and beta values later
  bool done = false;


  //std::vector<std::vector<std::vector<double>>> sym_sab( alpha.size(),
  //  std::vector<std::vector<double>> (beta.size(), 
  //  std::vector<double> ( ntempr, 0.0 ) ) );

  // This is only going to be used by coldh
  //std::vector<std::vector<std::vector<double>>> sym_sab_2( alpha.size(),
  //  std::vector<std::vector<double>> (beta.size(), 
  //  std::vector<double> ( ntempr, 0.0 ) ) );

  Eigen::Tensor<double,3> sym_sab_eigen(alpha.size(),beta.size(),ntempr);
  Eigen::Tensor<double,3> sym_sab_2_eigen(alpha.size(),beta.size(),ntempr);


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
       // std::cout << "we want to read in tempdependent parameters" << std::endl;
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


      //sym_sab_eigen.setZero();
      auto lambda_s_t_eff = contin( itemp, nphon, delta, tbeta, scaling, tev,
        sc, rho, alpha, beta, sym_sab_eigen);
      double lambda_s = std::get<0>(lambda_s_t_eff);


      //double t_contin = clock();

     // update the effective temperature list
      //t_eff_vec[itemp] = std::get<1>(lambda_s_t_eff) * temp;

 
      // Translational part of distribution, if any
      //std::cout << "\n-------- trans" << std::endl;
      if ( trans_weight > 0.0 ){
        trans( alpha, beta, trans_weight, delta, diffusion_const, sc, scaling,
          itemp, lambda_s, tbeta, t_eff_vec, temp_vec, sym_sab_eigen );
      }
      //double t_trans = clock();

      if ( oscEnergies.size() > 0 ){
      //std::cout << "\n-------- discre" << std::endl;
        std::vector<std::tuple<double,double>> oscEnergiesWeights(oscEnergies.size());
        for ( size_t i = 0; i < oscEnergies.size(); ++i ){
          oscEnergiesWeights[i] = std::make_tuple(oscEnergies[i],oscWeights[i]);
        }
        discre( itemp, sc, scaling, tev, lambda_s, trans_weight, tbeta, alpha,
          beta, temp_vec, oscEnergiesWeights, t_eff_vec, sym_sab_eigen );
      }
      //double t_discre = clock();

      std::cout << lambda_s << std::endl;

      /*
      std::cout << "contin: " << ((float)(t_contin-t))/CLOCKS_PER_SEC 
        << " seconds" << std::endl;
      std::cout << "trans:  " << ((float)(t_trans-t_contin))/CLOCKS_PER_SEC 
        << " seconds" << std::endl;
      std::cout << "discre: " << ((float)(t_discre-t_trans))/CLOCKS_PER_SEC 
        << " seconds" << std::endl;
      std::cout << "\n-------------------\nTotal: " << ((float)clock()-t)/CLOCKS_PER_SEC 
        << " seconds\n\n" << std::endl;
        */
      std::cout << "\n" << std::endl;
      //return;
      return sym_sab_eigen;

      if ( ncold != 0 ){
	bool free = false;
        coldh( itemp, temp, tev, ncold, trans_weight, tbeta, t_eff_vec, 
	  scaling, alpha, beta, dka, kappaVals, nbeta, lat, free, sym_sab_eigen, 
	  sym_sab_2_eigen );
      }


    }
    done = true;
  }

  //return;
    return sym_sab_eigen;
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

  return sym_sab_eigen;
}
