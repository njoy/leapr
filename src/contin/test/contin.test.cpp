#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "contin/contin.h"


template <typename A, typename B, typename C, typename D>
void checkSabLambdaTeff( const A& correctSab, 
    const B& output, 
    const C& sab,  
    const D& lambda, const D& teff ){

  REQUIRE( sab.dimension(0)*sab.dimension(1)*sab.dimension(2) == correctSab.size() );

  int l = 0;
  for ( int i = 0; i < sab.dimension(0); ++i ){
    for ( int j = 0; j < sab.dimension(1); ++j ){
      for ( int k = 0; k < sab.dimension(2); ++k ){
        REQUIRE( sab(i,j,k) == Approx(correctSab[l]).epsilon(1e-6) );
	l += 1;
      }
    }
  }

  REQUIRE( std::get<0>(output) == Approx(lambda).epsilon(1e-6) );
  REQUIRE( std::get<1>(output) == Approx(teff).epsilon(1e-6) );
}


TEST_CASE( "contin eigen" ){

  REQUIRE(true);
  int ntempr, nphon, itemp;
  double delta, tbeta, tev, sc, scaling, lambda_s, t_eff;
  std::vector<double> alpha, beta, rho, expected;

  GIVEN( "input values from input card and leapr subroutine" ){


    ntempr = 1; nphon = 3; itemp = 0;
    delta = 0.1; tbeta = 1.0; tev = 0.01723477; sc = 0.0253/tev; scaling = sc;
    alpha = { 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28};
    beta  = { 0.00, 0.15, 0.30, 0.60, 1.20 };
    rho   = { 0.002, 0.004, 0.02, 0.04, 0.2, 0.4 };


    WHEN( "3rd order expansion, with alpha & beta vals scaled by 0.0253/tev" ){

      Eigen::Tensor<double,3> symSab( alpha.size(), beta.size(), ntempr );

      auto output = contin( itemp, nphon, delta, tbeta, scaling, tev, sc, rho, 
          alpha, beta, symSab );

      THEN( "contin output matches expected value" ){
        lambda_s = 4.38473153E-2, t_eff = 12.279863466;
        expected = {6.474331963E-7, 7.658564709E-7, 8.842797456E-7, 
          1.121126294E-6, 1.594819393E-6, 1.294036143E-6, 1.530731302E-6, 
          1.767426462E-6, 2.240816781E-6, 3.187597419E-6, 2.584754483E-6, 
          3.057539823E-6, 3.530325163E-6, 4.475895844E-6, 6.367037205E-6, 
          5.156263261E-6, 6.099418664E-6, 7.042574067E-6, 8.928884872E-6, 
          1.270150648E-5, 1.025974716E-5, 1.213643384E-5, 1.401312052E-5, 
          1.776649387E-5, 2.527324059E-5, 2.030999488E-5, 2.402516585E-5, 
          2.774033682E-5, 3.517067875E-5, 5.003136263E-5, 3.979478204E-5, 
          4.707463949E-5, 5.435449694E-5, 6.891421183E-5, 9.803364162E-5, 
          7.638863877E-5, 9.036457859E-5, 1.043405184E-4, 1.322923980E-4, 
          1.881961573E-4};
        checkSabLambdaTeff( expected, output, symSab, lambda_s, t_eff );
      } // THEN
    } // WHEN

    WHEN( "6th order exp, alpha & beta vals scaled, and small grid space" ){
      nphon = 6; delta = 0.04; sc = 1.0; scaling = 1.0;
      alpha =  { 0.1, 0.2, 0.4, 0.8, 1.6 };

      Eigen::Tensor<double,3> symSab( alpha.size(), beta.size(), ntempr );

      auto output = contin( itemp, nphon, delta, tbeta, scaling, tev, sc, rho, 
          alpha, beta, symSab );

      THEN( "contin output matches expected value" ){
        lambda_s = 0.11157823, t_eff = 4.91699518;
        expected = {6.82096404E-5, 7.51470660E-5, 8.20844916E-5, 9.59593429E-5,
          1.23709045E-4, 1.34940719E-4, 1.48666045E-4, 1.62391372E-4, 
          1.89842025E-4, 2.44743331E-4, 2.64063030E-4, 2.90925221E-4, 
          3.17787412E-4, 3.71511794E-4, 4.78960557E-4, 5.05599918E-4, 
          5.57045626E-4, 6.08491333E-4, 7.11382748E-4, 9.17165579E-4, 
          9.26780041E-4, 1.02112863E-3, 1.11547723E-3, 1.30417442E-3, 
          1.68156881E-3};
	checkSabLambdaTeff( expected, output, symSab, lambda_s, t_eff );
      } // THEN
    } // WHEN

    WHEN( "6th order exp, user-defined normalizationand large grid space" ){
      delta = 4.; tbeta = 2.0; sc = 1.0; scaling = 1.0;
      alpha = { 0.1, 0.2, 0.4, 0.8, 1.6 };

      Eigen::Tensor<double,3> symSab( alpha.size(), beta.size(), ntempr );
      symSab.setZero();

      auto output = contin( itemp, nphon, delta, tbeta, scaling, tev, sc, rho, 
        alpha, beta, symSab );

      THEN( "contin output matches expected value" ){

        lambda_s = 2.179428E-3; t_eff = 491.1882921;
        expected = {1.37883996E-10, 1.58477481E-10, 1.79070966E-10, 
          2.20257936E-10, 3.02631876E-10, 2.75707903E-10, 3.16885898E-10, 
          3.58063894E-10, 4.40419886E-10, 6.05131869E-10, 5.51175522E-10, 
          6.33495628E-10, 7.15815735E-10, 8.80455949E-10, 1.20973637E-09, 
          1.10139053E-09, 1.26588730E-09, 1.43038407E-09, 1.75937760E-09, 
          2.41736468E-09, 2.19894405E-09, 2.52736456E-09, 2.85578506E-09, 
          3.51262608E-09, 4.82630810E-09};
	checkSabLambdaTeff( expected, output, symSab, lambda_s, t_eff );
      } // THEN
    } // WHEN
  } // GIVEN 
} // TEST CASE
