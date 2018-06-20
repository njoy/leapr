#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "discre/discre_util/sint.h"
#include "coldh/coldh_util/betaLoop.h"
#include <unsupported/Eigen/CXX11/Tensor>

void checkSab( const Eigen::Tensor<double,3>& sab,
  const std::vector<double>& correctSab ){

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
}


auto populateSymSab( const std::vector<double>& alpha, const std::vector<double>& beta ){
  Eigen::Tensor<double,3> sab(alpha.size(),beta.size(),1);
  int k = 1;
  for ( int i = 0; i < sab.dimension(0); ++i ){
    for ( int j = 0; j < sab.dimension(1); ++j ){
      sab(i,j,0) = k;
      k += 1;
    }
  }
  return sab;
}



auto populateSymSab( int a_size, int b_size, bool is_normal ){
  Eigen::Tensor<double,3> sab(a_size,b_size,1);
  int k = 1;
  for ( int i = 0; i < sab.dimension(0); ++i ){
    for ( int j = 0; j < sab.dimension(1); ++j ){
      if ( is_normal ){ sab(i,j,0) = i; }
      else {sab(i,j,0) = 0; }
      k += 1;
    }
  }
  return sab;
}


TEST_CASE( "beta loop helper function" ){
  std::vector<double> betan, rdbex, bex, sex, goodSymSab1(25), goodSymSab2(25); 
  double alpha, x, y, swe, swo, wt, tbart;
  int itemp, nbx, a, ncold;
  bool free;

  betan = { 0.1, 0.15, 0.30, 0.60, 1.20 };
  rdbex = { 1.6, 3.3, 6.6, 20, 5, 20, 6.6, 3.3, 1.6, 0, 0 };
  bex = { -.02, -.006, -.003, -.0015, -.001, .001, .0015, .003, .006, .02, 0 };
  sex = { 0.1, 0.2, 0.3, 0.5, 0.8, 1.3, 2.1, 3.4, 5.5, 8.9, 14.4 };

  alpha = 0.1, x = 3.2, y = 4.3, swe = 1.42, swo = 2.41, wt = 1.5, tbart = 820;
  itemp = 0, nbx = 10, a = 0, ncold = 1;

  auto sym_sab   = populateSymSab( 5, 5, true );
  auto sym_sab_2 = populateSymSab( 5, 5, false );


  GIVEN( "molecular translations are assumed to not be free" ){
    free = false;

    //------------------------------------------------------------------------
    goodSymSab1 = { 0.882913433, 0.649954077, 0.550346135, 0.394586983, 
      0.202842464, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 
      22, 23, 24, 25 };
    goodSymSab2 = { 0.8829134, 0.9175721, 1.029927, 1.297595, 2.059706, 0, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
    THEN( "output scattering laws are correct" ){
      betaLoop( betan, rdbex, bex, sex, alpha, wt, tbart, x, y, swe, swo, 
        itemp, nbx, a, ncold, free, sym_sab, sym_sab_2 );
      checkSab( sym_sab, goodSymSab1 );
      checkSab( sym_sab_2, goodSymSab2 );
    } // THEN
    //------------------------------------------------------------------------
    //------------------------------------------------------------------------
    
    bex = { -1.2, -0.6, -0.3, -0.15, -0.1, 0.1, 0.15, 0.3, 0.6, 1.2, 0.0 };

    goodSymSab1 = { 1.02763700, 0.395247529, 0.237150220, 0.158101736, 
      8.66802482E-2, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 
      21, 22, 23, 24, 25 }; 
    goodSymSab2 = { 1.02763700, 1.66002658, 2.68765983, 4.34768245, 
      6.44628349, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
    THEN( "output scattering laws are correct" ){
      betaLoop( betan, rdbex, bex, sex, alpha, wt, tbart, x, y, swe, swo, 
        itemp, nbx, a, ncold, free, sym_sab, sym_sab_2 );
      checkSab( sym_sab, goodSymSab1 );
      checkSab( sym_sab_2, goodSymSab2 );
    } // THEN
    //------------------------------------------------------------------------
    //------------------------------------------------------------------------

    sex = { 5, 4, 3, 2, 1, 1, 1, 2, 2, 1, 5 };
    x = 3.2; y = 4.3; swe = 1.42; swo = 2.41; wt = 1.5; tbart = 820;

    goodSymSab1 = { 0.790490930, 1.58097809, 2.37146522, 3.16195249, 
      3.70600979, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 
      22, 23, 24, 25 };
    goodSymSab2 = { 0.790490930, 0.790490910, 1.58097789, 1.58097781, 
      0.790490664, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
 
    THEN( "output scattering laws are correct" ){
      betaLoop( betan, rdbex, bex, sex, alpha, wt, tbart, x, y, swe, swo, 
        itemp, nbx, a, ncold, free, sym_sab, sym_sab_2 );
      checkSab( sym_sab, goodSymSab1 );
      checkSab( sym_sab_2, goodSymSab2 );
    } // THEN
    //------------------------------------------------------------------------
    //------------------------------------------------------------------------

    x = 0.85; y = 0.35; swe = 0.32; swo = 0.87; wt = 2.3; tbart = 950;

    goodSymSab1 = { 1.70316846, 3.37656089, 5.04886254, 6.69968895, 
      7.86359423, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 
      22, 23, 24, 25 };
    goodSymSab2 = { 1.70316846, 1.70258735, 3.37159930, 3.36200407, 
      1.68666467, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
 
    THEN( "output scattering laws are correct" ){
      betaLoop( betan, rdbex, bex, sex, alpha, wt, tbart, x, y, swe, swo, 
        itemp, nbx, a, ncold, free, sym_sab, sym_sab_2 );
      checkSab( sym_sab, goodSymSab1 );
      checkSab( sym_sab_2, goodSymSab2 );
    } // THEN
  } // GIVEN

  GIVEN( "molecular translations are assumed to be free" ){
    bool free = true;

    bex = { -1.2, -0.6, -0.3, -0.15, -0.1, 0.1, 0.15, 0.3, 0.6, 1.2, 0.0 };
    
    //------------------------------------------------------------------------
    goodSymSab1 =  { 0.518806235, 0.575763580, 0.554572267, 0.410837243, 
      9.16701822E-2, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 
      21, 22, 23, 24, 25 }; 
    goodSymSab2 = { 0.5188062, 0.4955643, 0.4108372, 0.2254724, 
      2.762366E-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    THEN( "output scattering laws are correct" ){
      betaLoop( betan, rdbex, bex, sex, alpha, wt, tbart, x, y, swe, swo, 
        itemp, nbx, a, ncold, free, sym_sab, sym_sab_2 );
      checkSab( sym_sab, goodSymSab1 );
      checkSab( sym_sab_2, goodSymSab2 );
    } // THEN
    //------------------------------------------------------------------------
    //------------------------------------------------------------------------

  } // GIVEN
} // TEST CASE
