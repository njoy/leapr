#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "../../discre/discre_util/sint.h"
#include "betaLoop.h"

void equal( double a, double b ){
  //std::cout << a << "    "  << b << std::endl;
  if (b == 0.0){ 
    REQUIRE( b-a < 1e-5 );
    return;
  }
  REQUIRE ( std::abs( (a-b)/(b) ) < 1e-5 );
}

void equal_vec_mega_vec( std::vector<std::vector<std::vector<double>>> a, 
  std::vector<double> b ){
  REQUIRE( a.size()*a[0].size()*a[0][0].size() == b.size() );
  int i = 0;
  for ( auto a1 : a ){
    for ( auto a2 : a1 ){
      for ( auto a3 : a2 ){
        equal( a3, b[i] );
        i += 1;
      }
    }
  }
}

auto populateSymSab( int a_size, int b_size, bool is_normal ){
  std::vector<std::vector<std::vector<double>>> sym_sab(a_size,
    std::vector<std::vector<double>>(b_size,std::vector<double>(1,0.0)));
  int i = 1;
  for ( auto a = 0; a < a_size; ++a ){
    for ( auto b = 0; b < b_size; ++b ){
      if ( is_normal ){ sym_sab[a][b][0] = i; }
      else {sym_sab[a][b][0] = 0; }
      i += 1;
    }
  }
  return sym_sab;
}


TEST_CASE( "beta loop helper function" ){
  std::vector<double> 
    betan { 0.1, 0.15, 0.30, 0.60, 1.20 },
    rdbex { 1.6, 3.3, 6.6, 20, 5, 20, 6.6, 3.3, 1.6, 0, 0 },
    bex { -.02, -.006, -.003, -.0015, -.001, .001, .0015, .003, .006, .02, 0 },
    sex { 0.1, 0.2, 0.3, 0.5, 0.8, 1.3, 2.1, 3.4, 5.5, 8.9, 14.4 },
    correct_sym_sab (25, 0.0), correct_sym_sab_2 (25, 0.0);

  double alpha = 0.1, x = 3.2, y = 4.3, swe = 1.42, swo = 2.41, wt = 1.5, 
         tbart = 820;
  int itemp = 0, nbx = 10, a = 0, ncold = 1;

  auto sym_sab   = populateSymSab( 5, 5, true );
  auto sym_sab_2 = populateSymSab( 5, 5, false );


  GIVEN( "molecular translations are assumed to not be free" ){
    bool free = false;

    //------------------------------------------------------------------------
    correct_sym_sab = { 0.882913433, 0.649954077, 0.550346135, 0.394586983, 
      0.202842464, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 
      22, 23, 24, 25  };
    correct_sym_sab_2 = { 0.8829134, 0.9175721, 1.029927, 1.297595, 2.059706, 0, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
    THEN( "output scattering laws are correct" ){
      betaLoop( betan, rdbex, bex, sex, alpha, wt, tbart, x, y, swe, swo, 
        itemp, nbx, a, ncold, free, sym_sab, sym_sab_2 );
      equal_vec_mega_vec( sym_sab, correct_sym_sab );
      equal_vec_mega_vec( sym_sab_2, correct_sym_sab_2 );
    } // THEN
    //------------------------------------------------------------------------
    //------------------------------------------------------------------------
    
    bex = { -1.2, -0.6, -0.3, -0.15, -0.1, 0.1, 0.15, 0.3, 0.6, 1.2, 0.0 };

    correct_sym_sab = { 1.02763700, 0.395247529, 0.237150220, 0.158101736, 
      8.66802482E-2, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 
      21, 22, 23, 24, 25 }; 
    correct_sym_sab_2 = { 1.02763700, 1.66002658, 2.68765983, 4.34768245, 
      6.44628349, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
    THEN( "output scattering laws are correct" ){
      betaLoop( betan, rdbex, bex, sex, alpha, wt, tbart, x, y, swe, swo, 
        itemp, nbx, a, ncold, free, sym_sab, sym_sab_2 );
      equal_vec_mega_vec( sym_sab, correct_sym_sab );
      equal_vec_mega_vec( sym_sab_2, correct_sym_sab_2 );
    } // THEN
    //------------------------------------------------------------------------
    //------------------------------------------------------------------------

    sex = { 5, 4, 3, 2, 1, 1, 1, 2, 2, 1, 5 };
    x = 3.2; y = 4.3; swe = 1.42; swo = 2.41; wt = 1.5; tbart = 820;

    correct_sym_sab = { 0.790490930, 1.58097809, 2.37146522, 3.16195249, 
      3.70600979, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 
      22, 23, 24, 25 };
    correct_sym_sab_2 = { 0.790490930, 0.790490910, 1.58097789, 1.58097781, 
      0.790490664, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
 
    THEN( "output scattering laws are correct" ){
      betaLoop( betan, rdbex, bex, sex, alpha, wt, tbart, x, y, swe, swo, 
        itemp, nbx, a, ncold, free, sym_sab, sym_sab_2 );
      equal_vec_mega_vec( sym_sab, correct_sym_sab );
      equal_vec_mega_vec( sym_sab_2, correct_sym_sab_2 );
    } // THEN
    //------------------------------------------------------------------------
    //------------------------------------------------------------------------

    x = 0.85; y = 0.35; swe = 0.32; swo = 0.87; wt = 2.3; tbart = 950;

    correct_sym_sab = { 1.70316846, 3.37656089, 5.04886254, 6.69968895, 
      7.86359423, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 
      22, 23, 24, 25 };
    correct_sym_sab_2 = { 1.70316846, 1.70258735, 3.37159930, 3.36200407, 
      1.68666467, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
 
    THEN( "output scattering laws are correct" ){
      betaLoop( betan, rdbex, bex, sex, alpha, wt, tbart, x, y, swe, swo, 
        itemp, nbx, a, ncold, free, sym_sab, sym_sab_2 );
      equal_vec_mega_vec( sym_sab, correct_sym_sab );
      equal_vec_mega_vec( sym_sab_2, correct_sym_sab_2 );
    } // THEN
  } // GIVEN

  GIVEN( "molecular translations are assumed to be free" ){
    bool free = true;

    bex = { -1.2, -0.6, -0.3, -0.15, -0.1, 0.1, 0.15, 0.3, 0.6, 1.2, 0.0 };
    
    //------------------------------------------------------------------------
    correct_sym_sab =  { 0.518806235, 0.575763580, 0.554572267, 0.410837243, 
      9.16701822E-2, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 
      21, 22, 23, 24, 25 }; 
    correct_sym_sab_2 = { 0.5188062, 0.4955643, 0.4108372, 0.2254724, 
      2.762366E-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    THEN( "output scattering laws are correct" ){
      betaLoop( betan, rdbex, bex, sex, alpha, wt, tbart, x, y, swe, swo, 
        itemp, nbx, a, ncold, free, sym_sab, sym_sab_2 );
      equal_vec_mega_vec( sym_sab, correct_sym_sab );
      equal_vec_mega_vec( sym_sab_2, correct_sym_sab_2 );
    } // THEN
    //------------------------------------------------------------------------
    //------------------------------------------------------------------------

  } // GIVEN
} // TEST CASE
