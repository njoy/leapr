#define CATCH_CONFIG_MAIN
#include "coh/coh_util/readWrite_util/loada.h"


TEST_CASE( "loada" ){
  GIVEN( "" ){

    WHEN( "na = 2" ){
      int i = 0, na = 2, nbuf = 10;
      std::fstream ntape("test_loada");
      std::vector<double> a(20,0.0), buf { 2, 4, 6, 8, 10, 12, 14, 16, 18, 20 };
      a[0] = 1e-5; a[1] = 2e-5; a[2] = 3e-5; a[3] = 4e-5;
      loada( i, na, ntape, nbuf, a, buf );

      std::vector<double> bufVals { 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 1e-5, 2e-5 };
      for ( size_t i = 0; i < bufVals.size(); ++i ){
        REQUIRE( bufVals[i] == Approx(buf[i]).epsilon(1e-6) );
      }
      ntape.close();
    } // WHEN
    WHEN( "na = 3" ){
      int i = 0, na = 3, nbuf = 10;
      std::fstream ntape("test_loada");
      std::vector<double> a(20,0.0), buf { 2, 4, 6, 8, 10, 12, 14, 16, 18, 20 };
      a[0] = 1e-5; a[1] = 2e-5; a[2] = 3e-5; a[3] = 4e-5;
      loada( i, na, ntape, nbuf, a, buf );

      std::vector<double> bufVals { 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 1e-5, 2e-5, 3e-5, 20.0 };
      for ( size_t i = 0; i < bufVals.size(); ++i ){
        REQUIRE( bufVals[i] == Approx(buf[i]).epsilon(1e-6) );
      }
      ntape.close();
    } // WHEN
    WHEN( "na = 3" ){
      int i = 11, na = 2, nbuf = 10;
      std::fstream ntape("test_loada");
      std::vector<double> a(20,0.0), buf { 2, 4, 6, 8, 10, 12, 14, 16, 18, 20 };
      a[0] = 1e-5; a[1] = 2e-5; a[2] = 3e-5; a[3] = 4e-5;
      loada( i, na, ntape, nbuf, a, buf );

      std::vector<double> bufVals { 1e-5, 2e-5, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0 };
      for ( size_t i = 0; i < bufVals.size(); ++i ){
        REQUIRE( bufVals[i] == Approx(buf[i]).epsilon(1e-6) );
      }
      ntape.close();
    } // WHEN
  } // GIVEN
  GIVEN( "" ){
    int i = 1, na = 3, nbuf = 10;
    std::fstream ntape("test_loada");
    std::vector<double> a(20,0.0), buf { 2, 4, 6, 8, 10, 12, 14, 16, 18, 20 };
    a[0] = 1e-5; a[1] = 2e-5; a[2] = 3e-5; a[3] = 4e-5;
    loada( i, na, ntape, nbuf, a, buf );

    std::vector<double> bufVals { 1e-5, 2e-5, 3e-5, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0 };
    for ( size_t i = 0; i < bufVals.size(); ++i ){
      REQUIRE( bufVals[i] == Approx(buf[i]).epsilon(1e-6) );
    }
    ntape.close();
  } // GIVEN
} // TEST CASE
