#include "catch.hpp"
#include "discre/discre_util/exts.h"
#include "testFuncs.h"



TEST_CASE( "exts" ){
  std::vector<double> beta { 0, 0.15, 0.3, 0.6, 1.2 },
    sexpb { 53.06972, 9.052783E-2, 2.393197E-2, 6.549547E-3, 1.876351E-3 },
    exb { 1, 0.9277434, 0.8607079, 0.7408182, 0.5488116, 0, 0, 0, 0, 0, 0 },
    sex;
  

  GIVEN( "First beta value is really small (<=1e-9)" ){
    beta = { 9e-10, 2.15, 3.30, 4.60, 5.20 };

    WHEN( "exb values are medium size" ){
      sexpb = { 11.9, 11.7, 11.1, 9, 6.3 };
      exb = { 1, 0.9936967987, 0.987433327, 0.975024577, 0.962771762, 0,
        0, 0, 0, 0, 0 };
      double correctSex[11] = { 6.3, 9, 11.1, 11.7, 11.9, 11.6262523, 10.9605103,
        8.77522119, 6.06546228, 0, 0 };

      THEN( "output sex vector is correct" ){
        auto sex = exts( sexpb, exb, beta );
        checkVectorAgainstRange(sex,correctSex);
      } // THEN
    } // WHEN

    WHEN( "exb values are smaller" ){
      sexpb = { 11.9, 11.7, 11.1, 9, 6.3 };
      exb = { 1.21, 0.25, 4.0e-2, 2.56e-2, 4.0e-4, 0, 0, 0, 0, 0, 0 };
      double correctSex[11] = { 6.3, 9, 11.1, 11.7, 11.9, 2.925, 0.444, 0.2303999, 
        2.52E-3, 0, 0 };

      THEN( "output sex vector is correct" ){
        auto sex = exts( sexpb, exb, beta );
        checkVectorAgainstRange(sex,correctSex);
      } // THEN
    } // WHEN
  } // GIVEN

  GIVEN( "First beta value is not that small (>1e-9)" ){
    WHEN( "exb values are medium size" ){
      beta  = { 0.1, 0.15, 0.3, 0.6, 1.2 };
      sexpb = { 11.9, 11.7, 11.1, 9, 6.3 };
      exb = { 1, 0.9936967987, 0.987433327, 0.975024577, 0.962771762, 0,
        0, 0, 0, 0, 0 };
      double correctSex[11] = { 6.3, 9, 11.1, 11.7, 11.9, 11.9, 11.6262523, 
        10.960510, 8.7752211, 6.06546228, 0 };

      THEN( "output sex vector is correct" ){
        auto sex = exts( sexpb, exb, beta );
        checkVectorAgainstRange(sex,correctSex);
      } // THEN

      beta = { 1.1, 2.15, 3.30, 4.60, 5.20 };

      THEN( "output sex vector is correct" ){
        auto sex = exts( sexpb, exb, beta );
        checkVectorAgainstRange(sex,correctSex);
      } // THEN
    } // WHEN

    WHEN( "exb values are smaller" ){
      sexpb = { 11.9, 11.7, 11.1, 9, 6.3 };
      beta = { 1.1, 2.15, 3.30, 4.60, 5.20 };
      exb = { 1.21, 0.25, 4.0E-2, 2.56E-2, 4.0E-4, 0, 0, 0, 0, 0, 0 };
      double correctSex[11] = { 6.3, 9, 11.1, 11.7, 11.9, 11.9, 2.92499995, 0.444, 
        0.2303999, 2.52E-3, 0 };

      THEN( "output sex vector is correct" ){
        auto sex = exts( sexpb, exb, beta );
        checkVectorAgainstRange(sex,correctSex);
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
