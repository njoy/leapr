#include "catch.hpp"
#include "endout/endout_util/tpidio.h"


TEST_CASE( "tpidio" ){
  int nin, nout, nscr, nb, nw, nsc, nsh, mth, mfh, math;
  nin = 0;
  nout = 24;
  nscr = 0;
  nb = 537;
  nw = 0;
  nsc = 0;
  nsh = 0;
  mth = 0;
  mfh = 0;
  math = 1;
  tpidio( nin, nout, nscr, nb, nw, nsc, nsh, mth, mfh, math );
  REQUIRE( true );
} // TEST CASE
