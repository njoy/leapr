

template <typename Vec, typename Array>
void checkAgainstArray( const Array& correct, const Vec& output ){
  // correct : C++ array (eg. double x[3]) that has the correct values
  // output  : vector (eg. std::vector<double>) that our function spit out.
  // If correct < output, we assume that all other values in output should be
  // equal to zero
  for ( size_t i = 0; i < output.size(); ++i ){
    if ( i < (sizeof(correct)/sizeof(*correct)) ) {
      REQUIRE( correct[i] == Approx(output[i]).epsilon(1e-6) );
    }
    else { REQUIRE( 0.0 == Approx(output[i]).epsilon(1e-6) ); }
  }
}


