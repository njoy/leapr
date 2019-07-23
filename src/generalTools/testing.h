
template <typename V>
void checkVec( V y1, V y2, float tol=1e-6 ){
  REQUIRE( y1.size() == y2.size() );
  for (size_t i = 0; i < y1.size(); ++i){
    REQUIRE( y2[i] == Approx(y1[i]).epsilon(tol) );
  }
}

template <typename V, typename Int>
void checkVec( V y1, V y2, Int end, float tol=1e-6 ){
  for (size_t i = 0; i < end; ++i){
    REQUIRE( y2[i] == Approx(y1[i]).epsilon(tol) );
  }
}



template <typename V>
void restAreZero(int n, const V& vec){
  for (size_t i = n; i < vec.size(); ++i){
    REQUIRE( 0.0 == Approx(vec[i]).epsilon(1e-6) );
  }
}


