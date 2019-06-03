
template <typename V>
void checkVec( V y1, V y2, float tol=1e-6 ){
  for (size_t i = 0; i < y1.size(); ++i){
    REQUIRE( y2[i] == Approx(y1[i]).epsilon(tol) );
  }
}


