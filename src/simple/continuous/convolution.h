#include "simple/continuous/interpolate.h"
#include <iostream>

void test(std::vector<double> z){
  std::cout << std::endl;
  for (auto& y : z){ std::cout << y << " "; }
  std::cout << std::endl;
}





template <typename V>
V reflect(V v){
  V fullVec( 2*v.size()-1 );
  std::reverse_copy (v.begin(), v.begin()+v.size(), fullVec.begin());
  std::copy (v.begin(), v.begin()+v.size(), fullVec.begin()+int(fullVec.size()/2));

  return fullVec;

}


template <typename V>
V convolve(V x_0, V y1_0, V y2_0){ 
  V x  = reflect(x_0);
  V y1 = reflect(y1_0);
  V y2 = reflect(y2_0);
  for ( int i = 0; i < int(x.size()/2); ++i ){ x[i] *= -1; }
  for ( int i = 0; i < int(y1.size()/2); ++i ){ y1[i] *= exp(-x[i]); }
  for ( int i = 0; i < int(y2.size()/2); ++i ){ y2[i] *= exp(-x[i]); }

  //test(x);
  //test(y1);
  //test(y2);
  V y3 (x_0.size(),0.0);

  for (size_t i = 0; i < x_0.size(); ++i ){
    for (size_t j = 0; j < x.size()-1; ++j ){
      auto delta = x[j+1]-x[j];
      //std::cout << j << "  beta  " << x[j] << "     delta  " << delta << std::endl;
      //
      //std::cout << interpolate(x, y1, x[j]) << "     " << y1[j] << std::endl;
      //std::cout << interpolate(x, y2, x[y2.size()-j-1]) << "     " << y2[y2.size()-j-1] << std::endl;
      //std::cout << std::endl;
      //y3[i] += delta*0.5*(y1[j]*y2[y2.size()-j-1]) + 
      //         delta*0.5*(y1[j+1]*y2[y2.size()-j-2]);
      y3[i] += delta*0.5*(interpolate(x,y1,x[j])*interpolate(x,y2,x[y2.size()-j-1])) + 
               delta*0.5*(interpolate(x,y1,x[j+1])*interpolate(x,y2,x[y2.size()-j-2]));

      std::cout << x[j] << "     " << x[y2.size()-j-1] << "      " << x[i] << "     " << x_0[i]-x[j] << std::endl;
      //std::cout << delta*0.5*(y1[j]*y2[y2.size()-j-1]) <<  "       " << \
      //  delta*0.5*(y1[j+1]*y2[y2.size()-j-2]) << std::endl;
      //test(y3);
    




    }
    break;
  }
  return y3;
  std::cout << y1.size()+y2.size()+y3.size() << std::endl;
}

