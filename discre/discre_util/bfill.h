#include <iostream>
#include <vector>


auto bfill( std::vector<double>& bex, std::vector<double>& rdbex, 
  std::vector<double>& betan ){

  int k = betan.size();
  for ( auto i = 0; i < betan.size(); ++i ){
    bex[i] = -betan[k-1];
    k = k - 1;
  }
  // If the first beta value is super small, set it to 1
  if ( betan[0] <= 1e-9 ){
    bex[betan.size()-1] = 0;
    k = betan.size() + 1;
  }else{
    // If it's not super small, copy this value into the next slot over
    k = betan.size() + 2;
    bex[betan.size()] = betan[0];
  }
  for( auto i = 1; i < betan.size(); ++i ){
    bex[k-1] = betan[i];
    k = k + 1;
  }
  
  for ( auto i = 0; i < k-2; ++i ){
    rdbex[i] = 1/(bex[i+1]-bex[i]);
  }
    
  return k - 1;
  // returning what leapr calls "nbx" 
}
