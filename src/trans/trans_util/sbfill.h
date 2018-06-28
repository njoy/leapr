#include <iostream>
#include <vector>
#include <cmath>
#include <range/v3/all.hpp>

void sbfill(std::vector<double>& sb, int nbt, double delta, double be, 
  std::vector<double>& s,std::vector<double>& betan, int ndmax){

  double bmin = -be - (nbt-1) * delta;
  double bmax = -be + (nbt-1) * delta + delta * 0.01;

  if ( 1 + (bmax-bmin) / delta > ndmax){ 
    throw std::exception();
  }
  
  int i = 0;
  unsigned int j = betan.size()-1;
  double current, toLeft, bet = bmin, slim = -225e0;
  //bool indexInRange = false; 
  
  while (bet < bmax){
    std::cout << bet << std::endl;
    double b = std::abs(bet);
    do {
      if ( b >  betan[j] and (j+1) != betan.size() )    { ++j; }
      if ( b <= betan[j] and  b < betan[j-1] and j > 1 ){ --j; }
    }  
    while (
         (b >  betan[j] or ( b < betan[j-1] and j > 1 ) ) and 
         (b <= betan[j] or ( j+1 != betan.size() ) ) );

    /*
    if ( b > betan[j] and  j + 1 == betan.size() ){ 
      indexInRange = b < 1.00001 * betan[j] ? true : false;
    }
    if ( b <= betan[j] ){ indexInRange = true; }
    */

    
    if ( (b > betan[j] and ( j+1!=betan.size() or b >= 1.00001*betan[j]) ) ){
      sb[i] = 0; 
    }
    else { 
      current = s[j]   < 0 ? slim : log( s[j]   );
      toLeft  = s[j-1] < 0 ? slim : log( s[j-1] );
      sb[i] = current + (b-betan[j])*(toLeft-current)/(betan[j-1]-betan[j]);
      if (bet > 0) { sb[i] = sb[i] - bet; }
      if ( sb[i] > slim ){ sb[i] = exp(sb[i]); }
    }

    i = i + 1;
    bet = bet + delta;
  }
  return;
  auto betRange = ranges::view::iota(0,int((bmax - bmin)/delta)+1) 
                | ranges::view::transform([bmin,delta](int i){
                    return bmin + (double(i)*delta); } );
  std::cout << betRange << std::endl;
}











