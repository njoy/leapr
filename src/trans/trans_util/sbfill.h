#include <iostream>
#include <vector>
#include <cmath>
#include <range/v3/all.hpp>

template <typename floatT, typename arrayT>
void sbfill(arrayT& sb, int nbt, floatT delta, floatT be, arrayT& s, 
  arrayT& betan, int ndmax){

  floatT bmin = -be - (nbt-1) * delta,
         bmax = -be + (nbt-1) * delta + delta * 0.01;

  if ( 1 + (bmax-bmin) / delta > ndmax){ 
    throw std::exception();
  }
  
  //int i = 0;
  unsigned int j = betan.size()-1;
  floatT current, toLeft, /*bet = bmin,*/ slim = -225e0;
  

  // We are going to have round up(bmax-bmin/delta) many iterations. 
  // So track the iteration number. And zip this with a beta range,
  // which goes from bmin --> bmax in steps of size delta
  auto iRange   = ranges::view::iota(0,int((bmax - bmin)/delta)+1);
  auto iBetaRange = ranges::view::zip( 
                      iRange, 
                      iRange | ranges::view::transform( [bmin,delta](int i){
                                 return bmin + (1.0*i*delta); } ) );


  for ( const auto& tupleEntry : iBetaRange ){
    int i = std::get<0>(tupleEntry);
    floatT bet = std::get<1>(tupleEntry);
    floatT b = std::abs(bet);

    do {
      if ( b >  betan[j] and j+1 != betan.size() )    { ++j; }
      if ( b <= betan[j] and  b < betan[j-1] and j > 1 ){ --j; }
    }  
    while (
         (b >  betan[j] or ( b < betan[j-1] and j > 1 ) ) and 
         (b <= betan[j] or ( j+1 != betan.size() ) ) 
    );

    
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
  }

}











