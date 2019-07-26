#include <iostream>
#include <vector>
#include <range/v3/all.hpp>

/* bfill prepares the bex and rdbex vectors to be used by sint. 
 * Given some betan vector [a,b,c,d,e], It will turn it into either
 * [-e,-d,-c,-b,-a,a,b,c,d,e,0] or [-e,-d,-c,-b,0,b,c,d,e,0,0], depending on
 * whether a is greater than or less than 1.0e-9 (respectively).
 */
template <typename Range>
auto bfill_new( int& bexSize, Range& rdbex, const Range& betan ){
  using std::begin; using std::end;
  // betan : [a,b,c,d,e] --> bex : [-e,-d,-c,-b,-a,0,0,0,0,0,0]
  Range bex(bexSize,0.0);
  std::reverse_copy(begin(betan),end(betan),begin(bex));
  bex = bex | ranges::view::transform([](auto x){return -x;});

  // If the first betan value is small, ignore it so we have a 
  // betan : [a,b,c,d,e] --> bex : [-e,-d,-c,-b,0,0,0,0,0,0,0]
  // If its not that small, we reflect keeping the first betan value, to be in
  // betan : [a,b,c,d,e] --> bex : [-e,-d,-c,-b,-a,a,0,0,0,0,0]
  if ( betan[0] <= 1e-9 ){ bex[betan.size()-1] = 0;        } 
  else {                   bex[betan.size()  ] = betan[0]; }

  // Copy the positive side of the beta values in
  // betan : [a,b,c,d,e] --> bex : [-e,-d,-c,-b,0,b,c,d,e,0,0]   or
  // betan : [a,b,c,d,e] --> bex : [-e,-d,-c,-b,-a,a,b,c,d,e,0]
  if (betan[0] <= 1e-9){
    std::copy(begin(betan)+1,end(betan),begin(bex)+betan.size());
  } else{
    std::copy(begin(betan),end(betan),begin(bex)+betan.size());
  }

  size_t k1 = (betan[0] <= 1e-9) ? 2*betan.size()-1 : 2*betan.size() ;
  for ( size_t i = 0; i < k1-1; ++i ){ rdbex[i] = 1.0 / (bex[i+1]-bex[i]); }
 
  return std::make_tuple(k1,bex);
  // returning what leapr calls "nbx" 
  // This is the number of entries in the bex and rdbex vectors before you get
  // to the 1 or 2 zeros at the end
}













template <typename Range>
auto bfill( Range& bex, Range& rdbex, const Range& betan ){
  using std::begin; using std::end;
  /* bfill prepares the bex and rdbex vectors to be used by sint. 
   * Given some betan vector [a,b,c,d,e], It will turn it into either
   * [-e,-d,-c,-b,-a,a,b,c,d,e,0] or [-e,-d,-c,-b,0,b,c,d,e,0,0], depending on
   * whether a is greater than or less than 1.0e-9 (respectively).
   */

  // Put the values of betan in the beginning of bex, in reverse order.
  // betan : [a,b,c,d,e] --> bex : [-e,-d,-c,-b,-a,0,0,0,0,0,0]
  Range bex2(bex.size(),0.0);
  std::reverse_copy(begin(betan),end(betan),begin(bex2));
  Range bex3 = bex2 | ranges::view::transform([](auto x){return -x;});

  // If the first betan value is small, ignore it so we have a 
  // betan : [a,b,c,d,e] --> bex : [-e,-d,-c,-b,0,0,0,0,0,0,0]
  // If its not that small, we reflect keeping the first betan value, to be in
  // betan : [a,b,c,d,e] --> bex : [-e,-d,-c,-b,-a,a,0,0,0,0,0]
  if ( betan[0] <= 1e-9 ){ bex3[betan.size()-1] = 0;        } 
  else {                   bex3[betan.size()  ] = betan[0]; }

  // Copy the positive side of the beta values in
  // betan : [a,b,c,d,e] --> bex : [-e,-d,-c,-b,0,b,c,d,e,0,0]   or
  // betan : [a,b,c,d,e] --> bex : [-e,-d,-c,-b,-a,a,b,c,d,e,0]
  if (betan[0] <= 1e-9){
    std::copy(begin(betan)+1,end(betan),begin(bex3)+betan.size());
  } else{
    std::copy(begin(betan),end(betan),begin(bex3)+betan.size());
  }

  size_t k1 = (betan[0] <= 1e-9) ? 2*betan.size()-1 : 2*betan.size() ;
  for ( size_t i = 0; i < k1-1; ++i ){
    rdbex[i] = 1.0 / (bex3[i+1]-bex3[i]);
  }
 
  for ( size_t i = 0; i < bex.size(); ++i ){ bex[i] = bex3[i]; }

  return std::make_tuple(k1,bex3);
  // returning what leapr calls "nbx" 
  // This is the number of entries in the bex and rdbex vectors before you get
  // to the 1 or 2 zeros at the end
}
