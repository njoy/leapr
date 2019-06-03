#include <iostream>
#include <vector>


auto bfill( std::vector<double>& bex, std::vector<double>& rdbex, 
  const std::vector<double>& betan ){
  /* bfill prepares the bex and rdbex vectors to be used by sint. 
   * Given some betan vector [a,b,c,d,e], It will turn it into either
   * [-e,-d,-c,-b,-a,a,b,c,d,e,0] or [-e,-d,-c,-b,0,b,c,d,e,0,0], depending on
   * whether a is greater than or less than 1.0e-9 (respectively).
   * rdbex  
   */
  int nbeta = betan.size();

  // Put the values of betan in the beginning of bex, in reverse order.
  // betan : [a,b,c,d,e] --> bex : [-e,-d,-c,-b,-a,0,0,0,0,0,0]
  int k = nbeta;
  for ( int i = 0; i < nbeta; ++i ){
    bex[i] = -betan[k-1];
    k = k - 1;
  }

  // If the first betan value is small, ignore it so we have a 
  // betan : [a,b,c,d,e] --> bex : [-e,-d,-c,-b,0,0,0,0,0,0,0]
  if ( betan[0] <= 1e-9 ){
    bex[nbeta-1] = 0;
    k = nbeta + 1;
  } 
  // If its not that small, we reflect keeping the first betan value, to be in
  // betan : [a,b,c,d,e] --> bex : [-e,-d,-c,-b,-a,a,0,0,0,0,0]
  else {
    k = nbeta + 2;
    bex[nbeta] = betan[0];
  }

  // Copy the positive side of the beta values in
  // betan : [a,b,c,d,e] --> bex : [-e,-d,-c,-b,0,b,c,d,e,0,0]
  // or
  // betan : [a,b,c,d,e] --> bex : [-e,-d,-c,-b,-a,a,b,c,d,e,0]
  for( int i = 1; i < nbeta; ++i ){
    bex[k-1] = betan[i];
    k = k + 1;
  }

  // Populate rdbex values
  for ( int i = 0; i < k-2; ++i ){
    rdbex[i] = 1/(bex[i+1]-bex[i]);
  }
    
  return k - 1;
  // returning what leapr calls "nbx" 
  // This is the number of entries in the bex and rdbex vectors before you get
  // to the 1 or 2 zeros at the end
}
