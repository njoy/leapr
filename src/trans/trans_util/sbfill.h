#include <iostream>
#include <vector>
#include <cmath>
#include <range/v3/all.hpp>

template <typename floatT, typename arrayT>
void sbfill(arrayT& sb, int nbt, floatT delta, floatT be, arrayT& s, 
  arrayT& betan, int ndmax){
  /* Overview
   * ------------------------------------------------------------------------
   * This is used, according to the manual pg. 711, to remap the current 
   * scattering law onto the same beta grid that the trans law is prepared on.
   * This is so that when trans later does convolution on the existant and 
   * translational scattering laws, it's easier. 
   *
   * Inputs
   * ------------------------------------------------------------------------
   * sb      : empty vector with ndmax ( max[ beta size, 1e6 ] )
   * nbt     : number of iterations that s_table_generation took to converge
   *           so this is basically the number of nonzero values in sd vector
   * delta   : spacing used while calculating sd vector (which in the manual
   *           is S_t(a,b)). This was calculated early on in trans.h
   * be      : beta value
   * s       : vector of existing S(a,b) values --> [S(a,b0), S(a,b1), ...]
   * betan   : beta vector
   * ndmax   : max[ beta size, 1e6 ]
   * 
   * Operations
   * ------------------------------------------------------------------------
   * * Establish the range bmin --> bmax across which we want to create a 
   *   unionized grid. We take -beta and go nbt*delta to the left and 
   *   nbt*delta to the right, since nbt is the number of S_t(a,b) values 
   *   that we computed in s_table_generation.
   * * We're going to march across the bmin, bmax range. bet is the value
   *   between bmin and bmax that we're currently at, and b is it's abs. val.
   * * Given a normalized beta grid betan [ beta0, beta1, betan2, ... ]
   *   find the pair ( beta_i, beta_(i+1) ) such that b lies between them. 
   *   This is done in the do while loop.
   * * Once we know where this s_table_generation's beta value roughly lies, 
   *   we log interpolate to find its corresponding value, and then
   *   have the S(a,b) value that corresponds to it. 
   *   So if I find that this value of b is in between beta3 and beta4, then
   *   I record into sb[i] 
   *       exp(  linearly interpolate log( S(a,b3) ) and log( S(a,b4) ) )
   * * March forward so that you do this for every bet/b value between bmin
   *   and bmax. If one of these bet/b values is outside betan range, set the
   *   sb[i] value equal to zero.
   * 
   * Outputs
   * ------------------------------------------------------------------------
   * * I have a set of beta values that interest me, that range from bmin to
   *   bmax. This range is centered about a specific beta value, and its size is
   *   determined by nbt (the number of interesting S_t(a,b) values from 
   *   s_table_generation). 
   * * For each beta value in this specific grid, I interpolate to find its
   *   approximate value in my existing S(a,b) table that was generated from
   *   contin. 
   * * So sb is a vector that through sbfill is populated with S(a,bi) where 
   *   bi is the ith beta value across this specific grid
   */

  // nbt is the number of S_t(a,b) values that were computed in 
  // s_table_generation. So we want to take -beta and look at the values that
  // lie in the +- nbt interval. Since s_table_generation used a spacing of 
  // delta (this is the delta that was computed in trans.h), we need to scale
  // our steps by that.
  floatT bmin = -be - (nbt-1) * delta,
         bmax = -be + (nbt-1) * delta + delta * 0.01;

  int numIters = (bmax-bmin)/delta + 1;
  // This is the same as 
  // if ( 2 * nbt - 0.99 > ndmax ){
  // which is effectively the same as 2*nbt > ndmax, since ndmax is an int
  if (numIters > ndmax){ 
    throw std::exception();
}
  
  unsigned int j = betan.size()-1;
  floatT slim = -225e0;
  

  // We are going to have round up(bmax-bmin/delta) many iterations. 
  // So track the iteration number. And zip this with a beta range,
  // which goes from bmin --> bmax in steps of size delta
  auto iBetaRange = ranges::view::zip( 
                      ranges::view::iota(0,numIters), 
                      ranges::view::iota(0,numIters) 
		    | ranges::view::transform( [bmin,delta](int i){
                        return bmin + (1.0*i*delta); } ) );


  for ( const auto& tupleEntry : iBetaRange ){
    int i = std::get<0>(tupleEntry);
    floatT bet = std::get<1>(tupleEntry);
    floatT b = std::abs(bet);

    do {
      // If desired point is to the right of current point, and I still have 
      // plenty of beta values to the right that I can explore, I'll just 
      // increase my index and keep searching. 
      // If I'm at last point and my desired point is ``basically'' 
      // where I'm at, then I'll say that j is the index of my correct point 
      // (indexInRange = true). If not then I know that desired point is 
      // between j and j + 1 ( where j + 1 not valid index ). Either way if
      // I'm running out of room, foundRange gets set to true since I've 
      // narrowed down my location to either valid or invalid.
      if ( b >  betan[j] and j+1 != betan.size() )      { ++j; }
      // If desired point is not to the right of current point, and also left of
      // (j - 1)th point, then I know I can just decrease my index. If I know
      // it's between j - 1 and j, then I have it narrowed down.
      // If I know that j == 1, and I know that desired point is not to the 
      // right of j - 1 (farthest left point, 0), then I know that my point 
      // has to be at 0. So for those latter two cases, I found the correct
      // index location and I have found the correct range.
      if ( b <= betan[j] and  b < betan[j-1] and j > 1 ){ --j; }
    }  
    while (
         (b >  betan[j] or ( b < betan[j-1] and j > 1 ) ) and 
         (b <= betan[j] or ( j+1 != betan.size() ) ) 
    );
    
    // If not in desired range
    if ( (b > betan[j] and ( j+1!=betan.size() or b >= 1.00001*betan[j]) ) ){
      sb[i] = 0; 
    }
    // If in desired range
    else { 
      floatT current = s[j]   < 0 ? slim : log( s[j]   );
      floatT toLeft  = s[j-1] < 0 ? slim : log( s[j-1] );
      sb[i] = current + (b-betan[j])*(toLeft-current)/(betan[j-1]-betan[j]);
      if (bet > 0) { sb[i] = sb[i] - bet; }
      if ( sb[i] > slim ){ sb[i] = exp(sb[i]); }
    }
  }

}











