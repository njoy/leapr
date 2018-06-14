
auto upstk( const double& e, const std::vector<double>& s, 
    std::vector<std::vector<double>>& stk, int nl, int nx, int& i ){
 /*-------------------------------------------------------------------
  * Update the linearization stack with energy e and cross
  * sections s.  Here, i is the current index to the stack in stk,
  * nl is the number of legendre orders in s, and nx is the
  * cycle length in the stack.
  *-------------------------------------------------------------------
  */

  stk[0][i]=stk[0][i-1];
  stk[0][i-1]=e;
  for ( int j = 0; j < nl; ++ j ){
    stk[1+j][i]=stk[1+j][i-1];
    stk[1+j][i-1]=s[j];
  }
  i = i + 1;
}


