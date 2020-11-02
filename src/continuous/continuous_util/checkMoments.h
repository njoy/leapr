#include <iostream>

template <typename F, typename A, typename A_of_ints>
auto checkMoments( const A& alpha, const A& beta, 
  const A_of_ints& maxt, const F& f0, const F& tbeta, 
  F tbar, A& ssm ){

  F ff1, ff2, be, ssct, ex, al, alw;

  // this definition is dumb
  int naint = 1;

  // int nbint = 1;
  // check the moments of s(alpha,beta)
  for ( int a = 0; a < int(alpha.size()); ++a ){
    if ( ( a % naint == 0 ) or ( a == int(alpha.size()) - 1) ){
      al = alpha[a];//sc/arat;
      alw = al*tbeta;
      F bel = 0, ff1l = 0, ff2l = 0, sum0 = 0, sum1 = 0;

      for ( int b = 0; b < int(beta.size()); ++b ){
        //int jprt=(b)%nbint+1;               // This doesn't seem to do
        //if (b == beta.size()-1) jprt=1;     // anything
        be = beta[b];//*sc;

        ex = -(alw-be)*(alw-be)/(4*alw*tbar);
        ssct = ex > -250.0 ? exp(ex)/sqrt(4*M_PI*alw*tbar) : 0;
        //if ( a==0 and b==0){ std::cout << "--  " <<  ssct << std::endl; }
        if (a+1 >= maxt[b]) {
          ssm[b+a*beta.size()] = ssct;
        }
        ff2 = ssm[b+a*beta.size()];
        ff1 = ssm[b+a*beta.size()]*exp(-be);
        if (b > 1) {
          sum0 = sum0+(be-bel)*(ff1l+ff2l+ff1+ff2)/2;
          sum1 = sum1+(be-bel)*(ff2l*bel+ff2*be-ff1l*bel-ff1*be)/2;
        }
        else {
          sum0 = 0;
          sum1 = 0;
        }
        ff1l = ff1;
        ff2l = ff2;
        bel = be;
      }
      sum0 = sum0/(1-exp(-al*f0));
      sum1 = sum1/al/tbeta;
    }
  }
}

