#include <iostream>
auto checkMoments( const double& sc, const std::vector<double>& alpha,
    const std::vector<double>& beta, const std::vector<int>& maxt,
    int itemp, double f0, double tbeta, double arat, double tbar, double explim,
    std::vector<std::vector<std::vector<double>>>& ssm ){

  double bel, ff0, ff1, ff1l, ff2, ff2l, sum0, sum1, be, alp, ssct, ex, al, alw;

  int naint = 1, nbint = 1;
  // check the moments of s(alpha,beta)
  for ( size_t a = 0; a < alpha.size(); ++a ){
      int iprt = ((a)%naint)+1;
      if (a == alpha.size()-1) iprt=1;
      if (iprt == 1) {
         al=alpha[a]*sc/arat;
         bel=0;
         ff1l=0;
         ff2l=0;
         sum0=0;
         sum1=0;
         for ( size_t b = 0; b < beta.size(); ++b ){
            int jprt=(b)%nbint+1;
            if (b == beta.size()-1) jprt=1;
            be=beta[b]*sc;
            alw=al*tbeta;
            alp=alw*tbar;
            ex=-(alw-be)*(alw-be)/(4*alp);
            ssct=0;
            if (ex > explim) ssct=exp(ex)/sqrt(4*M_PI*alp);
            //std::cout << a+1 << "   " << b+1 << "    " << maxt[b] << "    " << ssct << std::endl;
            if (int(a) >= maxt[b]-1) {
               ssm[a][b][itemp] = ssct;
            }
            ff2=ssm[a][b][itemp];
            ff1=ssm[a][b][itemp]*exp(-be);
            ff0=ssm[a][b][itemp]*exp(-be/2);
            if (b > 1) {
               sum0=sum0+(be-bel)*(ff1l+ff2l+ff1+ff2)/2;
               sum1=sum1+(be-bel)*(ff2l*bel+ff2*be-ff1l*bel-ff1*be)/2;
               ff1l=ff1;
               ff2l=ff2;
               bel=be;
            }
            else {
               bel=be;
               ff1l=ff1;
               ff2l=ff2;
               sum0=0;
               sum1=0;
            }
          }
         sum0=sum0/(1-exp(-al*f0));
         sum1=sum1/al/tbeta;
    }
  }
}

