#include <iostream>

auto checkMoments( const std::vector<double>& alpha, const std::vector<double>& beta, 
  std::vector<int>& maxt, const double& f0, const double& tbeta, 
  double tbar, std::vector<double>& ssm, std::ostream& output, int iprint ){

  double ff1, ff2, be, ssct, ex, al, alw;

  if (iprint == 2){
    output << std::endl << std::endl;
    output << "Checking S(a,b) with sum rule and normalization rule" << std::endl;
    output << "----------------------------------------------------" << std::endl;
    output << std::endl << std::endl;

    if (iprint != 0){
      for (int ibeta = 1; ibeta < beta.size(); ++ibeta){
        if (maxt[ibeta] > maxt[ibeta-1]){
            maxt[ibeta] = maxt[ibeta-1];
        }
      }
    }
  }

  // check the moments of s(alpha,beta)
  for ( int a = 0; a < int(alpha.size()); ++a ){
    al = alpha[a];
    alw = al*tbeta;
    double bel = 0, ff1l = 0, ff2l = 0, sum0 = 0, sum1 = 0;

    for ( int b = 0; b < int(beta.size()); ++b ){
      be = beta[b];//*sc;

      ex = -(alw-be)*(alw-be)/(4*alw*tbar);
      ssct = ex > -250.0 ? exp(ex)/sqrt(4*M_PI*alw*tbar) : 0;
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

    if (iprint == 2){
      output << std::setprecision(8)<< "   alpha    " << alpha[a] << "    (" << 
                                 a+1 << "/" << alpha.size() << ")" << std::endl;
      output << "             normalization check = " << sum0 << std::endl;
      output << "                  sum rule check = " << sum1 << std::endl << std::endl;
    }
  }

}

