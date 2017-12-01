#include <iostream>
#include <vector>
#include "../coldh/coldh_util/terpk.h"

auto skold( double cfrac, int itemp, double temp, 
  const std::vector<double>& alpha, const std::vector<double>& beta, 
  std::vector<double>& skappa,
  int ntempr, double awr, int lat, int nka, double dka, double scaling,
  std::vector<std::vector<std::vector<double>>>& symSab ){
  /* use skold approximation to add in the effects
   * of intermolecular coherence.
   */
  int i, j, k, kk, nal, ibeta;
  double tev, sc, amass, al, sk, ap, be, ss, s1, s2;
  double sum0, sum1, ff1l, ff2l, bel, ff1, ff2, waven;
  std::vector<double> scoh ( 1000, 0.0 );

  double amassn = 1.008664904; 
  double amu = 1.6605402e-24;
  double bk = 8.617385e-5;
  double hbar = 1.05457266e-27;
  double ev = 1.60217733e-12;

  
  // apply the skold approximation
  tev = bk * abs(temp);
  sc = 1;
  if (lat == 1) sc = 0.0253 / tev;
  amass = awr * amassn * amu;

  for ( auto b = 0; b < beta.size(); ++b ){
    for ( auto a = 0; a < alpha.size(); ++a ){
      al = alpha[b] * scaling;
      waven = 1.0e-8 * sqrt(2*amass*tev*ev*al)/hbar;
      sk = terpk(skappa,dka,waven);
      ap = alpha[b] / sk;
      for ( auto a2 = 0; a2 < alpha.size(); ++a2 ){
        kk = a2;
        if (ap < alpha[a2]) break;
      }
      if (kk == 0) kk = 1;
      call terp1(alpha(kk-1),ssm[kk-1][b][itemp],alpha(kk),ssm[kk][b][itemp],ap,scoh(a),5);
      scoh[a] = scoh[a] * sk;
    }
    for ( auto a = 0; a < nalpha; ++a ){
      ssm[a][b][itemp] = (1-cfrac)*ssm[a][b][itemp]+cfrac*scoh[a];
    }
  }
/*

  !--report the results
  for ( auto nal = 0; nal < nalpha; ++nal ){
    al=alpha(nal)*sc/arat
    for ( auto i = 0; i < beta.size(); ++i ){
      be=beta(i)*sc
      ss=ssm(i,nal,itemp)
      s1=ss*exp(-be/2)
      s2=ss*exp(-be)
    }
  }
      if (iprt == 1) then
         sum0=0
         sum1=0
         ff1l=0
         ff2l=0
         bel=0
         do ibeta=1,beta.size()
            be=beta(ibeta)
            ff2=ssm(ibeta,nal,itemp)
            ff1=ssm(ibeta,nal,itemp)*exp(-be)
            if (ibeta > 1) then
               sum0=sum0+(be-bel)*(ff1l+ff2l+ff1+ff2)/2
               sum1=sum1+(be-bel)*(ff2l*bel+ff2*be-ff1l*bel-ff1*be)/2
               ff1l=ff1
               ff2l=ff2
               bel=be
            else
               bel=be
               ff1l=ff1
               ff2l=ff2
               sum0=0
               sum1=0
            endif
         enddo
         sum1=sum1/al
         if (iprint == 2) then
            write(nsyso,'(''     normalization check ='',f8.4)') sum0
            write(nsyso,'(''          sum rule check ='',f8.4)') sum1
         else if (iprint == 1) then
            write(nsyso,'(1x,f10.4,2f10.4)') al,sum0,sum1
         endif
      endif
   enddo
   return
   end subroutine skold
*/
}

