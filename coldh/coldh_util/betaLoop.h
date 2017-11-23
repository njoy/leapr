#include <iostream>
#include <vector>

auto betaLoop(int nbeta, std::vector<double>& betan ){
  //--loop over all beta values
  //    results for positive beta go into ssp
  //    results for negative beta go into sym_sab
  int jjmax = 2 * nbeta - 1;
  int k;
  for ( auto jj = 0; jj < jjmax; ++jj ){
    if (jj < nbeta)  k = nbeta - jj + 1;
    if (jj >= nbeta) k = jj - nbeta + 1;
    double be = betan[k-1];
    if (jj < nbeta) be = -be;
    int sn=0;
    int total=0;

    //--loop over all oscillators
    // para-h2: j=0,2,....; ortho-h2: j=1,3,....
    // ortho-d2: j=0,2,....; para-d2: j=1,3,....
    int ipo=1;
    if (law == 2.or.law == 5) ipo=2;
    int jt1 = 2 * 3; // THIS IS REALLY WEIRD originally 2 * jterm but jterm
                     // always equals 3? Please check out
    if (ipo == 2) jt1 = jt1 + 1;
    for ( auto l = ipo; l < jt1; l = l + 2; ){
      j=l-1;

   //   call bt(j,pj,x)

      //--sum over even values of j-prime
      snlg=0;
      do lp=1,10,2
        int jp=lp-1;
        auto betap=(-j*(j+1)+jp*(jp+1))*x/2;
        int  tmp=(2*jp+1)*pj*swe*4*sumh(j,jp,y);
        if (jj == 1.and.tmp >= small) {
          total=total+tmp;
        }
        bn=be+betap;
        if (ifree == 1) {
          ex=-(alp-abs(bn))**2/(4*alp);
          if (bn > zero) ex=ex-bn;
          add=exp(ex)/sqrt(4*pi*alp);
        } else {
          add=sint(bn,bex,rdbex,sex,nbx,al,wt,tbart, betan,nbeta,maxbb);
        }
        snlg=snlg+tmp*add;
      }

      //--sum over the odd values of j-prime
      snlk=0;
      for ( auto lp = 2; lp < 10; lp = lp + 2; ){
        jp=lp-1;
        betap=(-j*(j+1)+jp*(jp+1))*x/2;
        tmp=(2*jp+1)*pj*swo*4*sumh(j,jp,y);
        if (jj == 1.and.tmp >= small) {
          total=total+tmp;
        }
        bn = be + betap;
        if (ifree == 1) {
          ex=-(alp-abs(bn))**2/(4*alp);
          if (bn > zero) ex=ex-bn;
          add=exp(ex)/sqrt(4*pi*alp);
        else
          add=sint(bn,bex,rdbex,sex,nbx,al,wt,tbart,betan,nbeta,maxbb);
        }
        snlk = snlk + tmp * add;
      }

      //--continue the j loop
      sn=sn+snlg+snlk;
    }

    //--continue the beta loop
    if (jj <= nbeta) sym_sab[nal][k][itemp] = sn;
    if (jj >= nbeta) ssp[nal][k][itemp] = sn;
  }

}
