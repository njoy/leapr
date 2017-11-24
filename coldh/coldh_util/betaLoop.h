#include <iostream>
#include <vector>

auto betaLoop( std::vector<double>& betan, double alp, double x, int itemp,
    int a, int law, std::vector<std::vector<std::vector<double>>>& sym_sab ){
  //--loop over all beta values
  //    results for positive beta go into ssp
  //    results for negative beta go into sym_sab
  double pj, swo, swe;
  int ifree = 0;
  int jjmax = 2 * betan.size() - 1;
  int k, jp, j;
  double add, betap, be;
  double pi = 3.14159265358979;
  for ( auto jj = 0; jj < jjmax; ++jj ){
    if ( jj < betan.size() ){
      k = betan.size() - jj;
      be = -betan[k-1];
    } else {
      k = jj - betan.size() + 2;
      be = betan[k-1];
    }
    int sn = 0;
    int total = 0;

    //--loop over all oscillators
    // para-h2: j=0,2,....; ortho-h2: j=1,3,....
    // ortho-d2: j=0,2,....; para-d2: j=1,3,....
    int ipo = law == 2 or law == 5 ? 2 : 1;

    //int jt1 = 2 * 3; // THIS IS REALLY WEIRD origialy 2 * jterm but jterm
                     // always equals 3? Please check out
    int jt1 = law == 2 or law == 5 ? 7 : 6;

    for ( auto l = ipo; l < jt1; l = l + 2 ){
      j = l - 1;

   //   call bt(j,pj,x)

      //--sum over even values of j-prime
      double snlg = 0.0;
      for ( auto lp = 1; lp < 10; lp = lp + 2 ){
        jp = lp - 1;
        betap=(-j*(j+1)+jp*(jp+1))*x/2;
        int  tmp=(2*jp+1)*pj*swe*4;//*sumh(j,jp,y);
        if (jj == 1 and tmp >= 1.0e-6) {
          total += tmp;
        }
        double bn=be+betap;
        if (ifree == 1) {
          double ex=-(alp-abs(bn))*(alp-abs(bn))/(4*alp);
          if (bn > 0.0) ex=ex-bn;
          add = exp(ex) / sqrt(4*pi*alp);
        }// else {
         // add=sint(bn,bex,rdbex,sex,nbx,al,wt,tbart, betan,betan.size(),maxbb);
       // }
        snlg = snlg + tmp * add;
      }
/*
      //--sum over the odd values of j-prime
      double snlk=0;
      for ( auto lp = 2; lp < 10; lp = lp + 2 ){
        jp=lp-1;
        betap=(-j*(j+1)+jp*(jp+1))*x/2;
        int tmp=(2*jp+1)*pj*swo*4;//*sumh(j,jp,y);
        if (jj == 1 and tmp >= 1.0e-6) {
          total += tmp;
        }
        double bn = be + betap;
        if (ifree == 1) {
          double ex=-(alp-abs(bn))*(alp-abs(bn))/(4*alp);
          if (bn > 0.0) ex=ex-bn;
          add = exp(ex) / sqrt(4*pi*alp);
//        else
  ///        add=sint(bn,bex,rdbex,sex,nbx,al,wt,tbart,betan,betan.size(),maxbb);
        }
        snlk = snlk + tmp * add;
      }
u

    */
      //--continue the j loop
    //  sn=sn+snlg+snlk;
    }

    //--continue the beta loop
  //  std::cout << k << std::endl;
    //if (jj <= betan.size()) sym_sab[a][k-1][itemp] = sn;
  //  if (jj >= betan.size()) ssp[a][k][itemp] = sn;
  }

}
