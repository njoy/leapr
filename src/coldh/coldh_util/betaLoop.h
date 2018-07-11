#include <iostream>
#include <vector>
#include <cmath>
#include "betaLoop_util/jprimeLoop.h"
#include "betaLoop_util/bt.h"


template <typename A, typename F, typename R1, typename R2>
auto betaLoop( const A& betan, const A& rdbex, const A& bex, const A& sex, 
  const F& alphaVal, const F& wt, const F& tbart, const F& x, const F& y, 
  const F& evenSumConst, const F& oddSumConst, int itemp, int nbx, const int a, 
  int ncold, bool free, R1& sym_sab, R2& sym_sab_2 ){

  unsigned int jjmax = 2 * betan.size() - 1; 
  int k;
  double pj, be, snlg, snlk, sn, total;


  auto jjR = ranges::view::iota((unsigned) 0,jjmax);
  auto kR  = jjR | ranges::view::transform( [s=betan.size()](size_t jj){
                     return jj < s ? s - jj : jj - s + 2; } );
  auto beR = ranges::view::zip(jjR,kR) 
           | ranges::view::transform( [s=betan.size(),betan](auto t){
               size_t jj = std::get<0>(t); int k = std::get<1>(t);
               return jj < s - 1 ? -betan[k-1] : betan[k-1]; } );

  auto jR = ranges::view::iota(0,3) 
          | ranges::view::transform([ncold](auto i){
              if ( ncold == 1 or ncold == 4 ){
                return 2*i + 1;
              } else { return 2*i; } } );

  auto pjR = jR | ranges::view::transform([x](auto j){ return bt(j,x); });

  auto jj_k_be = ranges::view::zip(jjR,kR,beR);
  auto j_pj = ranges::view::zip(jR,pjR);
  auto bothLoops = ranges::view::cartesian_product(jj_k_be,j_pj);

  auto snlg_snlk_R = bothLoops | ranges::view::transform(
      [&total,x,evenSumConst,oddSumConst,bex,rdbex,sex,betan,alphaVal,wt,tbart,y,nbx,free]
                (auto t){
                 auto jj_k_be_val = std::get<0>(t);
                 auto j_pj_val    = std::get<1>(t);
                 auto jj   = std::get<0>(jj_k_be_val);
                 auto k    = std::get<1>(jj_k_be_val);
                 double be = std::get<2>(jj_k_be_val); 

                 auto j    = std::get<0>(j_pj_val);
                 auto pj   = std::get<1>(j_pj_val);
                 if (j == 1 or j == 0){ total = 0.0; }
                 auto snlg = jPrime( total, j, be, x, evenSumConst, pj, jj, bex, 
                   rdbex, sex, betan, alphaVal, wt, tbart, y, nbx, false, free );

                 auto snlk = jPrime( total, j, be, x, oddSumConst, pj, jj, bex, 
                   rdbex, sex, betan, alphaVal, wt, tbart, y, nbx, true, free );

                 //return snlg;
                 return std::make_pair(snlg,snlk);
               } );

  auto snR = snlg_snlk_R 
           | ranges::view::chunk(3) 
           | ranges::view::transform([](auto range){ 
              return ranges::accumulate(range,0.0,[](auto l,auto r){
                       return l+std::get<0>(r)+std::get<1>(r);}); } );

  
  auto snR2 = snR 
            | ranges::view::slice(0,int(betan.size()));
 
  std::cout << snR << std::endl;

  auto snR2_2  = snR | ranges::view::slice(int(betan.size())-1,ranges::end);

  auto snR2_rev = ranges::view::iota(0,int(betan.size())) 
                | ranges::view::transform([snR2](auto i){ 
                    return snR2[snR2.size()-i-1]; } );

  int size_a = sym_sab.size()/betan.size();
  int size_b = betan.size();
  auto output_sym_sab_L = sym_sab | ranges::view::slice( 0, a*size_b );
  auto output_sym_sab_R = sym_sab 
                        | ranges::view::slice( (a+1)*size_b, ranges::end );
  auto output_sym_sab = ranges::view::concat(output_sym_sab_L,snR2_rev,output_sym_sab_R);
  auto output_sym_sab_2_L = sym_sab_2 | ranges::view::slice( 0, a*size_b );
  auto output_sym_sab_2_R = sym_sab_2 
                        | ranges::view::slice( (a+1)*size_b, ranges::end );
  auto output_sym_sab_2 = ranges::view::concat(output_sym_sab_2_L,snR2_2,output_sym_sab_2_R);



  std::cout << output_sym_sab << std::endl;
  std::cout << output_sym_sab_2 << std::endl;

//if (jj < betan.size())    sym_sab(a,k-1,itemp) = sn;
//if (jj >= betan.size()-1) sym_sab_2(a,k-1,itemp) = sn;



  int counter = 0;

  for ( size_t jj = 0; jj < jjmax; ++jj ){

    k  = jj < betan.size()     ?  betan.size() - jj : jj - betan.size() + 2;
    be = jj < betan.size() - 1 ? -betan[k-1]        : betan[k-1];

    int start = ncold == 1 or ncold == 4 ? 1 : 0;
    int end   = ncold == 1 or ncold == 4 ? 6 : 5;

    total = 0; sn = 0;

    for ( auto j = start; j < end; j = j + 2 ){
      pj = bt(j,x);

      snlg = jPrime( total, j, be, x, evenSumConst, pj, jj, bex, rdbex, sex, betan, 
        alphaVal, wt, tbart, y, nbx, false, free );
      snlk = jPrime( total, j, be, x, oddSumConst, pj, jj, bex, rdbex, sex, betan, 
        alphaVal, wt, tbart, y, nbx, true,  free );
      auto rangeBoi = snlg_snlk_R[counter];
      //std::cout << snlg << "   " << snlk << std::endl;
      //std::cout << snlg << "   " << std::get<0>(rangeBoi) << "   " << snlk << "   " << std::get<1>(rangeBoi) << std::endl;

      ++counter;
      sn = sn + snlg + snlk;
    }
    //std::cout << sn << std::endl;
    //std::cout << "   " << std::endl;
      //if (jj < betan.size())    sym_sab(a,k-1,itemp) = sn;
      //if (jj >= betan.size()-1) sym_sab_2(a,k-1,itemp) = sn;
  }
}
