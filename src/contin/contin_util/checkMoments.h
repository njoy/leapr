#include <iostream>
#include <unsupported/Eigen/CXX11/Tensor>
#include <range/v3/all.hpp>

auto checkMoments( const double& sc, const std::vector<double>& alpha,
  const std::vector<double>& beta, const std::vector<int>& maxt,
  int itemp, const double& f0, const double& tbeta, const double& arat, 
  double tbar, Eigen::Tensor<double,3>& ssm ){

  double /*ff0,*/ ff1, ff2, be, ssct, ex, al, alw;


  auto alphaBeta = ranges::view::cartesian_product( 
                     ranges::view::iota(0,int(alpha.size())), 
                     ranges::view::iota(0,int(beta.size())));

  auto sab1 = ranges::view::repeat_n(0,int(alpha.size()*beta.size()));
	  
  auto sab2 = ranges::view::repeat_n(0,int(alpha.size()*beta.size()))
            | ranges::view::chunk(1)   

//          | ranges::view::chunk(1)          Making it so that I only
//                                            have S(a,b) not S(a,b,T).
//                                            This will probs make ranges
//                                            easier. I'll just make a new
//                                            one of these for each iteration
//                                            through my happy temp loop
             | ranges::view::chunk(beta.size());

  auto sab3 = ranges::view::zip(
                ranges::view::repeat_n(0,int(alpha.size()*beta.size())),
                alphaBeta) ;
	 // | ranges::view::chunk(beta.size());
  
  
  auto al2 =alpha|ranges::view::transform([sc,iarat=1.0/arat](auto x){return x*scr*iarat;});
  auto be2 = beta|ranges::view::transform([sc]     (auto x){return x*sc;     });
  auto alw2= al2 |ranges::view::transform([tbeta]  (auto x){return x*tbeta;  });


  auto ssct2 = ranges::view::zip(
                 ranges::view::for_each(alw2,[beta](double d){  
                   return ranges::yield_from(
                     ranges::view::repeat_n(d,beta.size())); }),
                 be2 | ranges::view::cycle 
                     | ranges::view::take_exactly(alpha.size()*beta.size()))
             | ranges::view::transform( [tbar](auto x){ 
                 double alw = std::get<0>(x), be = std::get<1>(x);
                 double ex = -std::pow(alw-be,2) / (4*alw*tbar );
                 return ex>-250.0 ? exp(ex)/sqrt(4*M_PI*alw*tbar) : 0.0; });


  // ex2 = zip [ a1, a1, a1, ... a1, a2, a2, ... , ... , aN ]
  //                               with
  //           [ b1, b2, b3, ... bN, b1, b2, ... , ... , bN ]
  //  to get [ (a1,b1), (a1,b2), (a1,b3), ... (a1,bN), (a2,b1), ... , (aN,bN) ]
  //  and then plus that into ex = -(alw-be)^2 / (4*alw*tbar)
  //  where ai represents alw = a*sc*tbeta/arat and bi = b*sc
  //  and then int ssct form where 
  //      ssct = ex > -250.0 ? exp(ex)/sqrt(4*M_PI*alw*tbar) : 0;


  std::cout << sab1 << std::endl;  
  std::cout << sab2 << std::endl;  
//  std::cout << sab3 << std::endl;  
  //std::cout << ssct2[1] << std::endl;  
  std::cout << ssct2 << std::endl;  





  // this definition is dumb
  int naint = 1;

  // int nbint = 1;
  // check the moments of s(alpha,beta)
  //int counter = 0; 
  for ( size_t a = 0; a < alpha.size(); ++a ){
    if ( ( a % naint == 0 ) or ( a == alpha.size() - 1) ){
      al = alpha[a]*sc/arat;
      alw = al*tbeta;
      double bel = 0, ff1l = 0, ff2l = 0, sum0 = 0, sum1 = 0;

      for ( size_t b = 0; b < beta.size(); ++b ){
        //int jprt=(b)%nbint+1;               // This doesn't seem to do
        //if (b == beta.size()-1) jprt=1;     // anything
        be = beta[b]*sc;
        ex = -(alw-be)*(alw-be)/(4*alw*tbar);
        ssct = ex > -250.0 ? exp(ex)/sqrt(4*M_PI*alw*tbar) : 0;
	std::cout << ssct << "     " <<  std::endl;
        if (int(a)+1 >= maxt[b]) {
          ssm(a,b,itemp) = ssct;
        }
        ff2 = ssm(a,b,itemp);
        ff1 = ssm(a,b,itemp)*exp(-be);
        //ff0 = ssm[a][b][itemp]*exp(-be/2);   // This isn't used either
        if (b > 0) {
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
  return;



  std::cout << sab2 << std::endl;
  std::cout << al2 << std::endl;
  std::cout << std::endl;
  std::cout << be2 << std::endl;
  std::cout << std::endl;
  std::cout << alw2 << std::endl;
  std::cout << std::endl;
  //std::cout << ( alpha | ranges::view::all )<< std::endl;
  //std::cout << std::endl;
  std::cout << ssct2 << std::endl;


  RANGES_FOR(auto entry , alphaBeta){
    std::cout <<  std::get<0>(entry) << "      " << 
	          std::get<1>(entry) << std::endl;
  }   

  RANGES_FOR(auto entry , sab3){
    std::cout << std::get<0>(entry) << "    " << 
                 std::get<0>(std::get<1>(entry)) << "  " << 
                 std::get<1>(std::get<1>(entry)) << std::endl;
  }   




}

