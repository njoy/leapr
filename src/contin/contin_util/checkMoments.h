#include <iostream>
//#include <unsupported/Eigen/CXX11/Tensor>
#include <range/v3/all.hpp>

template <typename rangeT>
auto checkMoments( const double& sc, const std::vector<double>& alpha,
  const std::vector<double>& beta, const std::vector<int>& maxt,
  /*int itemp,*/ const double& f0, const double& tbeta, const double& arat, 
  double tbar, rangeT ssm ){
  
  double /*ff0,*/ ff1, ff2, be, ssct, ex, al, alw;


  auto alphaBeta = ranges::view::cartesian_product( 
                     ranges::view::iota(0,int(alpha.size())), 
                     ranges::view::iota(0,int(beta.size())));
  
  /*
		     */

	  
  //int numTemps = 3;
  auto sab2 = ranges::view::repeat_n(0,int(alpha.size()*beta.size()))

//	    | ranges::view::chunk(numTemps)
//          | ranges::view::chunk(1)          Making it so that I only
//                                            have S(a,b) not S(a,b,T).
//                                            This will probs make ranges
//                                            easier. I'll just make a new
//                                            one of these for each iteration
//                                            through my happy temp loop
             | ranges::view::chunk(beta.size());

  /*
  auto sab3 = ranges::view::zip(
                ranges::view::repeat_n(0,int(alpha.size()*beta.size())),
                alphaBeta) ;
	*/
	 // | ranges::view::chunk(beta.size());
  

 
 
  auto al2 = alpha | ranges::view::transform([sc,iarat=1.0/arat](auto x){
                       return x*sc*iarat; });
  auto be2 = beta  | ranges::view::transform([sc]   (auto x){return x*sc;    });
  auto alw2= al2   | ranges::view::transform([tbeta](auto x){return x*tbeta; });
                   //| ranges::view::transform([beta](auto x){ return ranges::view::iota(0,int(beta.size())) 
                   //| ranges::view::transform([x](auto ){return x;}); }) 
                   //| ranges::view::join;

  //auto alw3 = ranges::view::iota(0,int(beta.size()*alpha.size())) | ranges::view::transform([alw2](auto x){ return 1.0*x; } );

  auto ex2a = ranges::view::cartesian_product(alw2,be2);
  auto ex2 = ex2a
             | ranges::view::transform( [tbar](auto x){ 
                 double alw = std::get<0>(x), 
                         be = std::get<1>(x);
                 return -std::pow(alw-be,2) / (4*alw*tbar ); });

  //int counter = 0;
  auto ssct2 = ranges::view::zip_with(
                 [tbar](auto ex,auto alw_be_tuple){
                   auto alw = std::get<0>(alw_be_tuple);
                   return ex>-250.0 ? exp(ex)/sqrt(4*M_PI*alw*tbar) : 0.0; }
		  ,ex2, ex2a);


  auto ff2_2 = ranges::view::zip_with( [maxt](auto ssct,auto ab_tuple,auto ssmVal){
    if (int(std::get<0>(ab_tuple))+1 >= maxt[int(std::get<1>(ab_tuple))]) { 
      return ssct; }
    return ssmVal;  }, ssct2, alphaBeta, ssm );


  auto ff1_2 = ranges::view::zip( ff2_2, alphaBeta )
             | ranges::view::transform([sc,beta](auto t){
                 auto ff2 = std::get<0>(t);
                 auto ab = std::get<1>(t);
                 auto be = beta[std::get<1>(ab)]*sc;
                 //return be;});//*exp(-be); });
                 return ff2*exp(-be); });
//RANGES_FOR( auto t, ff1_2 ){ std::cout << t <<std::endl; }

std::cout << std::endl;
  auto be_bel = ranges::view::concat(ranges::view::single(0.0),be2 )
              | ranges::view::sliding(2);
  auto ff1l_2 = ranges::view::concat(ranges::view::single(0.0),ff1_2 )
              | ranges::view::sliding(2);

  auto ff2l_2 = ranges::view::concat(ranges::view::single(0.0),ff2_2 )
              | ranges::view::sliding(2);

  auto sum0_2 = ranges::view::zip_with(                  
		  [alpha,beta,sc,f0,arat](auto ab, auto ff1l, auto ff2l){ 
                    auto bel = beta[std::get<1>(ab)-1]*sc;
                    auto be  = beta[std::get<1>(ab)]*sc;
                    auto al  = alpha[std::get<0>(ab)]*sc/arat;
                    if (std::get<1>(ab) == 0){ return 0.0; }
                    return (be-bel)*(ff1l[0]+ff1l[1]+ff2l[0]+ff2l[1])*0.5/(1.0-exp(-al*f0));},
                    //return (1.0-exp(-al*f0));
      //sum0 = sum0/(1-exp(-al*f0));
                  alphaBeta,ff1l_2,ff2l_2);
  auto sum0_2b = ranges::view::cartesian_product(ff1l_2,be_bel);
//RANGES_FOR( auto t, sum0_2b ){ std::cout << std::get<1>(t)[0] << "   " << std::get<1>(t)[1]  << "    " << std::get<0>(t)[0] << "    " << std::get<0>(t)[1]  <<std::endl; }
		  //[](auto be, auto ff1l, auto ff2l){ 
                  //  return (be[1]-be[0])*(ff1l[0]+ff1l[1]+ff2l[0]+ff2l[1])*0.5;},
                  //be_bel,ff1l_2,ff2l_2);



  auto sum1_2 = ranges::view::zip_with(                  
		  [](auto be, auto ff1l, auto ff2l){ 
                    return (be[1]-be[0])*(-ff1l[0]*be[0]-ff1l[1]*be[1]+
                                           ff2l[0]*be[0]+ff2l[1]*be[1]);},
                  be_bel,ff1l_2,ff2l_2);

  //RANGES_FOR( auto b, ff1l_2 ){ std::cout << b[0] << "     " << b[1] << std::endl;}

double countDouble = 0.0;
RANGES_FOR( auto x, sum0_2 ){ 
  countDouble += x;
}

/*
RANGES_FOR( auto x, ranges::view::zip( ff1l_2,ff2l_2 )){
    auto ff1l_val = std::get<0>(x);
    auto ff2l_val = std::get<1>(x);
    //std::cout << ff2l_val << std::endl;
    std::cout << ff1l_val[0] << "  " << ff1l_val[1] << "  " << ff2l_val[0] << "   " << ff2l_val[1] << std::endl;
  }
  */

  {
  std::vector<std::vector<double>> ssm (alpha.size(),std::vector<double>(beta.size(),0.0));


  std::cout << std::endl;
  std::cout << ff2_2 << std::endl;
  std::cout << std::endl;

  // this definition is dumb
  int naint = 1;

  // int nbint = 1;
  // check the moments of s(alpha,beta)
  double sum0 = 0, sum1 = 0;
  int counter = 0; 
  for ( size_t a = 0; a < alpha.size(); ++a ){
    if ( ( a % naint == 0 ) or ( a == alpha.size() - 1) ){
      al = alpha[a]*sc/arat;
      alw = al*tbeta;
      double bel = 0, ff1l = 0, ff2l = 0;// sum0 = 0, sum1 = 0;

      for ( size_t b = 0; b < beta.size(); ++b ){
        //int jprt=(b)%nbint+1;               // This doesn't seem to do
        //if (b == beta.size()-1) jprt=1;     // anything
        be = beta[b]*sc;
        ex = -(alw-be)*(alw-be)/(4*alw*tbar);
        ssct = ex > -250.0 ? exp(ex)/sqrt(4*M_PI*alw*tbar) : 0;
        if (int(a)+1 >= maxt[b]) {
          ssm[a][b] = ssct;
        }
        ff2 = ssm[a][b];
        ff1 = ssm[a][b]*exp(-be);
        //std::cout << ff1l << "      " << ff1  << std::endl;
        //std::cout << bel << "   " << be << "    " << ff1l << "   " << ff1 << std::endl;//"   " << ff2l << "   " << ff2 << std::endl;
        //std::cout << ff1 << "     " << ff1_2[counter++] << std::endl;
        //ff0 = ssm[a][b][itemp]*exp(-be/2);   // This isn't used either
        if (b > 0) {
          sum0 = sum0+(be-bel)*(ff1l+ff2l+ff1+ff2)/2.0;
          sum1 = sum1+(be-bel)*(ff2l*bel+ff2*be-ff1l*bel-ff1*be)/2.0;
        }
        else {
          sum0 += 0;
          sum1 += 0;
        }
        std::cout << sum0 << std::endl;
      sum0 = sum0/(1.0-exp(-al*f0));
        std::cout << sum0 << "     " << ((be-bel)*(ff1l+ff2l+ff1+ff2)/2.0)/(1.0-exp(-al*f0))<< "       " << sum0_2[counter++] << std::endl;
        //std::cout << ((be-bel)*(ff1l+ff2l+ff1+ff2)/2.0)/(1-exp(-al*f0))<< "       " << sum0_2[counter++] << std::endl;
        ff1l = ff1;
        ff2l = ff2;
        bel = be;
      }
      //std::cout << "sum1:  " << sum1 << std::endl;
      
      //sum0 = sum0/(1.0-exp(-al*f0));
      //std::cout << "sum0:  " << (1-exp(-al*f0)) << "     alpha:  " << a << std::endl;
      //sum1 = sum1/al/tbeta;
      //sum0 = 0; sum1= 0;
    }
  }
  std::cout << sum0 << std::endl;
  ++counter;
  }
 
  std::cout << std::endl;
  std::cout << sum0_2 << std::endl;
  std::cout << std::endl;
auto sum = ranges::accumulate(sum0_2, 0.0);
  std::cout << sum << std::endl;
  return ff2_2;
  RANGES_FOR( auto x, ranges::view::zip(alpha,sum0_2)){ std::cout << std::get<1>(x) << "      " << std::get<0>(x) << std::endl;
      }

  std::cout << sum1_2 << std::endl;

  std::cout << f0 << std::endl;


  std::cout << sab2 << std::endl;
  //std::cout << al2 << std::endl;
  std::cout << std::endl;
  std::cout << be2 << std::endl;
  std::cout << std::endl;
  std::cout << alw2 << std::endl;
  std::cout << std::endl;
  //std::cout << ( alpha | ranges::view::all )<< std::endl;
  //std::cout << std::endl;
  //std::cout << ssct2 << std::endl;


  RANGES_FOR(auto entry , alphaBeta){
    std::cout <<  std::get<0>(entry) << "      " << 
	          std::get<1>(entry) << std::endl;
  }   
return ff2_2;
  /*
  */

  /*
  RANGES_FOR(auto entry , sab3){
    std::cout << std::get<0>(entry) << "    " << 
                 std::get<0>(std::get<1>(entry)) << "  " << 
                 std::get<1>(std::get<1>(entry)) << std::endl;
  }   
  */

  



}

