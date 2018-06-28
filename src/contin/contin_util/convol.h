
#include <range/v3/all.hpp>

template <typename A, typename B>
auto convol( const A& t1, const A& t2, const B& delta ){

  auto iRange = ranges::view::iota(0,18);
  auto jRange = ranges::view::iota(0,6);
  auto ijRange = ranges::view::cartesian_product(iRange,jRange);
  std::cout << std::get<0>(ijRange[5]) << "    " << std::get<1>(ijRange[5]) << std::endl;
  std::cout << std::get<0>(ijRange[6]) << "    " << std::get<1>(ijRange[6]) << std::endl;
  std::cout << std::get<0>(ijRange[7]) << "    " << std::get<1>(ijRange[7]) << std::endl;
  std::cout << std::endl;
  std::cout << std::get<0>(ijRange[11]) << "    " << std::get<1>(ijRange[11]) << std::endl;
  std::cout << std::get<0>(ijRange[12]) << "    " << std::get<1>(ijRange[12]) << std::endl;
  std::cout << std::get<0>(ijRange[13]) << "    " << std::get<1>(ijRange[13]) << std::endl;
  int i;
  RANGES_FOR( auto x, ijRange ){ std::cout << i++ << "    " << std::get<0>(x) << "  " << std::get<1>(x) << std::endl;}
  std::vector<double> t3(t2.size(),0.0);


  int i1, i2, len_t1 = int(t1.size()), len_t2 = int(t2.size());
  double f1, f2;



  auto i1Range = ijRange 
               | ranges::view::transform([](auto t){
                   return std::get<0>(t) + std::get<1>(t);
                 });
  auto i2Range = ijRange
               | ranges::view::transform([](auto t){
                   return std::get<0>(t) - std::get<1>(t);
                 });

  auto i12 = ranges::view::zip(i1Range,i2Range);

  
  /*
  auto f1Range = ranges::view::zip(i12,ijRange)
               | ranges::view::transform([t1,t2,delta](auto t){
		   auto i1 = std::get<0>(std::get<0>(t));
		   auto j  = std::get<1>(std::get<1>(t));
		   return ( i1 - 1 > int(t1.size())) ? 0 : t2[i1]*exp(-j*delta);});
  auto f2Range = ranges::view::zip(i12)
               | ranges::view::transform([t1,t2,delta](auto t){
		   //auto i1 = std::get<0>(std::get<0>(t));
		   auto i2 = std::get<1>(std::get<0>(t));
                   if      ( i2 >= 0 and  i2-1 < int(t1.size()) ){ return t2[ i2];                  }
                   else if ( i2 <  0 and -i2-3 < int(t1.size()) ){ return t2[-i2] * exp( i2*delta );}
                   else                                  { return 0.0;                        }
		   } );
   
   */

  
  auto t3Range = ranges::view::zip(i12,ijRange)
               | ranges::view::transform([t1,t2,delta](auto t){
		   auto i1 = std::get<0>(std::get<0>(t));
		   auto i2 = std::get<1>(std::get<0>(t));
		   auto j  = std::get<1>(std::get<1>(t));
		   auto f1 = ( i1 - 1 > int(t1.size())) ? 0 : t2[i1]*exp(-j*delta);
                   double f2;
                   if      ( i2 >= 0 and  i2-1 < int(t1.size()) ){ f2 = t2[ i2];                  }
                   else if ( i2 <  0 and -i2-3 < int(t1.size()) ){ f2 = t2[-i2] * exp( i2*delta );}
                   else                                  { f2 = 0.0;                        }
		   return ( j == 0 or j == int(t1.size())-1 ) ?
                     t1[j] * ( f1+f2 ) * 0.5 :
                     t1[j] * ( f1+f2 );
		 } );

		   
                
  /*
  RANGES_FOR( auto entry, ijRange){
    std::cout << std::get<0>(entry) << "  " << 
		 std::get<1>(entry) << std::endl;
  }
  */
	       /*
  RANGES_FOR( auto entry, f1Range){
    std::cout << std::get<0>(entry) << "  " << 
	         std::get<0>(std::get<1>(entry)) << "   " << 
		 std::get<0>(std::get<1>(entry)) << std::endl;
  }
  */




  
  int counter = 0;
  for ( int i = 0; i < len_t2; ++i ){    // i iterates through t3
    for ( int j = 0; j < len_t1; ++j ){  // j iterates through t1, t2
      i1 = i + j;
      i2 = i - j;
      if ( t1[j] > 0 ){
        // Convolution will only be significant if t1[j] is not zero
    
        f1 = ( i1 - 1 > len_t1 ) ? 0 : t2[i1]*exp(-j*delta);

        if      ( i2 >= 0 and  i2-1 < len_t1 ){ f2 = t2[ i2];                  }
        else if ( i2 <  0 and -i2-3 < len_t1 ){ f2 = t2[-i2] * exp( i2*delta );}
        else                                  { f2 = 0;                        }

        // If at one of the endpoints, only give half contribution
	//auto value = ( j == 0 or j == len_t1 - 1 ) ? t1[j] * ( f1 + f2 ) * 0.5 :
        //                                         t1[j] * ( f1 + f2 );


        t3[i] += ( j == 0 or j == len_t1 - 1 ) ? t1[j] * ( f1 + f2 ) * 0.5 :
                                                 t1[j] * ( f1 + f2 );
	//std::cout << "----- " << value << "     " << t3Range[counter] << 
	//"      " << i << " " << j << "       " << 
	//std::get<0>(ijRange[counter]) << " " << 
	//std::get<1>(ijRange[counter])<<  std::endl;
	//std::cout << "----- " << i << " " << j << "       " << 
	//std::get<0>(ijRange[counter]) << " " << 
	//std::get<1>(ijRange[counter])<<  std::endl;
	++counter;
      } // if

    } // for

    t3[i] = ( t3[i] * delta < 1e-30 ) ? 0 : t3[i] * delta;

  } // for
  
  return t3;
}


