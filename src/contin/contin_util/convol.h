
#include <range/v3/all.hpp>

template <typename arrayT, typename floatT>
auto convol( const arrayT& t1, const arrayT& t2, const floatT& delta ){

  auto ijRange = ranges::view::cartesian_product(
                 ranges::view::iota(0,int(t2.size())),
                 ranges::view::iota(0,int(t1.size())));

  auto i1Range = ijRange 
               | ranges::view::transform([](auto t){
                   return std::get<0>(t) + std::get<1>(t);
                 });
  auto i2Range = ijRange
               | ranges::view::transform([](auto t){
                   return std::get<0>(t) - std::get<1>(t);
                 });

  auto t3Range = 
      ranges::view::zip(i1Range,i2Range,ijRange)
    | ranges::view::transform([t1,t2,delta](auto t){
        int i1 = std::get<0>(t), i2 = std::get<1>(t),
            j  = std::get<1>(std::get<2>(t));
        floatT f1 = ( i1-1 > int(t1.size())) ? 0 : t2[i1]*exp(-j*delta), f2 = 0;
        if (i2 >= 0 and  i2-1 < int(t1.size())){ f2 = t2[ i2];                 }
        if (i2 <  0 and -i2-3 < int(t1.size())){ f2 = t2[-i2] * exp(i2*delta); }

        return ( j == 0 or j == int(t1.size())-1 ) ? t1[j] * (f1+f2) * 0.5 :
                                                     t1[j] * (f1+f2); }) 
    | ranges::view::chunk(int(t1.size())) 
    | ranges::view::transform([delta](auto range){
        floatT rangeSum = ranges::accumulate(range,0.0);
        return ( rangeSum * delta < 1e-30 ) ? 0 : rangeSum * delta; });
		   
  return t3Range;
                
}

