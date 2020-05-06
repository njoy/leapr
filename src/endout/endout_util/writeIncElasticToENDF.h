#include "ENDFtk.hpp"



template <typename Range, typename Float>
auto writeIncElasticToENDF( Float sigma_b, Range temperatures, 
  Range debyeWallerVec ){//, int za, Float awr ){

  using namespace njoy::ENDFtk;
  using IncoherentElastic = section::Type<7,2>::IncoherentElastic;
  //using Elastic = section::Type<7,2>;

  if (temperatures.size() == 1){
      temperatures.resize(2);
      debyeWallerVec.resize(2);
      temperatures[1] = temperatures[0];
      debyeWallerVec[1] = debyeWallerVec[0];
  }

  std::vector<long> boundaries   = { long(temperatures.size()) },
                    interpolants = { 2 };
  IncoherentElastic incohElastic( sigma_b,
                                  std::move( boundaries ),
                                  std::move( interpolants ),
                                  std::move( temperatures ),
                                  std::move( debyeWallerVec ) ); 
  //Elastic chunk( za, awr, incohElastic );
  return incohElastic;
}
