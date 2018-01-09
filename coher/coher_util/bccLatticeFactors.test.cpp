#define CATCH_CONFIG_MAIN
#include <iostream>
#include "../../catch.hpp"
#include "bccLatticeFactors.h"

void equal( double a, double b ){
  if (b == 0.0){ 
    REQUIRE( std::abs(b-a) < 1e-6 );
    return;
  }
  REQUIRE ( std::abs( (a-b)/(b) ) < 1e-6 );
}

void equal_vec( std::vector<double> a, std::vector<double> b ){
  REQUIRE( a.size() == b.size() );
  for ( int i = 0; i < a.size(); ++i ){
    equal( a[i], b[i] );
  }
}

void equalBcc( std::tuple<int,int,int,double>& a, double& b ){
  equal( taubcc(std::get<0>(a),std::get<1>(a),std::get<2>(a),std::get<3>(a) ), b );
}


TEST_CASE( "taubcc" ){
  GIVEN( "inputs" ){
    std::vector<std::tuple<int,int,int,double>> inputs
      { {0,0,0,2}, {1,0,0,2}, {0,1,0,2}, {0,0,1,2}, {1,1,0,2},
        {1,0,1,2}, {0,1,1,2}, {1,1,1,2}, {1,2,3,4}, {5,3,6,4},
        {8,7,9,.1} };
    std::vector<double> output { 0, 78.956835, 78.956835, 78.956835, 236.870506,
      236.870506, 236.870506, 473.741011, 3947.84176, 21002.518165, 1519.9191 };
    
    for ( auto i = 0; i < output.size(); ++i ){
      equalBcc( inputs[i], output[i] );
    }
    
  } // GIVEN
} // TEST CASE



TEST_CASE( "Function to Compute BCC Lattice Factors" ){
  double ulim = 9.6e19, wint = 0, t2 = 5.7e-6, 
    a = 2.8e-8, c1 = 1.5e15;
  int ifl = 1, lat = 6;
  std::vector<double> b (60000, 0.0);
  GIVEN( "inputs" ){
    THEN( "outputs" ){
      int imax = bccLatticeFactors( ulim, b, ifl, wint, t2, lat, a, c1 );
      std::vector<double> b_0_to_100 { 7.97850038E+19, 4.47815741E-10, 
        7.62981036E+19, 4.57934242E-10, 7.29294035E+19, 4.68391128E-10, 
        6.96789033E+19, 4.79191756E-10, 6.65466032E+19, 4.90339686E-10, 
        6.35325030E+19, 5.01836225E-10, 6.06366029E+19, 5.13679885E-10, 
        5.78589027E+19, 5.25865760E-10, 5.51994026E+19, 5.38384815E-10, 
        5.26581025E+19, 5.51223073E-10, 5.02350024E+19, 5.64360718E-10, 
        4.79301023E+19, 5.77771098E-10, 4.57434021E+19, 5.91419663E-10, 
        4.36749021E+19, 6.05262853E-10, 4.17246020E+19, 6.19246969E-10, 
        3.98925019E+19, 6.33307095E-10, 3.81786018E+19, 6.47366141E-10, 
        3.65829017E+19, 6.61334099E-10, 3.51054016E+19, 6.75107632E-10, 
        3.37461016E+19, 6.88570141E-10, 3.25050015E+19, 7.01592432E-10, 
        3.13821015E+19, 7.14034144E-10, 3.03774014E+19, 7.25746050E-10, 
        2.94909014E+19, 7.36573292E-10, 2.87226013E+19, 7.46359569E-10,
        2.80725013E+19, 7.54952166E-10, 2.75406013E+19, 7.62207612E-10, 
        2.71269013E+19, 7.67997661E-10, 2.68314012E+19, 7.72215143E-10, 
        2.66541012E+19, 7.74779229E-10, 2.65950012E+19, 7.75639617E-10, 
        7.62981036E+19, 4.57934242E-10, 7.28703035E+19, 4.68581029E-10, 
        6.95607033E+19, 4.79598713E-10, 6.63693031E+19, 4.90994200E-10, 
        6.32961030E+19, 5.02772488E-10, 6.03411029E+19, 5.14936135E-10, 
        5.75043027E+19, 5.27484643E-10, 5.47857026E+19, 5.40413728E-10, 
        5.21853025E+19, 5.53714490E-10, 4.97031023E+19, 5.67372448E-10, 
        4.73391022E+19, 5.81366472E-10, 4.50933021E+19, 5.95667590E-10, 
        4.29657020E+19, 6.10237700E-10, 4.09563019E+19, 6.25028215E-10, 
        3.90651018E+19, 6.39978686E-10, 3.72921017E+19, 6.55015475E-10, 
        3.56373017E+19, 6.70050575E-10, 3.41007016E+19, 6.84980698E-10, 
        3.26823015E+19, 6.99686790E-10 };

      for ( auto i = 0; i < 100; ++i ){ equal( b[i], b_0_to_100[i] ); }
      std::vector<double> b_1199_to_1300 { 1.01846471E-9, 1.42431006E+19, 
        1.05988249E-9, 1.31793006E+19, 1.10182801E-9, 1.22337005E+19, 
        1.14361822E-9, 1.14063005E+19, 1.18437051E-9, 1.06971005E+19, 
        1.22300139E-9, 1.01061004E+19, 1.25825361E-9, 9.63330046E+18, 
        1.28876115E-9, 9.27870044E+18, 1.31315627E-9, 9.04230043E+18, 
        1.33021096E-9, 8.92410042E+18, 1.33899133E-9, 8.92410042E+18, 
        1.33899133E-9, 9.04230043E+18, 1.33021096E-9, 9.27870044E+18, 
        1.31315627E-9, 9.63330046E+18, 1.28876115E-9, 1.01061004E+19, 
        1.25825361E-9, 1.06971005E+19, 1.22300139E-9, 1.14063005E+19, 
        1.18437051E-9, 1.22337005E+19, 1.14361822E-9, 1.31793006E+19, 
        1.10182801E-9, 1.42431006E+19, 1.05988249E-9, 3.25050015E+19, 
        7.01592432E-10, 3.02001014E+19, 7.27873302E-10, 2.80134013E+19, 
        7.55748109E-10, 2.59449012E+19, 7.85297075E-10, 2.39946011E+19, 
        8.16588432E-10, 2.21625010E+19, 8.49670630E-10, 2.04486009E+19, 
        8.84561821E-10, 1.88529009E+19, 9.21235995E-10, 1.73754008E+19, 
        9.59605161E-10, 1.60161007E+19, 9.99497230E-10, 1.47750007E+19, 
        1.04062974E-9, 1.36521006E+19, 1.08258062E-9, 1.26474006E+19, 
        1.12475868E-9, 1.17609005E+19, 1.16637902E-9, 1.09926005E+19, 
        1.20645122E-9, 1.03425004E+19, 1.24379044E-9, 9.81060047E+18, 
        1.27706262E-9, 9.39690045E+18, 1.30487129E-9, 9.10140043E+18, 
        1.32588506E-9, 8.92410042E+18, 1.33899133E-9, 8.86500042E+18, 
        1.34344722E-9, 8.92410042E+18, 1.33899133E-9, 9.10140043E+18, 
        1.32588506E-9, 9.39690045E+18, 1.30487129E-9, 9.81060047E+18, 
        1.27706262E-9, 1.03425004E+19, 1.24379044E-9, 1.09926005E+19, 
        1.20645122E-9, 1.17609005E+19, 1.16637902E-9, 1.26474006E+19, 
        1.12475868E-9, 1.36521006E+19, 1.08258062E-9 };
      for ( auto i = 0; i < 100; ++i ){ equal( b[1199+i], b_1199_to_1300[i] ); }

      std::vector<double> b_52799_to_52900 { 1.18437051E-9, 1.25883006E+19, 
        1.12739586E-9, 1.38885006E+19, 1.07332763E-9, 1.53069007E+19, 
        1.02238944E-9, 1.68435008E+19, 9.74639055E-10, 1.84983008E+19, 
        9.30023819E-10, 2.02713009E+19, 8.88421746E-10, 2.21625010E+19, 
        8.49670630E-10, 2.41719011E+19, 8.13588097E-10, 2.62995012E+19, 
        7.79984971E-10, 2.85453013E+19, 7.48673868E-10, 3.09093014E+19, 
        7.19474483E-10, 1.11699005E+19, 1.19683791E-9, 1.01652004E+19, 
        1.25459056E-9, 9.27870044E+18, 1.31315627E-9, 8.51040040E+18, 
        1.37115008E-9, 7.86030037E+18, 1.42672548E-9, 7.32840035E+18, 
        1.47759494E-9, 6.91470033E+18, 1.52115444E-9, 6.61920031E+18, 
        1.55473805E-9, 6.44190030E+18, 1.57598830E-9, 6.38280030E+18, 
        1.58326773E-9, 6.44190030E+18, 1.57598830E-9, 6.61920031E+18, 
        1.55473805E-9, 6.91470033E+18, 1.52115444E-9, 7.32840035E+18, 
        1.47759494E-9, 7.86030037E+18, 1.42672548E-9, 8.51040040E+18, 
        1.37115008E-9, 9.27870044E+18, 1.31315627E-9, 1.01652004E+19, 
        1.25459056E-9, 1.11699005E+19, 1.19683791E-9, 1.22928005E+19, 
        1.14086583E-9, 1.35339006E+19, 1.08729777E-9, 1.48932007E+19, 
        1.03649203E-9, 1.63707007E+19, 9.88613088E-10, 1.79664008E+19, 
        9.43690207E-10, 1.96803009E+19, 9.01662740E-10, 2.15124010E+19, 
        8.62413504E-10, 2.34627011E+19, 8.25792600E-10, 2.55312012E+19, 
        7.91633869E-10, 2.77179013E+19, 7.59765937E-10, 3.00228014E+19, 
        7.30019371E-10, 3.24459015E+19, 7.02231114E-10, 1.10517005E+19, 
        1.20322109E-9, 1.01061004E+19, 1.25825361E-9, 9.27870044E+18, 
        1.31315627E-9, 8.56950041E+18, 1.36641379E-9, 7.97850038E+18, 
        1.41611771E-9, 7.50570036E+18, 1.46003874E-9, 7.15110034E+18, 
        1.49580009E-9, 6.91470033E+18, 1.52115444E-9 };
      for ( auto i = 0; i < 100; ++i ){ equal( b[52799+i], b_52799_to_52900[i] ); }

      std::vector<double> b_59549_to_59590 { 6.33307095E-10, 4.17246020E+19, 
        6.19246969E-10, 4.36749021E+19, 6.05262853E-10, 4.57434021E+19, 
        5.91419663E-10, 4.79301023E+19, 5.77771098E-10, 5.02350024E+19, 
        5.64360718E-10, 5.26581025E+19, 5.51223073E-10, 5.51994026E+19, 
        5.38384815E-10, 5.78589027E+19, 5.25865760E-10, 6.06366029E+19, 
        5.13679885E-10, 6.35325030E+19, 5.01836225E-10, 6.65466032E+19, 
        4.90339686E-10, 6.96789033E+19, 4.79191756E-10, 7.29294035E+19, 
        4.68391128E-10, 7.62981036E+19, 4.57934242E-10, 7.97850038E+19, 
        4.47815741E-10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

      for ( auto i = 0; i < 50; ++i ){ equal( b[59549+i], b_59549_to_59590[i] ); }

      REQUIRE( imax == 29789 );

    } // THEN
  } // GIVEN
} // TEST CASE

