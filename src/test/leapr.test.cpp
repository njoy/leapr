#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "leapr.cpp"
#include <range/v3/all.hpp>
#include "generalTools/testing.h"
#include <variant>
#include "tsl_h2o_tape24.h"
#include "tsl_sic_tape24.h"
#include "tsl_be_tape24.h"
#include "tsl_ortho_h_tape24.h"


using Tabulated = section::Type< 7, 4 >::Tabulated;
using CoherentElastic = section::Type< 7, 2 >::CoherentElastic;
using IncoherentElastic = section::Type< 7, 2 >::IncoherentElastic;







template <typename MF7>
void check_MF7(MF7 my_MF7, MF7 good_MF7){ 
      REQUIRE( my_MF7.hasMT(4) == good_MF7.hasMT(4));
      REQUIRE( my_MF7.hasMT(2) == good_MF7.hasMT(2));

      if (my_MF7.hasMT(2)){
        njoy::ENDFtk::section::Type<7,2>   my_MT2 = my_MF7.MT(2_c);
        njoy::ENDFtk::section::Type<7,2> good_MT2 = my_MF7.MT(2_c);
        REQUIRE( my_MT2.ZA() == Approx(good_MT2.ZA()) );
        REQUIRE( my_MT2.AWR() == Approx(good_MT2.AWR()) );
        REQUIRE( my_MT2.LTHR() == good_MT2.LTHR() );
        REQUIRE( my_MT2.elasticScatteringType() == good_MT2.elasticScatteringType() );

        if (my_MT2.elasticScatteringType() == 1){
          auto   my_law = std::get<CoherentElastic>(  my_MT2.scatteringLaw());
          auto good_law = std::get<CoherentElastic>(good_MT2.scatteringLaw());

          //REQUIRE( my_law.ZA() == Approx(good_law.ZA()) );
          //REQUIRE( my_law.AWR() == Approx(good_law.AWR()) );
          REQUIRE( my_law.LTHR() == good_law.LTHR() );
          REQUIRE( my_law.temperatureDependenceFlag() == good_law.temperatureDependenceFlag() );
          REQUIRE( my_law.numberTemperatures() == good_law.numberTemperatures() );
          REQUIRE( my_law.numberBraggEdges() == good_law.numberBraggEdges() );
          REQUIRE( my_law.LT() == good_law.LT() );
          REQUIRE( my_law.NP() == good_law.NP() );
          checkVec(my_law.boundaries(),good_law.boundaries());
          checkVec(my_law.interpolants(),good_law.interpolants());
          checkVec(my_law.temperatures(),good_law.temperatures());
          checkVec(my_law.energies(),good_law.energies());
          checkVec(my_law.LI(),good_law.LI());
          checkVec(my_law.temperatureInterpolants(),good_law.temperatureInterpolants());
          REQUIRE(my_law.thermalScatteringValues().size() == good_law.thermalScatteringValues().size());
          for (size_t j = 0; j < my_law.thermalScatteringValues().size(); ++j){
            checkVec(my_law.thermalScatteringValues()[j],good_law.thermalScatteringValues()[j]);
          }



        }
        else if (my_MT2.elasticScatteringType() == 2){
          auto   my_law = std::get<IncoherentElastic>(  my_MT2.scatteringLaw());
          auto good_law = std::get<IncoherentElastic>(good_MT2.scatteringLaw());
          REQUIRE( my_law.LTHR() == good_law.LTHR() );
          REQUIRE( my_law.elasticScatteringType() == good_law.elasticScatteringType() );
          REQUIRE( my_law.SB() == good_law.SB() );
          REQUIRE( my_law.SB() == good_law.SB() );
          REQUIRE( my_law.NP() == good_law.NP() );
          REQUIRE( my_law.numberTemperatures() == good_law.numberTemperatures() );
          REQUIRE( my_law.NR() == good_law.NR() );
          checkVec(my_law.interpolants(),good_law.interpolants());
          checkVec(my_law.boundaries(),good_law.boundaries());
          checkVec(my_law.temperatures(),good_law.temperatures());
          checkVec(my_law.debyeWallerValues(),good_law.debyeWallerValues());

        }

      } // if mt2


      if (my_MF7.hasMT(4)){

        njoy::ENDFtk::section::Type<7,4> my_MT4   =   my_MF7.MT(4_c);
        njoy::ENDFtk::section::Type<7,4> good_MT4 = good_MF7.MT(4_c);

        //std::string buffer;
        //auto output = std::back_inserter(buffer);
        //my_MT4.print(output,27,7);
        //std::cout << buffer << std::endl;
 
 
        REQUIRE( my_MT4.MT() == good_MT4.MT() );
        REQUIRE( my_MT4.ZA() == good_MT4.ZA() );
        REQUIRE( my_MT4.AWR() == good_MT4.AWR() );
        REQUIRE( my_MT4.LAT() == good_MT4.LAT() );
        REQUIRE( my_MT4.temperatureOption() == good_MT4.temperatureOption() );
        REQUIRE( my_MT4.LASYM() == good_MT4.LASYM() );
        REQUIRE( my_MT4.symmetryOption() == good_MT4.symmetryOption() );
        REQUIRE( my_MT4.NC() == good_MT4.NC() );
  
        auto   my_barray =   my_MT4.constants(),
             good_barray = good_MT4.constants();
        REQUIRE( my_barray.LLN()             == good_barray.LLN()             );
        REQUIRE( my_barray.sabStorageType()  == good_barray.sabStorageType()  );
        REQUIRE( my_barray.NI()              == good_barray.NI()              );
        REQUIRE( my_barray.numberConstants() == good_barray.numberConstants() );
        REQUIRE( my_barray.NS()              == good_barray.NS()              );
        REQUIRE(   my_barray.numberNonPrincipalScatterers() == 
                   good_barray.numberNonPrincipalScatterers() );

        REQUIRE( my_barray.epsilon() == Approx(good_barray.epsilon()).epsilon(1e-6) );
        REQUIRE(   my_barray.upperEnergyLimit() == 
          Approx(good_barray.upperEnergyLimit()).epsilon(1e-6) );
        checkVec(  my_barray.totalFreeCrossSections(),
                 good_barray.totalFreeCrossSections());
        checkVec(  my_barray.atomicWeightRatios(),
                 good_barray.atomicWeightRatios());
        checkVec(  my_barray.numberAtoms(),
                 good_barray.numberAtoms());
        checkVec(  my_barray.analyticalFunctionTypes(),
                 good_barray.analyticalFunctionTypes());
   
        auto   my_table = std::get<Tabulated>(  my_MT4.scatteringLaw()),
             good_table = std::get<Tabulated>(good_MT4.scatteringLaw());
  
        REQUIRE( my_table.NR() == good_table.NR() );
        REQUIRE( my_table.NB() == good_table.NB() );
        REQUIRE( my_table.numberBetas() == good_table.numberBetas() );

        checkVec(my_table.boundaries(),good_table.boundaries());
        checkVec(my_table.interpolants(),good_table.interpolants());
  
        for (int ibeta = 0; ibeta < my_table.numberBetas(); ++ibeta){
          auto my_value   =   my_table.betas()[ibeta];
          auto good_value = good_table.betas()[ibeta];
          REQUIRE( my_value.beta() == Approx(good_value.beta()).epsilon(1e-6) );
          REQUIRE( my_value.LT() == good_value.LT() );
          REQUIRE(   my_value.temperatureDependenceFlag() == 
                   good_value.temperatureDependenceFlag() );
          REQUIRE( my_value.NT() == good_value.NT() );
          REQUIRE( my_value.numberTemperatures() == good_value.numberTemperatures() );
  
          REQUIRE( my_value.NR() == good_value.NR() );
          REQUIRE( my_value.NA() == good_value.NA() );
          REQUIRE( my_value.numberAlphas() == good_value.numberAlphas() );
          checkVec( my_value.boundaries(),good_value.boundaries());
          checkVec( my_value.interpolants(),good_value.interpolants());
    
          checkVec( my_value.temperatures(),good_value.temperatures());
          checkVec( my_value.alphas(),good_value.alphas());
          checkVec( my_value.temperatureInterpolants(),good_value.temperatureInterpolants());
          REQUIRE(   my_value.thermalScatteringValues().size() == 
                   good_value.thermalScatteringValues().size());
          for ( size_t i = 0; i < my_value.thermalScatteringValues().size(); ++i){
            checkVec(  my_value.thermalScatteringValues()[i],
                     good_value.thermalScatteringValues()[i]);
          }
        } // beta loop
      } // if mf7 is present

}











TEST_CASE("full func"){
  int nphon, ncold, lat, iel, npr, nss, b7, mat, isabt, ilog, nsk, mss, za;
  double sps, awr, aws, delta, twt, c, tbeta, smin, spr;
  std::vector<double> alphas, betas, temps;

  /*


  GIVEN( "H in H2O ENDF-B/VIII.0 input" ) {
    WHEN( "continuous, translational, and discrete oscillator options used" ) {
      nphon = 100;
      mat = 1; za = 1001; isabt = 0; ilog = 0; smin = 2e-38;
      awr = 0.9991673; spr = 20.43608; npr = 2; iel = 0; ncold = 0; nsk = 0;
      nss = 1; b7 = 1; aws = 15.85751; sps = 3.7939; mss = 1;
      lat = 1;
      alphas = {1e-5, 12.0};
      betas  = {0.0, 1.0, 9.0};
      std::vector<int> generalInfo {nphon,mat,za,isabt,ilog,lat};
      std::vector<int> scatterControl {npr,iel,ncold,nsk,nss,b7,mss};
      std::vector<double> scatterInfo {awr,spr,aws,sps};


      temps = { 283.6, 550.0 };
      std::vector<double> rho_283 
      { 0.000000E+00, 4.005204E-01, 1.119462E+00, 1.698629E+00, 2.010136E+00, 
        2.072919E+00, 1.898872E+00, 1.580587E+00, 1.322490E+00, 1.206897E+00, 
        1.171350E+00, 1.193990E+00, 1.252900E+00, 1.308722E+00, 1.372691E+00, 
        1.458393E+00, 1.544730E+00, 1.632884E+00, 1.737451E+00, 1.839542E+00, 
        1.926015E+00, 2.027550E+00, 2.141893E+00, 2.232498E+00, 2.288891E+00, 
        2.340017E+00, 2.381759E+00, 2.383409E+00, 2.343025E+00, 2.309382E+00, 
        2.315281E+00, 2.363291E+00, 2.456813E+00, 2.608183E+00, 2.803355E+00, 
        3.023636E+00, 3.287621E+00, 3.596284E+00, 3.917294E+00, 4.268025E+00, 
        4.680325E+00, 5.108988E+00, 5.509916E+00, 5.935525E+00, 6.396799E+00, 
        6.827404E+00, 7.218567E+00, 7.614806E+00, 7.963608E+00, 8.212350E+00, 
        8.407257E+00, 8.565375E+00, 8.611701E+00, 8.551322E+00, 8.491609E+00, 
        8.456424E+00, 8.382308E+00, 8.252262E+00, 8.082789E+00, 7.892225E+00, 
        7.709093E+00, 7.552676E+00, 7.385949E+00, 7.188253E+00, 6.997005E+00, 
        6.834275E+00, 6.658732E+00, 6.452670E+00, 6.251437E+00, 6.078947E+00, 
        5.906580E+00, 5.721030E+00, 5.545298E+00, 5.382871E+00, 5.199881E+00, 
        4.998420E+00, 4.809347E+00, 4.637318E+00, 4.460711E+00, 4.285237E+00, 
        4.119970E+00, 3.957456E+00, 3.802720E+00, 3.667038E+00, 3.537688E+00, 
        3.392823E+00, 3.237844E+00, 3.083621E+00, 2.938483E+00, 2.794379E+00, 
        2.618337E+00, 2.398675E+00, 2.174743E+00, 1.973006E+00, 1.767542E+00, 
        1.559136E+00, 1.388105E+00, 1.245415E+00, 1.106893E+00, 9.899717E-01, 
        9.134314E-01, 8.420908E-01, 7.698313E-01, 7.170319E-01, 6.787253E-01, 
        6.379317E-01, 6.056904E-01, 5.856277E-01, 5.618792E-01, 5.377198E-01, 
        5.222762E-01, 5.089940E-01, 4.943599E-01, 4.838341E-01, 4.732544E-01, 
        4.588599E-01, 4.485080E-01, 4.417946E-01, 4.331451E-01 },
      rho_550 
      { 0.000000E+00, 4.324438E-01, 1.149558E+00, 1.832510E+00, 2.240700E+00,
        2.487719E+00, 2.709575E+00, 2.874102E+00, 2.991449E+00, 3.119467E+00,
        3.233615E+00, 3.348222E+00, 3.487496E+00, 3.640314E+00, 3.791807E+00,
        3.963414E+00, 4.135150E+00, 4.291491E+00, 4.451899E+00, 4.642881E+00,
        4.819324E+00, 4.985605E+00, 5.123753E+00, 5.316881E+00, 5.490964E+00,
        5.629008E+00, 5.779357E+00, 5.934624E+00, 6.069339E+00, 6.214651E+00,
        6.328005E+00, 6.433464E+00, 6.565717E+00, 6.626268E+00, 6.671826E+00,
        6.788980E+00, 6.783930E+00, 6.845387E+00, 6.873965E+00, 6.874528E+00,
        6.823557E+00, 6.769115E+00, 6.716161E+00, 6.661474E+00, 6.617604E+00,
        6.522212E+00, 6.395003E+00, 6.329291E+00, 6.214303E+00, 6.122996E+00,
        6.002230E+00, 5.893373E+00, 5.767771E+00, 5.630944E+00, 5.495732E+00,
        5.382765E+00, 5.241958E+00, 5.089607E+00, 4.982199E+00, 4.846551E+00,
        4.738482E+00, 4.579191E+00, 4.449045E+00, 4.332928E+00, 4.195822E+00,
        4.038285E+00, 3.910536E+00, 3.780957E+00, 3.657666E+00, 3.521597E+00,
        3.388630E+00, 3.246092E+00, 3.096439E+00, 2.962325E+00, 2.813099E+00,
        2.667498E+00, 2.523798E+00, 2.387101E+00, 2.249823E+00, 2.108490E+00,
        1.995867E+00, 1.880804E+00, 1.767768E+00, 1.665855E+00, 1.574605E+00,
        1.478388E+00, 1.397494E+00, 1.322468E+00, 1.263067E+00, 1.193514E+00,
        1.135380E+00, 1.081620E+00, 1.035438E+00, 9.859974E-01, 9.494125E-01,
        9.155681E-01, 8.800506E-01, 8.444502E-01, 8.097923E-01, 7.877114E-01,
        7.617647E-01, 7.387365E-01, 7.203094E-01, 7.016640E-01, 6.885593E-01,
        6.740936E-01, 6.574666E-01, 6.430402E-01, 6.322819E-01, 6.192831E-01,
        6.127931E-01, 6.031065E-01, 5.967796E-01, 5.921670E-01, 5.861204E-01,
        5.801695E-01, 5.781248E-01, 5.760360E-01, 5.787275E-01 };
      std::vector<std::vector<double>> rhoVec {rho_283,rho_550}; 
      std::vector<double> rho_dx_vec = {0.001265,0.001265}; // Spacing in eV of rho for different temperatures
                                //   283           550
      std::vector<double>   twt_vec { 6.9210E-03 , 2.1290E-02 };
      std::vector<double>     c_vec { 3.5910E+00 , 2.2343E+01 };
      std::vector<double> tbeta_vec { 5.2308E-01 , 5.0871E-01 };
      std::vector<std::vector<double>> transInfo {twt_vec,c_vec,tbeta_vec};
      
      std::vector<double> oscE_283   = { 0.205,    0.415}; 
      std::vector<double> oscW_283   = { 0.15667, 0.31333 }; 
      std::vector<double> oscE_550   = { 0.205,    0.415}; 
      std::vector<double> oscW_550   = { 0.15667, 0.31333 }; 
      std::vector<std::vector<double>> oscE_vec = {oscE_283,oscE_550};
      std::vector<std::vector<double>> oscW_vec = {oscW_283,oscW_550};
  
      auto my_MF7 = full_LEAPR(generalInfo, scatterControl, scatterInfo, temps, alphas, betas,
                            rhoVec, rho_dx_vec, transInfo, oscE_vec, oscW_vec, smin);
  
     
      auto begin = h2o.begin(), end = h2o.end();
      long lineNumber = 1;
      StructureDivision division(begin,end,lineNumber);
      njoy::ENDFtk::file::Type<7> h2o_MF7(division,begin,end,lineNumber);
  
      check_MF7(my_MF7,h2o_MF7);

    } // WHEN
  } // GIVEN


  GIVEN( "Si in SiC ENDF-B/VIII.0 input" ) {
    WHEN( "" ) {
      nphon = 100;
      mat = 43; za = 143; isabt = 0; ilog = 0; smin = 1e-75;
      awr = 27.844; spr = 2.042; npr = 1; iel = 0; ncold = 0; nsk = 0;
      nss = 0; b7 = 0; aws = 0; sps = 0; mss = 0;
      lat = 1;
      alphas = {0.5, 2.5};
      betas  = {1.0, 2.0, 6.0};
      std::vector<int> generalInfo {nphon,mat,za,isabt,ilog,lat};
      std::vector<int> scatterControl {npr,iel,ncold,nsk,nss,b7,mss};
      std::vector<double> scatterInfo {awr,spr,aws,sps};


      temps = { 300.0 };
      std::vector<double> rho_300 {
        0.000000E+00, 4.450636E-03, 1.456132E-02, 2.924358E-02, 5.227066E-02,
        8.221577E-02, 1.220293E-01, 1.759206E-01, 2.323762E-01, 2.994495E-01,
        3.817621E-01, 4.804249E-01, 6.073396E-01, 7.250425E-01, 8.596940E-01,
        1.028602E+00, 1.225708E+00, 1.483551E+00, 1.760605E+00, 2.073063E+00,
        2.444496E+00, 2.896588E+00, 3.413903E+00, 3.945487E+00, 4.500293E+00,
        5.141146E+00, 5.841038E+00, 6.654291E+00, 7.603319E+00, 8.658871E+00,
        9.969959E+00, 1.168247E+01, 1.437138E+01, 1.669329E+01, 1.772658E+01,
        1.858584E+01, 1.951784E+01, 2.057091E+01, 2.169448E+01, 2.289473E+01,
        2.419614E+01, 2.552416E+01, 2.692603E+01, 2.858531E+01, 3.087500E+01,
        2.936516E+01, 2.460824E+01, 2.145302E+01, 1.886116E+01, 1.639432E+01,
        1.539083E+01, 1.569239E+01, 1.599049E+01, 1.625384E+01, 1.655882E+01,
        1.703456E+01, 1.269941E+01, 6.601566E+00, 5.448713E+00, 6.157348E+00,
        7.068402E+00, 8.408568E+00, 1.031724E+01, 1.324826E+01, 1.894926E+01,
        2.742393E+01, 3.172214E+01, 2.669816E+01, 1.913398E+01, 1.464408E+01,
        1.174114E+01, 9.469814E+00, 7.628035E+00, 5.864667E+00, 3.898618E+00,
        2.482198E+00, 1.713215E+00, 6.918037E-01, 3.168440E-02, 0.000000E+00,
        0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00, 0.000000E+00,
        0.000000E+00, 0.000000E+00, 0.000000E+00, 5.172309E-03, 3.464289E-01,
        3.641323E+00, 1.485216E+01, 2.669359E+01, 3.211585E+01, 2.519812E+01,
        1.278351E+01, 5.693553E+00, 1.124512E+00, 0.000000E+00, 2.044555E-02,
        2.178784E-01, 8.314270E-01, 2.457587E+00, 6.139903E+00, 6.256002E+00,
        3.370903E+00, 2.669282E+00, 2.364179E+00, 2.200114E+00, 2.076445E+00,
        1.994763E+00, 1.954955E+00, 1.929053E+00, 1.937071E+00, 2.001562E+00,
        2.196842E+00, 2.097044E+00, 1.486082E+00, 9.119468E-01, 4.707854E-01 };

      std::vector<double> rho_dx_vec {0.001},
                             twt_vec { 0.0 },
                               c_vec { 0.0 },
                           tbeta_vec { 1.0 },
                            oscE_300 { },
                            oscW_300 { }; 

      std::vector<std::vector<double>> rhoVec {rho_300},
                                       transInfo {twt_vec,c_vec,tbeta_vec},
                                       oscE_vec {oscE_300},
                                       oscW_vec {oscW_300};
  
      auto my_MF7 = full_LEAPR(generalInfo, scatterControl, scatterInfo, temps, alphas, betas,
                            rhoVec, rho_dx_vec, transInfo, oscE_vec, oscW_vec, smin);
     
      auto begin = sic.begin(), end = sic.end();
      long lineNumber = 1;
      StructureDivision division(begin,end,lineNumber);
      njoy::ENDFtk::file::Type<7> sic_MF7(division,begin,end,lineNumber);

      //std::string buffer;
      //auto output = std::back_inserter(buffer);
      //my_MF7.print(output,27);
      //std::cout << buffer << std::endl;
 

      check_MF7(my_MF7,sic_MF7);
      
    } // WHEN
  } // GIVEN


  GIVEN( "Be-metal ENDF-B/VIII.0 input" ) {
    WHEN( "" ) {
      nphon = 100;
      mat = 26; za = 126; isabt = 0; ilog = 0; smin = 1e-75;
      awr = 8.93478; spr = 6.153875; npr = 1; iel = 2; ncold = 0; nsk = 0;
      nss = 0; b7 = 0; aws = 0; sps = 0; mss = 0;
      lat = 1;
      alphas = {10, 20};
      betas  = {15, 25};
      std::vector<int> generalInfo {nphon,mat,za,isabt,ilog,lat};
      std::vector<int> scatterControl {npr,iel,ncold,nsk,nss,b7,mss};
      std::vector<double> scatterInfo {awr,spr,aws,sps};


      temps = { 296.0, 400.0, 500.0 };
      std::vector<double> rho {
        0.0000E+00, 7.2477E-04, 3.7084E-03, 8.0087E-03, 1.0642E-02, 1.5897E-02, 
        2.7372E-02, 4.1843E-02, 5.0214E-02, 6.5036E-02, 8.3674E-02, 9.9329E-02, 
        1.1977E-01, 1.4296E-01, 1.6484E-01, 1.8945E-01, 2.1887E-01, 2.3537E-01, 
        2.6166E-01, 3.0003E-01, 3.4054E-01, 3.8728E-01, 4.2481E-01, 4.7598E-01, 
        5.1890E-01, 5.7400E-01, 6.2970E-01, 6.5754E-01, 7.2042E-01, 7.9118E-01, 
        8.6756E-01, 9.2948E-01, 1.0030E+00, 1.1163E+00, 1.2048E+00, 1.2870E+00, 
        1.4139E+00, 1.5249E+00, 1.6221E+00, 1.7638E+00, 1.8924E+00, 2.0388E+00, 
        2.2056E+00, 2.3709E+00, 2.5558E+00, 2.7595E+00, 3.0108E+00, 3.2603E+00, 
        3.5066E+00, 3.7442E+00, 4.0067E+00, 4.3677E+00, 4.7164E+00, 5.0820E+00, 
        5.5881E+00, 6.0898E+00, 6.5510E+00, 7.0877E+00, 7.5931E+00, 8.0736E+00, 
        8.6232E+00, 9.2283E+00, 9.9334E+00, 1.0613E+01, 1.1278E+01, 1.1973E+01, 
        1.2784E+01, 1.3744E+01, 1.4739E+01, 1.5918E+01, 1.7654E+01, 1.9834E+01, 
        2.1455E+01, 2.2574E+01, 2.3744E+01, 2.4900E+01, 2.6227E+01, 2.7931E+01, 
        2.9747E+01, 2.9884E+01, 2.7358E+01, 2.4817E+01, 2.3690E+01, 2.3242E+01, 
        2.3624E+01, 2.3473E+01, 2.2368E+01, 2.1447E+01, 2.0724E+01, 2.1121E+01, 
        2.4240E+01, 2.7607E+01, 2.7643E+01, 2.5431E+01, 2.3755E+01, 2.3377E+01, 
        2.3410E+01, 2.3504E+01, 2.3647E+01, 2.3681E+01, 2.3805E+01, 2.3714E+01, 
        2.3385E+01, 2.3050E+01, 2.2244E+01, 2.1008E+01, 1.9536E+01, 1.8341E+01, 
        1.8075E+01, 1.8606E+01, 1.9599E+01, 2.1037E+01, 2.3193E+01, 2.4016E+01, 
        2.3573E+01, 2.5664E+01, 3.0187E+01, 3.1256E+01, 2.7257E+01, 2.2765E+01, 
        1.4893E+01, 6.8192E+00, 3.8444E+00, 2.4718E+00, 1.3358E+00, 3.5968E-01, 0
      };

      std::vector<double> rho_dx_vec {0.00069552,0.00069552,0.00069552},
                             twt_vec { 0.0, 0.0, 0.0 },
                               c_vec { 0.0, 0.0, 0.0 },
                           tbeta_vec { 1.0, 1.0, 1.0 },
                                oscE { },
                                oscW { }; 

      std::vector<std::vector<double>> rhoVec {rho, rho, rho},
                                       transInfo {twt_vec,c_vec,tbeta_vec},
                                       oscE_vec {oscE, oscE, oscE},
                                       oscW_vec {oscW, oscW, oscW};
  

      auto my_MF7 = full_LEAPR(generalInfo, scatterControl, scatterInfo, temps, alphas, betas,
                            rhoVec, rho_dx_vec, transInfo, oscE_vec, oscW_vec, smin);
     
      auto begin = be.begin(), end = be.end();
      long lineNumber = 1;
      StructureDivision division(begin,end,lineNumber);
      njoy::ENDFtk::file::Type<7> be_MF7(division,begin,end,lineNumber);

      check_MF7(my_MF7,be_MF7);
      
    } // WHEN
  } // GIVEN

  */


  GIVEN( "Ortho-H ENDF-B/VIII.0 input" ) {
    WHEN( "" ) {
      nphon = 100;
      mat = 3; za = 1001.0; isabt = 0; ilog = 0; smin = 1e-75;
      awr = 0.99917; spr = 20.478; npr = 2; iel = 0; ncold = 1; nsk = 0;
      nss = 0; b7 = 0; aws = 0; sps = 0; mss = 0;
      lat = 0;
      alphas = {3.0, 10.0};
      betas  = {4.0,  8.0};
      std::vector<int> generalInfo {nphon,mat,za,isabt,ilog,lat};
      std::vector<int> scatterControl {npr,iel,ncold,nsk,nss,b7,mss};
      std::vector<double> scatterInfo {awr,spr,aws,sps};


      temps = { 20.0 };
      std::vector<double> rho { 
        0.0, .01563, .0625, .141, .25, .391, .5625, .766, 1., 1.266, 1.5625, 
        1.89, 2.25, 2.64, 3.0625, 3.52, 4., 4.6, 5.5, 7.0, 8.5, 9.2, 9.5, 9.4, 
        9.2, 8.9, 8.5, 8.0, 7.5, 7.05, 6.7, 6.4, 6.2, 6.1, 6.2, 6.45, 6.7, 6.95, 
        7.1, 6.55, 5.5, 3., 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
      };

      std::vector<double> kappa {
        9.29815e-02, 9.29815e-02, 9.31478e-02, 9.43837e-02, 9.54580e-02, 
        9.68871e-02, 9.86644e-02, 1.00813e-01, 1.03358e-01, 1.06334e-01, 
        1.09777e-01, 1.13730e-01, 1.18249e-01, 1.23395e-01, 1.29241e-01, 
        1.35875e-01, 1.43399e-01, 1.51932e-01, 1.61618e-01, 1.72622e-01, 
        1.85143e-01, 1.99416e-01, 2.15717e-01, 2.34376e-01, 2.55779e-01, 
        2.80385e-01, 3.08731e-01, 3.41443e-01, 3.79240e-01, 4.22934e-01, 
        4.73417e-01, 5.31620e-01, 5.98441e-01, 6.74614e-01, 7.98670e-01, 
        1.11800e+00, 1.40870e+00, 1.66330e+00, 1.87990e+00, 2.05540e+00, 
        2.18890e+00, 2.28030e+00, 2.33060e+00, 2.34180e+00, 2.31640e+00, 
        2.25810e+00, 2.17060e+00, 2.05850e+00, 1.92640e+00, 1.77900e+00, 
        1.62130e+00, 1.45810e+00, 1.29380e+00, 1.13270e+00, 9.78790e-01, 
        8.35450e-01, 7.05650e-01, 5.91800e-01, 4.95760e-01, 4.18820e-01, 
        3.61720e-01, 3.24640e-01, 3.07270e-01, 3.08840e-01, 3.28120e-01, 
        3.63550e-01, 4.13250e-01, 4.75100e-01, 5.46810e-01, 6.26010e-01, 
        7.10230e-01, 7.97080e-01, 8.84210e-01, 9.69420e-01, 1.05070e+00, 
        1.12610e+00, 1.19430e+00, 1.25380e+00, 1.30370e+00, 1.34320e+00, 
        1.37190e+00, 1.38970e+00, 1.39680e+00, 1.39350e+00, 1.38050e+00, 
        1.35860e+00, 1.32890e+00, 1.29240e+00, 1.25050e+00, 1.20450e+00, 
        1.15570e+00, 1.10550e+00, 1.05510e+00, 1.00590e+00, 9.58990e-01, 
        9.15410e-01, 8.76080e-01, 8.41740e-01, 8.12980e-01, 7.90250e-01, 
        7.73800e-01, 7.63720e-01, 7.59950e-01, 7.62260e-01, 7.70320e-01, 
        7.83640e-01, 8.01640e-01, 8.23670e-01, 8.48980e-01, 8.76810e-01, 
        9.06370e-01, 9.36840e-01, 9.67450e-01, 9.97440e-01, 1.02610e+00, 
        1.05290e+00, 1.07710e+00, 1.09840e+00, 1.11630e+00, 1.13060e+00, 
        1.14100e+00, 1.14760e+00, 1.15040e+00, 1.14930e+00, 1.14470e+00, 
        1.13690e+00, 1.12600e+00, 1.11260e+00, 1.09710e+00, 1.07990e+00, 
        1.06160e+00, 1.04260e+00, 1.02330e+00, 1.00440e+00, 9.86140e-01, 
        9.69010e-01, 9.53360e-01, 9.39490e-01, 9.27660e-01, 9.18060e-01, 
        9.10820e-01, 9.06020e-01, 9.03650e-01, 9.03670e-01, 9.05980e-01, 
        9.10430e-01, 9.16820e-01, 9.24910e-01, 9.34440e-01, 9.45130e-01, 
        9.56670e-01, 9.68760e-01, 9.81090e-01, 9.93350e-01, 1.00530e+00, 
        1.01660e+00, 1.02700e+00, 1.03640e+00, 1.04450e+00, 1.05130e+00, 
        1.05660e+00, 1.06030e+00, 1.06250e+00, 1.06310e+00, 1.06220e+00, 
        1.05990e+00, 1.05630e+00, 1.05150e+00, 1.04570e+00, 1.03910e+00, 
        1.03180e+00, 1.02410e+00, 1.01610e+00, 1.00800e+00, 1.00010e+00, 
        9.92500e-01, 9.85380e-01, 9.78870e-01, 9.73110e-01, 9.68200e-01, 
        9.64220e-01, 9.61220e-01, 9.59240e-01, 9.58280e-01, 9.58320e-01, 
        9.59320e-01, 9.61220e-01, 9.63940e-01, 9.67390e-01, 9.71450e-01, 
        9.76020e-01, 9.80960e-01, 9.86140e-01, 9.91450e-01, 9.96740e-01, 
        1.00190e+00, 1.00680e+00, 1.01140e+00, 1.01550e+00, 1.01920e+00, 
        1.02220e+00, 1.02460e+00, 1.02640e+00, 1.02750e+00, 1.02790e+00, 
        1.02770e+00, 1.02680e+00, 1.02540e+00, 1.02340e+00, 1.02100e+00, 
        1.01820e+00, 1.01510e+00, 1.01170e+00, 1.00820e+00, 1.00470e+00, 
        1.00120e+00, 9.97740e-01, 9.94500e-01, 9.91510e-01, 9.88820e-01, 
        9.86480e-01, 9.84540e-01, 9.83020e-01, 9.81930e-01, 9.81300e-01, 
        9.81110e-01, 9.81350e-01, 9.82000e-01, 9.83030e-01, 9.84410e-01, 
        9.86080e-01, 9.88010e-01, 9.90130e-01, 9.92400e-01, 9.94750e-01, 
        9.97130e-01, 9.99490e-01, 1.00180e+00, 1.00390e+00, 1.00590e+00, 
        1.00770e+00, 1.00920e+00, 1.01050e+00, 1.01150e+00, 1.01220e+00, 
        1.01260e+00, 1.01270e+00, 1.01250e+00, 1.01200e+00, 1.01130e+00, 
        1.01040e+00, 1.00920e+00, 1.00790e+00, 1.00650e+00, 1.00490e+00, 
        1.00340e+00, 1.00170e+00, 1.00010e+00, 9.98600e-01, 9.97150e-01, 
        9.95800e-01, 9.94600e-01, 9.93550e-01, 9.92690e-01, 9.92010e-01, 
        9.91540e-01, 9.91260e-01, 9.91190e-01, 9.91320e-01, 9.91630e-01, 
        9.92110e-01, 9.92740e-01, 9.93510e-01, 9.94400e-01, 9.95370e-01, 
        9.96410e-01, 9.97490e-01, 9.98590e-01, 9.99670e-01, 1.00070e+00, 
        1.00170e+00, 1.00260e+00, 1.00350e+00, 1.00420e+00, 1.00480e+00, 
        1.00530e+00, 1.00560e+00, 1.00580e+00, 1.00580e+00, 1.00580e+00, 
        1.00560e+00, 1.00530e+00, 1.00480e+00, 1.00430e+00, 1.00370e+00, 
        1.00310e+00, 1.00240e+00, 1.00160e+00, 1.00090e+00, 1.00010e+00, 9.99400e-01 };

      double dka = 0.05;

      std::tuple<std::vector<double>,double> kappaInfo {kappa,dka};
      std::vector<double> rho_dx_vec {0.00025},
                             twt_vec { 0.025 },
                               c_vec { 40.0  },
                           tbeta_vec { 0.475 },
                                oscE { },
                                oscW { }; 

      std::vector<std::vector<double>> rhoVec {rho, rho, rho},
                                       transInfo {twt_vec,c_vec,tbeta_vec},
                                       oscE_vec {oscE, oscE, oscE},
                                       oscW_vec {oscW, oscW, oscW};
  
      full_LEAPR(generalInfo, scatterControl, scatterInfo, temps, alphas, betas,
                            rhoVec, rho_dx_vec, transInfo, oscE_vec, oscW_vec, smin, kappaInfo);
      //auto my_MF7 = full_LEAPR(generalInfo, scatterControl, scatterInfo, temps, alphas, betas,
      //                      rhoVec, rho_dx_vec, transInfo, oscE_vec, oscW_vec, smin);
     
      auto begin = ortho_h.begin(), end = ortho_h.end();
      long lineNumber = 1;
      StructureDivision division(begin,end,lineNumber);
      njoy::ENDFtk::file::Type<7> ortho_h_MF7(division,begin,end,lineNumber);

      //check_MF7(my_MF7,ortho_h_MF7);
      
    } // WHEN
  } // GIVEN




} // TEST CASE





/*

TEST_CASE( "leapr" ){
  int nphon, ncold, lat, iel, npr, nss, b7;
  double sps, awr, aws, delta, twt, c, tbeta, dka, cfrac;
  std::vector<double> alpha, beta, temp, rho, oscE, oscW, kappa;

  GIVEN( "coarse alpha, beta grids (for quick testing)" ) {
    WHEN( "continuous, translational, and discrete oscillator options used" ) {
      nphon = 100;
      awr   = 0.99917;  iel = 0; npr = 2.0; ncold = 0; 
      aws   = 15.85316; sps = 3.8883; 
      lat   = 0;
      nss   = 0; b7 = 0;
      alpha = {.01008,   .033,     0.050406, 0.504060, 0.554466, 1.873790, 
               2.055660, 5.671180, 7.472890, 20.52030, 32.46250, 50.0};
      beta  = {0.000000, 0.006375, 0.038250, 0.06575, 0.564547, 1.77429,
               3.422470, 8.258620, 18.72180, 25.0};
      temp  = { 296.0 };                                          
      delta = 0.00255;                                           
      rho = {0, 0.0005, 0.001, 0.002, 0.0035, 0.005, 0.0075, 0.01, 0.013, 0.0165, 
             0.02, 0.0245, 0.029, 0.034, 0.0395, 0.045, 0.0506, 0.0562, 0.0622, 
             0.0686, 0.075, 0.083, 0.091, 0.099, 0.107, 0.115, 0.1197, 0.1214, 
             0.1218, 0.1195, 0.1125, 0.1065, 0.1005, 0.09542, 0.09126, 0.0871, 
             0.0839, 0.0807, 0.07798, 0.07574, 0.0735, 0.07162, 0.06974, 0.06804, 
             0.06652, 0.065, 0.0634, 0.0618, 0.06022, 0.05866, 0.0571, 0.05586, 
             0.05462, 0.0535, 0.0525, 0.0515, 0.05042, 0.04934, 0.04822, 0.04706, 
             0.0459, 0.04478, 0.04366, 0.04288, 0.04244, 0.042, 0.0, };
      twt    = 0.055556;   c = 0.0;   tbeta = 0.444444;
      oscE   = { 0.205,    0.48};                                 
      oscW   = { 0.166667, 0.333333 };                            
      dka    = 0.0;
      cfrac  = 0.0;
      auto secondaryScatterInput = std::make_tuple(nss,b7,std::vector<double>(0));
      auto oscEnergiesWeights = ranges::view::zip(oscE,oscW);

      AND_WHEN ( "alpha and beta scaling not needed (lat = 0)" ){
        lat   = 0;            
        auto out = leapr( nphon, awr, iel, npr, ncold, aws, lat, alpha, beta, temp, delta, 
                          rho, twt, c, tbeta, oscEnergiesWeights, dka, kappa, cfrac, secondaryScatterInput );
  
        std::vector<double> sabCorrect { 1.188670E+1, 1.168376E+1, 6.300891E+0, 
        1.781912E+0, 3.162210E-4, 4.692060E-4, 3.081720E-4, 2.514578E-5, 
        1.701154E-4, 6.09094E-10, 6.527907E+0, 6.500191E+0, 5.437262E+0, 
        3.734029E+0, 1.050239E-3, 1.518856E-3, 9.981392E-4, 1.995450E-4, 
        1.200483E-3, 1.282587E-8, 5.256809E+0, 5.241938E+0, 4.691417E+0, 
        3.687751E+0, 1.615068E-3, 2.303674E-3, 1.516466E-3, 3.464544E-4, 
        1.860820E-3, 3.807045E-8, 1.478889E+0, 1.480819E+0, 1.486305E+0, 
        1.470555E+0, 1.288961E-1, 2.017464E-2, 1.408493E-2, 7.462729E-3, 
        1.027139E-2, 1.373864E-5, 1.392511E+0, 1.394460E+0, 1.401735E+0, 
        1.388500E+0, 1.549325E-1, 2.190231E-2, 1.538571E-2, 8.382694E-3, 
        1.082756E-2, 1.748288E-5, 5.552458E-1, 5.564724E-1, 5.626823E-1, 
        5.679902E-1, 3.628562E-1, 5.450169E-2, 4.333900E-2, 2.265130E-2, 
        1.672152E-2, 3.550178E-4, 5.089281E-1, 5.100789E-1, 5.158600E-1, 
        5.210213E-1, 3.557392E-1, 5.783144E-2, 4.633066E-2, 2.381095E-2, 
        1.691639E-2, 4.420672E-4, 1.455618E-1, 1.459800E-1, 1.481698E-1, 
        1.499734E-1, 1.599368E-1, 9.192485E-2, 7.299081E-2, 3.854093E-2, 
        1.533922E-2, 3.925552E-3, 9.083275E-2, 9.109755E-2, 9.245460E-2, 
        9.367509E-2, 1.060097E-1, 8.399707E-2, 7.119046E-2, 4.314996E-2, 
        1.459694E-2, 6.322492E-3, 6.735360E-3, 6.756551E-3, 6.863877E-3, 
        6.958397E-3, 8.741588E-3, 1.328617E-2, 1.925801E-2, 3.149450E-2, 
        2.194943E-2, 1.818368E-2, 9.378210E-4, 9.408032E-4, 9.559292E-4, 
        9.692905E-4, 1.230320E-3, 2.054059E-3, 3.601626E-3, 1.031145E-2, 
        1.933435E-2, 1.935657E-2, 6.228358E-5, 6.248242E-5, 6.348621E-5, 
        6.436519E-5, 8.208659E-5, 1.425013E-4, 2.777029E-4, 1.245381E-3, 
        6.606240E-3, 1.044380E-2};
  
        std::vector<double> effectiveTemps {1397.671178314},
                            debyeWaller    {0.2352041964  };

        THEN( "the scattering law is correctly returned on a scaled grid" ){
          REQUIRE( ranges::equal(std::get<0>(out),sabCorrect,    equal) );
          REQUIRE( ranges::equal(std::get<1>(out),effectiveTemps,equal) );
          REQUIRE( ranges::equal(std::get<2>(out),debyeWaller,   equal) );
        } // THEN
      } // AND WHEN
      AND_WHEN ( "alpha and beta scaling is requested (lat = 1)" ){
        lat   = 1;            
        auto out = leapr( nphon, awr, iel, npr, ncold, aws, lat, alpha, beta, temp, delta, 
                          rho, twt, c, tbeta, oscEnergiesWeights, dka, kappa, cfrac, secondaryScatterInput );
        std::vector<double> sabCorrect { 1.193557E+1, 1.173342E+1, 6.361294E+0, 
          1.816275E+0, 3.124950E-4, 4.639237E-4, 3.122224E-4, 7.882191E-5, 
          1.011487E-5,7.888274E-10, 6.555078E+0, 6.527275E+0, 5.467981E+0, 
          3.766171E+0, 1.038090E-3, 1.501995E-3, 1.011017E-3, 5.830163E-4, 
          8.176285E-5, 1.597145E-8, 5.278888E+0, 5.263947E+0, 4.715237E+0, 
          3.714425E+0, 1.596591E-3, 2.278303E-3, 1.535826E-3, 9.525713E-4, 
          1.452123E-4, 4.673899E-8, 1.486402E+0, 1.488315E+0, 1.493827E+0, 
          1.477974E+0, 1.316516E-1, 1.997813E-2, 1.423771E-2, 1.018858E-2, 
          4.169200E-3, 1.556686E-5, 1.399711E+0, 1.401644E+0, 1.408946E+0, 
          1.395614E+0, 1.581786E-1, 2.169138E-2, 1.555017E-2, 1.108814E-2, 
          4.735061E-3, 1.974070E-5, 5.593066E-1, 5.605274E-1, 5.667146E-1, 
          5.720510E-1, 3.663735E-1, 5.415787E-2, 4.367399E-2, 2.386446E-2, 
          1.252101E-2, 3.835034E-4, 5.127808E-1, 5.139310E-1, 5.196984E-1, 
          5.248412E-1, 3.591896E-1, 5.751228E-2, 4.667544E-2, 2.488494E-2, 
          1.293548E-2, 4.758096E-4, 1.473396E-1, 1.477655E-1, 1.499618E-1, 
          1.517737E-1, 1.616866E-1, 9.247084E-2, 7.336139E-2, 3.867818E-2, 
          1.372647E-2, 4.048636E-3, 9.213506E-2, 9.240120E-2, 9.376556E-2, 
          9.499839E-2, 1.073307E-1, 8.467699E-2, 7.158413E-2, 4.326283E-2, 
          1.351358E-2, 6.455035E-3, 6.930727E-3, 6.952350E-3, 7.061863E-3, 
          7.158157E-3, 8.973773E-3, 1.358021E-2, 1.959090E-2, 3.179135E-2, 
          2.203982E-2, 1.822503E-2, 9.778850E-4, 9.809690E-4, 9.966260E-4, 
          1.010448E-3, 1.279979E-3, 2.127476E-3, 3.710551E-3, 1.050751E-2, 
          1.954142E-2, 1.944808E-2, 6.624281E-5, 6.645256E-5, 6.751134E-5, 
          6.843838E-5, 8.710059E-5, 1.505165E-4, 2.916502E-4, 1.290367E-3, 
          6.740754E-3, 1.056903E-2};
        std::vector<double> effectiveTemps {1397.671178314},
                            debyeWaller    {0.2352041964  };

        THEN( "the scattering law is correctly returned on a scaled grid" ){
          REQUIRE( ranges::equal(std::get<0>(out),sabCorrect,    equal) );
          REQUIRE( ranges::equal(std::get<1>(out),effectiveTemps,equal) );
          REQUIRE( ranges::equal(std::get<2>(out),debyeWaller,   equal) );
        } // THEN
      } // AND WHEN
      AND_WHEN( "material is cold" ){
        temp  = { 100.0 };                                          
        auto out = leapr( nphon, awr, iel, npr, ncold, aws, lat, alpha, beta, temp, delta, 
                          rho, twt, c, tbeta, oscEnergiesWeights, dka, kappa, cfrac, secondaryScatterInput );

        std::vector<double> sabCorrect { 1.191059E+1, 1.170723E+1, 6.313288E+0, 
          1.785018E+0, 1.591388E-5, 1.977783E-5, 3.162744E-5, 5.709636E-5, 
          8.766729E-6, 2.004754E-9, 6.570177E+0, 6.542268E+0, 5.472154E+0, 
          3.757436E+0, 5.258778E-5, 6.470506E-5, 1.031839E-4, 1.860698E-4, 
          2.855014E-5, 3.566762E-8, 5.308403E+0, 5.293371E+0, 4.737147E+0, 
          3.723126E+0, 8.070512E-5, 9.875910E-5, 1.572592E-4, 2.834417E-4, 
          4.352879E-5, 8.395187E-8, 1.616967E+0, 1.619044E+0, 1.624822E+0, 
          1.607225E+0, 1.252489E-1, 9.603589E-4, 1.504800E-3, 2.708091E-3, 
          4.368001E-4, 3.862403E-5, 1.535347E+0, 1.537456E+0, 1.545270E+0, 
          1.530291E+0, 1.539718E-1, 1.052854E-3, 1.647889E-3, 2.966128E-3, 
          4.811426E-4, 4.881106E-5, 7.499578E-1, 7.516101E-1, 7.599277E-1, 
          7.668914E-1, 4.630478E-1, 4.192283E-3, 4.983860E-3, 9.056953E-3, 
          1.687982E-3, 5.798369E-4, 7.055258E-1, 7.071121E-1, 7.150962E-1, 
          7.220613E-1, 4.659443E-1, 5.252341E-3, 5.387029E-3, 9.805895E-3, 
          1.860422E-3, 7.168131E-4, 3.175353E-1, 3.183719E-1, 3.225950E-1, 
          3.262726E-1, 3.279453E-1, 6.996016E-2, 1.131200E-2, 2.105390E-2, 
          5.444883E-3, 6.518448E-3, 2.396012E-1, 2.402501E-1, 2.435240E-1, 
          2.463885E-1, 2.625363E-1, 9.412094E-2, 1.390007E-2, 2.456444E-2, 
          7.245434E-3, 9.318142E-3, 5.184109E-2, 5.199249E-2, 5.275752E-2, 
          5.342847E-2, 6.426167E-2, 6.634176E-2, 3.319597E-2, 2.811754E-2, 
          1.693344E-2, 1.658023E-2, 1.635642E-2, 1.640538E-2, 1.665411E-2, 
          1.687429E-2, 2.082778E-2, 2.691698E-2, 2.342642E-2, 1.982073E-2, 
          1.878867E-2, 1.665706E-2, 3.443731E-3, 3.454243E-3, 3.507459E-3, 
          3.554265E-3, 4.448901E-3, 6.591632E-3, 8.198762E-3, 9.342037E-3, 
          1.423867E-2, 1.370927E-2};

        std::vector<double> effectiveTemps {1366.18152086 },
                            debyeWaller    {5.68713599E-2 };

        THEN( "the scattering law is correctly returned on a scaled grid" ){
          REQUIRE( ranges::equal(std::get<0>(out),sabCorrect,    equal) );
          REQUIRE( ranges::equal(std::get<1>(out),effectiveTemps,equal) );
          REQUIRE( ranges::equal(std::get<2>(out),debyeWaller,   equal) );
        } // THEN

      } // AND WHEN

      AND_WHEN( "material is hot" ){
        temp  = { 600.0 };                                          
        auto out = leapr( nphon, awr, iel, npr, ncold, aws, lat, alpha, beta, temp, delta, 
                          rho, twt, c, tbeta, oscEnergiesWeights, dka, kappa, cfrac, secondaryScatterInput );

        std::vector<double> sabCorrect { 1.182050E+1, 1.161877E+1, 6.267602E+0, 
          1.775094E+0, 2.538069E-3, 1.274758E-3, 3.453815E-6, 2.589951E-7, 
          2.71834E-15, 0.000000E+0, 6.414739E+0, 6.387587E+0, 5.344998E+0, 
          3.674163E+0, 8.254109E-3, 4.092687E-3, 3.690549E-5, 2.840638E-6, 
          5.29688E-13, 3.88683E-20, 5.121307E+0, 5.106933E+0, 4.572468E+0, 
          3.597790E+0, 1.253030E-2, 6.194625E-3, 8.554847E-5, 6.672075E-6, 
          2.151270E-6, 1.84458E-13, 1.207943E+0, 1.209742E+0, 1.214528E+0, 
          1.204356E+0, 1.856763E-1, 5.413021E-2, 8.248944E-3, 8.973797E-4, 
          1.418359E-4, 8.605015E-8, 1.119380E+0, 1.120947E+0, 1.126974E+0, 
          1.119192E+0, 2.084083E-1, 5.868859E-2, 9.891214E-3, 1.101023E-3, 
          1.650710E-4, 1.300575E-7, 3.296997E-1, 3.309046E-1, 3.366573E-1, 
          3.392579E-1, 2.810407E-1, 1.296941E-1, 5.873684E-2, 8.093005E-3, 
          8.179071E-4, 2.181782E-5, 2.937824E-1, 2.948494E-1, 3.003364E-1, 
          3.026732E-1, 2.655510E-1, 1.336042E-1, 6.440555E-2, 9.156526E-3, 
          9.122884E-4, 3.149881E-5, 6.493049E-2, 6.513407E-2, 6.616440E-2, 
          6.709328E-2, 8.129466E-2, 9.848076E-2, 9.553444E-2, 3.456580E-2, 
          4.882244E-3, 1.157556E-3, 3.851669E-2, 3.863880E-2, 3.925635E-2, 
          3.976768E-2, 4.933831E-2, 6.852243E-2, 7.992601E-2, 4.483372E-2, 
          8.770187E-3, 2.590310E-3, 1.913698E-3, 1.919809E-3, 1.950662E-3, 
          1.977712E-3, 2.520742E-3, 4.266602E-3, 7.722683E-3, 2.172338E-2, 
          2.834815E-2, 2.092619E-2, 1.642067E-4, 1.647312E-4, 1.673778E-4, 
          1.696941E-4, 2.169625E-4, 3.781273E-4, 7.515546E-4, 3.461227E-3, 
          1.537487E-2, 2.126633E-2, 5.110181E-6, 5.126508E-6, 5.208941E-6, 
          5.281167E-6, 6.762964E-6, 1.198200E-5, 2.513165E-5, 1.576405E-4, 
          2.127419E-3, 5.732186E-3 };

        std::vector<double> effectiveTemps {1507.296344984},
                            debyeWaller    {0.782347162724};

        THEN( "the scattering law is correctly returned on a scaled grid" ){
          REQUIRE( ranges::equal(std::get<0>(out),sabCorrect,    equal) );
          REQUIRE( ranges::equal(std::get<1>(out),effectiveTemps,equal) );
          REQUIRE( ranges::equal(std::get<2>(out),debyeWaller,   equal) );
        } // THEN

      } // AND WHEN
    } // WHEN
  } // GIVEN

  GIVEN( "test9 (simple H in H2O)" ) {
    WHEN( "continuous, translational, and discrete oscillator options used" ) {
      nphon = 100;
      awr   = 0.99917;  iel = 0; ncold = 0; 
      aws   = 15.85316; sps = 3.8883; 
      lat   = 1;
      nss   = 0; b7 = 0;
      alpha = { 0.01008, 0.015, 0.0252, 0.033, 0.050406, 0.0756, 0.100812, 
      0.151218, 0.201624, 0.252030, 0.302436, 0.352842, 0.403248, 0.453654, 
      0.504060, 0.554466, 0.609711, 0.670259, 0.736623, 0.809349, 0.889061, 
      0.976435, 1.072130, 1.177080, 1.292110, 1.418220, 1.556330, 1.707750, 
      1.873790, 2.055660, 2.255060, 2.473520, 2.712950, 2.975460, 3.263080, 
      3.578320, 3.923900, 4.302660, 4.717700, 5.172560, 5.671180, 6.217580, 
      6.816500, 7.472890, 8.192280, 8.980730, 9.844890, 10.79190, 11.83030, 
      12.96740, 14.21450, 15.58150, 17.07960, 18.72080, 20.52030, 22.49220, 
      24.65260, 27.02160, 29.61750, 32.46250, 35.58160, 38.99910, 42.74530, 
      46.85030, 50.0 };
      beta = { 0.000000, 0.006375, 0.012750, 0.025500, 0.038250, 0.051000, 
      0.065750, 0.0806495, 0.120974, 0.161299, 0.241949, 0.322598, 0.403248, 
      0.483897, 0.564547, 0.645197, 0.725846, 0.806496, 0.887145, 0.967795, 
      1.048440, 1.129090, 1.209740, 1.290390, 1.371040, 1.451690, 1.532340, 
      1.612990, 1.693640, 1.774290, 1.854940, 1.935590, 2.016240, 2.096890, 
      2.177540, 2.258190, 2.338840, 2.419490, 2.500140, 2.580790, 2.669500, 
      2.767090, 2.874450, 2.992500, 3.122350, 3.265300, 3.422470, 3.595360, 
      3.785490, 3.994670, 4.224730, 4.477870, 4.756310, 5.062580, 5.399390, 
      5.769970, 6.177660, 6.626070, 7.119240, 7.661810, 8.258620, 8.915110, 
      9.637220, 10.43200, 11.30510, 12.26680, 13.32430, 14.48670, 15.76600, 
      17.17330, 18.72180, 20.42450, 22.29760, 24.35720, 25.0 };


      temp  = { 296.0 };                                          
      delta = 0.00255;                                           
      rho = {0, 0.0005, 0.001, 0.002, 0.0035, 0.005, 0.0075, 0.01, 0.013, 0.0165, 
             0.02, 0.0245, 0.029, 0.034, 0.0395, 0.045, 0.0506, 0.0562, 0.0622, 
             0.0686, 0.075, 0.083, 0.091, 0.099, 0.107, 0.115, 0.1197, 0.1214, 
             0.1218, 0.1195, 0.1125, 0.1065, 0.1005, 0.09542, 0.09126, 0.0871, 
             0.0839, 0.0807, 0.07798, 0.07574, 0.0735, 0.07162, 0.06974, 0.06804, 
             0.06652, 0.065, 0.0634, 0.0618, 0.06022, 0.05866, 0.0571, 0.05586, 
             0.05462, 0.0535, 0.0525, 0.0515, 0.05042, 0.04934, 0.04822, 0.04706, 
             0.0459, 0.04478, 0.04366, 0.04288, 0.04244, 0.042, 0.0, };
      twt    = 0.055556;   c = 0.0;   tbeta = 0.444444;
      oscE   = { 0.205,    0.48};                                 
      oscW   = { 0.166667, 0.333333 };                            
      dka    = 0.0;
      cfrac  = 0.0;
      auto oscEnergiesWeights = ranges::view::zip(oscE,oscW);
      auto secondaryScatterInput = std::make_tuple(nss,b7,std::vector<double>(0));


      auto out = leapr( nphon, awr, iel, npr, ncold, aws, lat, alpha, beta, temp, delta, 
                        rho, twt, c, tbeta, oscEnergiesWeights, dka, kappa, cfrac, secondaryScatterInput );
      std::vector<double> sab = std::get<0>(out),
      sab_0_99 { 1.193557E+1, 1.173342E+1, 1.115293E+1, 9.042512E+0, 6.361295E+0, 
      3.862215E+0, 1.816276E+0, 6.976093E-1, 1.990600E-2, 5.694252E-4, 
      3.233509E-4, 2.982178E-4, 2.996330E-4, 2.961039E-4, 3.053551E-4, 
      3.177549E-4, 3.265274E-4, 3.366876E-4, 3.482198E-4, 3.577151E-4, 
      3.690864E-4, 3.825873E-4, 3.943943E-4, 4.067551E-4, 4.197706E-4, 
      4.318199E-4, 4.422418E-4, 4.513514E-4, 4.592862E-4, 4.676470E-4, 
      4.766943E-4, 4.858306E-4, 4.954565E-4, 5.093038E-4, 5.232462E-4, 
      5.360422E-4, 5.478218E-4, 5.586641E-4, 5.670666E-4, 5.670803E-4, 
      5.574108E-4, 5.388335E-4, 5.099740E-4, 4.642119E-4, 4.116225E-4, 
      3.618202E-4, 3.183550E-4, 2.801315E-4, 2.468202E-4, 2.184691E-4, 
      1.934883E-4, 1.714009E-4, 1.502529E-4, 1.302423E-4, 1.135101E-4, 
      9.827482E-5, 8.024491E-5, 1.829132E-5, 2.490368E-7, 1.936146E-7, 
      3.065288E-7, 1.477736E-7, 1.342006E-7, 1.349326E-7, 8.946933E-8, 
      4.469996E-8, 2.540386E-8, 8.788727E-9, 2.848220E-9, 1.400958E-8, 
      4.429853E-8, 7.649074E-8, 6.102628E-8, 2.022216E-8, 1.531153E-8, 
      9.770919E+0, 9.675363E+0, 9.353680E+0, 8.135114E+0, 6.428072E+0, 
      4.614642E+0, 2.783247E+0, 1.467645E+0, 1.338837E-1, 5.248679E-3, 
      4.879553E-4, 4.447880E-4, 4.446880E-4, 4.413456E-4, 4.539861E-4, 
      4.719887E-4, 4.855082E-4, 5.005741E-4, 5.173926E-4, 5.318829E-4, 
      5.487859E-4, 5.685026E-4, 5.863043E-4, 6.046558E-4, 6.238830E-4},
      sab_500_599 { 1.933185E-3, 1.713349E-3, 1.505545E-3, 1.310263E-3, 
      1.142633E-3, 9.898227E-4, 7.678806E-4, 2.175561E-4, 2.742972E-5, 
      2.130389E-5, 2.803759E-3, 1.459825E-5, 1.326542E-5, 1.330020E-5, 
      8.946200E-6, 4.492583E-6, 2.576754E-6, 9.327469E-7, 3.202771E-7, 
      1.394743E-6, 3.535918E-4, 7.481312E-6, 6.082111E-6, 2.035847E-6, 
      1.493183E-6, 2.966986E+0, 2.966118E+0, 2.965242E+0, 2.940778E+0, 
      2.894150E+0, 2.813548E+0, 2.692185E+0, 2.543587E+0, 2.042263E+0, 
      1.491382E+0, 5.977390E-1, 1.654652E-1, 3.429444E-2, 8.198898E-3, 
      4.822959E-3, 4.631539E-3, 4.751903E-3, 4.899810E-3, 5.053373E-3, 
      5.212293E-3, 5.378678E-3, 5.552414E-3, 5.730512E-3, 5.909620E-3, 
      6.085454E-3, 6.252787E-3, 6.407194E-3, 6.548134E-3, 6.679658E-3, 
      6.808422E-3, 6.940981E-3, 7.083453E-3, 7.240706E-3, 7.413686E-3, 
      7.595304E-3, 7.774139E-3, 7.936868E-3, 8.066715E-3, 8.141638E-3, 
      8.138855E-3, 8.031090E-3, 7.785123E-3, 7.376013E-3, 6.797479E-3, 
      6.102597E-3, 5.392798E-3, 4.745504E-3, 4.178315E-3, 3.686777E-3, 
      3.264617E-3, 2.897048E-3, 2.569013E-3, 2.260247E-3, 1.970561E-3, 
      1.719532E-3, 1.490394E-3, 1.144583E-3, 3.557205E-4, 6.217376E-5, 
      6.932661E-5, 4.867260E-3, 3.264047E-5, 2.965699E-5, 2.967388E-5, 
      2.011452E-5, 1.013746E-5, 5.839988E-6, 2.179936E-6, 7.984628E-7, 
      3.130141E-6, 1.080462E-3, 1.663127E-5, 1.365138E-5, 4.595530E-6, 
      3.350952E-6 },
      sab_2000_2099 {2.82225E-2, 2.54145E-2, 2.28959E-2, 2.05665E-2, 1.83361E-2, 
      1.59629E-2, 1.28038E-2, 8.69996E-3, 6.09553E-3, 1.30761E-2, 2.35583E-2, 
      7.47693E-3, 2.71252E-3, 2.52189E-3, 1.93672E-3, 1.12279E-3, 7.04719E-4, 
      3.81741E-4, 3.30999E-4, 3.49701E-4, 1.33180E-2, 1.37823E-3, 1.25240E-3, 
      5.09148E-4, 3.82105E-4, 6.03086E-1, 6.04395E-1, 6.05716E-1, 6.08336E-1, 
      6.10986E-1, 6.13671E-1, 6.16068E-1, 6.16343E-1, 6.17192E-1, 6.10888E-1, 
      5.87412E-1, 5.45918E-1, 4.91873E-1, 4.30847E-1, 3.66619E-1, 3.03499E-1, 
      2.46086E-1, 1.96338E-1, 1.54737E-1, 1.22182E-1, 9.78521E-2, 8.03188E-2, 
      6.84620E-2, 6.09220E-2, 5.64906E-2, 5.42060E-2, 5.33368E-2, 5.33679E-2, 
      5.39279E-2, 5.47784E-2, 5.57793E-2, 5.68357E-2, 5.78748E-2, 5.88597E-2, 
      5.97560E-2, 6.05148E-2, 6.11068E-2, 6.15032E-2, 6.16633E-2, 6.15652E-2, 
      6.11320E-2, 6.02531E-2, 5.88112E-2, 5.67071E-2, 5.38964E-2, 5.04132E-2, 
      4.64293E-2, 4.21842E-2, 3.80052E-2, 3.41351E-2, 3.06999E-2, 2.76853E-2, 
      2.49871E-2, 2.24926E-2, 2.00920E-2, 1.75198E-2, 1.41606E-2, 9.88724E-3, 
      7.24035E-3, 1.42990E-2, 2.43807E-2, 8.87663E-3, 3.24062E-3, 2.97052E-3, 
      2.29914E-3, 1.35337E-3, 8.55207E-4, 4.72120E-4, 4.03303E-4, 4.27195E-4, 
      1.37061E-2, 1.65404E-3, 1.47705E-3, 6.12595E-4, 4.62908E-4};


      std::vector<double> effectiveTemps {1397.671178},
                          debyeWaller    {0.235204196 };

      THEN( "the scattering law is correctly returned on a scaled grid" ){
        REQUIRE( ranges::equal(std::get<1>(out),effectiveTemps,equal) );
        REQUIRE( ranges::equal(std::get<2>(out),debyeWaller,   equal) );
        checkPartOfVec( sab, sab_0_99,         0 );
        checkPartOfVec( sab, sab_500_599,    500 );
        checkPartOfVec( sab, sab_2000_2099, 2000 );
      } // THEN

    } // WHEN
  } // GIVEN
*/
/*

  GIVEN( "Beryllium metal" ) {
    WHEN( "continuous and coherent elastic options used" ) {
      nphon = 100;
      awr   = 8.93478;  iel = 2; npr = 1; ncold = 0; 
      aws   = 0; sps = 0; 
      lat   = 1;
      nss   = 0; b7 = 0;
      alpha = { 3.504421E-3, 4.022632E-3, 6.983712E-3, 8.016418E-3, 1.391734E-2, 
        1.597535E-2, 2.773489E-2, 3.183614E-2, 5.527088E-2, 6.344398E-2, 
        1.101454E-1, 1.264330E-1, 2.195010E-1, 2.519593E-1, 4.374279E-1, 
        5.021119E-1, 8.717190E-1, 1.000623E+0, 1.737187E+0, 1.994071E+0, 
        3.461916E+0, 3.973842E+0, 6.899007E+0, 7.919187E+0, 1.374854E+1, 
        1.578158E+1, 2.739849E+1, 3.145000E+1, 5.460050E+1, 6.267447E+1, 
        7.707822E+1, 8.847604E+1 };

      beta = { 0.000000E+0, 1.091104E-1, 4.091640E-1, 5.182744E-1, 8.183280E-1, 
               9.274384E-1, 1.227492E+0, 1.336602E+0, 1.636656E+0, 1.745766E+0, 
               2.045820E+0, 2.154930E+0, 2.454984E+0, 2.564094E+0, 2.864148E+0, 
               2.973259E+0, 3.273312E+0, 3.382423E+0, 5.663993E+0, 7.072014E+0, 
               1.302265E+1, 1.625997E+1, 2.994166E+1, 3.738491E+1, 6.884185E+1, 
               8.595535E+1, 1.582811E+2, 1.976285E+2 };

      temp  = { 296.0 };                                          
      delta = 0.00069552;
      rho = {0.0000E+0, 7.2477E-4, 3.7084E-3, 8.0087E-3, 1.0642E-2, 1.5897E-2, 
             2.7372E-2, 4.1843E-2, 5.0214E-2, 6.5036E-2, 8.3674E-2, 9.9329E-2, 
             1.1977E-1, 1.4296E-1, 1.6484E-1, 1.8945E-1, 2.1887E-1, 2.3537E-1, 
             2.6166E-1, 3.0003E-1, 3.4054E-1, 3.8728E-1, 4.2481E-1, 4.7598E-1, 
             5.1890E-1, 5.7400E-1, 6.2970E-1, 6.5754E-1, 7.2042E-1, 7.9118E-1, 
             8.6756E-1, 9.2948E-1, 1.0030E+0, 1.1163E+0, 1.2048E+0, 1.2870E+0, 
             1.4139E+0, 1.5249E+0, 1.6221E+0, 1.7638E+0, 1.8924E+0, 2.0388E+0, 
             2.2056E+0, 2.3709E+0, 2.5558E+0, 2.7595E+0, 3.0108E+0, 3.2603E+0, 
             3.5066E+0, 3.7442E+0, 4.0067E+0, 4.3677E+0, 4.7164E+0, 5.0820E+0, 
             5.5881E+0, 6.0898E+0, 6.5510E+0, 7.0877E+0, 7.5931E+0, 8.0736E+0, 
             8.6232E+0, 9.2283E+0, 9.9334E+0, 1.0613E+1, 1.1278E+1, 1.1973E+1, 
             1.2784E+1, 1.3744E+1, 1.4739E+1, 1.5918E+1, 1.7654E+1, 1.9834E+1, 
             2.1455E+1, 2.2574E+1, 2.3744E+1, 2.4900E+1, 2.6227E+1, 2.7931E+1, 
             2.9747E+1, 2.9884E+1, 2.7358E+1, 2.4817E+1, 2.3690E+1, 2.3242E+1, 
             2.3624E+1, 2.3473E+1, 2.2368E+1, 2.1447E+1, 2.0724E+1, 2.1121E+1, 
             2.4240E+1, 2.7607E+1, 2.7643E+1, 2.5431E+1, 2.3755E+1, 2.3377E+1, 
             2.3410E+1, 2.3504E+1, 2.3647E+1, 2.3681E+1, 2.3805E+1, 2.3714E+1, 
             2.3385E+1, 2.3050E+1, 2.2244E+1, 2.1008E+1, 1.9536E+1, 1.8341E+1, 
             1.8075E+1, 1.8606E+1, 1.9599E+1, 2.1037E+1, 2.3193E+1, 2.4016E+1, 
             2.3573E+1, 2.5664E+1, 3.0187E+1, 3.1256E+1, 2.7257E+1, 2.2765E+1, 
             1.4893E+1, 6.8192E+0, 3.8444E+0, 2.4718E+0, 1.3358E+0, 3.5968E-1,0};
      twt    = 0.0;   c = 0.0;   tbeta = 1.0;
      oscE   = {};                                 
      oscW   = {};                            
      dka    = 0.0;
      cfrac  = 0.0;
 
      auto oscEnergiesWeights = ranges::view::zip(oscE,oscW);
      auto secondaryScatterInput = std::make_tuple(nss,b7,std::vector<double>(0));


      auto out = leapr( nphon, awr, iel, npr, ncold, aws, lat, alpha, beta, temp, delta, 
                        rho, twt, c, tbeta, oscEnergiesWeights, dka, kappa, cfrac, secondaryScatterInput );

      std::vector<double> sab_alpha_0 { 8.64142E-5, 8.45470E-5, 1.22065E-4, 
      1.26150E-4, 1.66740E-4, 1.89049E-4, 2.77557E-4, 3.32110E-4, 5.68091E-4, 
      6.79697E-4, 1.21621E-3, 1.39873E-3, 8.78432E-4, 9.43453E-4, 7.28456E-4, 
      5.75077E-4, 6.30556E-4, 6.67230E-5, 3.01255E-7, 3.22515E-10, 5.89579E-17, 
      3.16696E-21, 0, 0, 0, 0, 0, 0 },
      sab_alpha_1 { 9.91904E-5, 9.70488E-5, 1.40103E-4, 1.44791E-4, 1.91367E-4, 
      2.16967E-4, 3.18530E-4, 3.81129E-4, 6.51920E-4, 7.79990E-4, 1.39565E-3, 
      1.60510E-3, 1.00804E-3, 1.08266E-3, 8.35956E-4, 6.59953E-4, 7.23628E-4, 
      7.66242E-5, 3.96861E-7, 4.87658E-10, 1.06910E-16, 6.86348E-21},
      sab_alpha_10 {2.70046E-3, 2.65063E-3, 3.76450E-3, 3.88592E-3, 5.07818E-3, 
      5.73432E-3, 8.33423E-3, 9.93737E-3, 1.68898E-2, 2.01808E-2, 3.600775E-2, 
      4.13947E-2, 2.60717E-2, 2.79984E-2, 2.16924E-2, 1.71854E-2, 1.887945E-2, 
      2.28032E-3, 2.85967E-4, 9.47584E-6, 5.62397E-10, 1.80854E-12, 3.57178E-24, 
      0, 0, 0, 0, 0},
      sab_alpha_20 {3.75140E-2, 3.83295E-2, 4.51205E-2, 4.62253E-2, 5.21572E-2, 
      5.53431E-2, 6.70168E-2, 7.38845E-2, 1.04935E-1, 1.19801E-1, 1.90161E-1, 
      2.14452E-1, 1.54299E-1, 1.64337E-1, 1.43548E-1, 1.26969E-1, 1.42634E-1, 
      7.61906E-2, 7.98099E-2, 5.53260E-2, 4.12479E-3, 6.43247E-4, 3.05770E-8, 
      4.80870E-11, 4.46480E-25, 0, 0, 0},
      sab_alpha_30 {1.49048E-9, 1.57328E-9, 1.82453E-9, 1.92517E-9, 2.23031E-9, 
      2.35245E-9, 2.72251E-9, 2.87052E-9, 3.31865E-9, 3.49776E-9, 4.03965E-9, 
      4.25607E-9, 4.91037E-9, 5.17150E-9, 5.96040E-9, 6.27502E-9, 7.22481E-9, 
      7.60333E-9, 2.16204E-8, 4.03234E-8, 4.68935E-7, 1.57899E-6, 1.09317E-4, 
      6.22095E-4, 2.37926E-2, 2.16143E-2, 1.67615E-7, 3.61578E-16},
      bragg_0_49 { 1.59285E-3, 0.00000E+0, 5.21980E-3, 8.97802E-3, 6.37138E-3, 
      1.08350E-2, 6.81265E-3, 4.71521E-2, 1.15912E-2, 1.20496E-2, 1.43356E-2, 
      0.00000E+0, 1.56594E-2, 2.07339E-2, 1.72522E-2, 0.00000E+0, 1.95554E-2, 
      2.78308E-2, 2.08792E-2, 4.48901E-3, 2.20308E-2, 3.49609E-2, 2.24720E-2, 
      2.59620E-2, 2.54855E-2, 5.41751E-3, 2.72506E-2, 7.85868E-3, 2.99950E-2, 
      0.00000E+0, 3.07053E-2, 7.40339E-3, 3.52148E-2, 2.07394E-2, 3.65386E-2, 
      6.78675E-3, 3.81315E-2, 3.98609E-2, 3.98211E-2, 0.00000E+0, 4.11449E-2, 
      2.55823E-2, 4.29100E-2, 1.25253E-2, 4.50409E-2, 1.83381E-2, 4.63647E-2, 
      6.02482E-3, 4.69782E-2, 1.19707E-2},
      bragg_500_549 { 4.706568, 1.660490E-1, 4.728253, 6.495063E-1, 4.818357, 
      1.398833, 5, 1.398833, 5.168748E19, 1.043008E-6, 5.428632E19, 2.198066E-8, 
      5.434033E19, 1.033517E-6, 5.685315E19, 1.197387E-6, 5.970882E19, 
      4.658899E-9, 5.971524E19, 3.881319E-8, 5.978069E19, 1.086285E-8, 
      5.982474E19, 9.308766E-9, 5.982846E19, 1.551413E-9, 5.983052E19, 
      6.205273E-9, 5.984922E19, 3.023788E-8, 5.993128E19, 9.226990E-7, 
      6.226481E19, 1.332173E-6, 6.538844E19, 2.373716E-8, 6.546593E19, 
      9.527974E-7, 6.792245E19, 1.411565E-6, 7.132786E19, 7.103531E-9, 
      7.136292E19, 5.111142E-8, 7.151167E19, 1.915366E-8, 7.154989E19, 
      1.056123E-7, 7.175038E19, 8.802859E-7 };
        
      std::vector<double> effectiveTemps {433.38329725},
                          debyeWaller    {0.603338907 };

      THEN( "the scattering law is correctly returned on a scaled grid" ){
        REQUIRE( ranges::equal(std::get<1>(out),effectiveTemps,equal) );
        REQUIRE( ranges::equal(std::get<2>(out),debyeWaller,   equal) );

        auto braggOut = std::get<3>(out);
        if (auto* bragg = std::get_if<std::vector<double>>(&braggOut)) {
          checkPartOfVec( *bragg, bragg_0_49,         0 );
          checkPartOfVec( *bragg, bragg_500_549,    500 );
        }
        else{ REQUIRE(false); }
        
        checkPartOfVec( std::get<0>(out), sab_alpha_0,     0  );
        checkPartOfVec( std::get<0>(out), sab_alpha_1,  28*1  );
        checkPartOfVec( std::get<0>(out), sab_alpha_10, 28*10 );
        checkPartOfVec( std::get<0>(out), sab_alpha_20, 28*20 );
        checkPartOfVec( std::get<0>(out), sab_alpha_30, 28*30 );

      } // THEN

    } // WHEN
  } // GIVEN

  GIVEN( "Heavy water" ) {
    WHEN( "continuous, translational, discrete osc, and skold options used"){
      nphon = 200;
      awr = 1.9968; npr = 2; iel = 0; ncold = 0; 
      aws = 0.0;
      lat = 1;
      alpha = { 1.637E-4, 5.012E-4, 5.892E-3, 7.4571E-3, 4.1899E-2, 4.59191E-2, 
      2.12114E-1, 2.21044E-1, 1.43, 1.79, 52.5, 56.72, 94, 113, 892.9, 922.9, 
      1192.9, 1252.9, 1572.9 };
      beta  = { 0, 5E-6, 1.5E-5, 0.02, 0.05, 0.1, 2.96, 3.11, 3.26, 20.7, 22.1, 
      23.5, 34, 36.6, 127, 134, 166.1, 174.1, 286.1, 294.1, 366.1, 374.1 };
      
      delta = 1.265e-3;

      rho = { 0.0000, 7.3032E-1, 2.0992, 3.2419, 3.8133, 3.8734, 3.4510, 2.8342, 
      2.3973, 2.2189, 2.2095, 2.3078, 2.4769, 2.6574, 2.8489, 3.0774, 3.3251, 
      3.6013, 3.9209, 4.2384, 4.5240, 4.8234, 5.1446, 5.4343, 5.6693, 5.8538, 
      6.0082, 6.1852, 6.4578, 6.8818, 7.4595, 8.1641, 8.9778, 9.8244, 1.0544E1, 
      1.1020E1, 1.1301E1, 1.1507E1, 1.1599E1, 1.1458E1, 1.1202E1, 1.0937E1, 
      1.0651E1, 1.0335E1, 1.0025E1, 9.6916, 9.3215, 9.0176, 8.7790, 8.4502, 
      8.0246, 7.6439, 7.3074, 6.9159, 6.5075, 6.1742, 5.8941, 5.5955, 5.2722, 
      4.9620, 4.6809, 4.4273, 4.2036, 4.0025, 3.7957, 3.5501, 3.2763, 2.9786, 
      2.6354, 2.2564, 1.9197, 1.6535, 1.4280, 1.2293, 1.0914, 9.9702E-1, 
      9.0640E-1, 8.2919E-1, 7.7757E-1, 7.3396E-1, 6.9195E-1, 6.6491E-1, 
      6.4252E-1, 6.1313E-1, 5.9026E-1, 5.7573E-1, 5.5665E-1, 5.3679E-1, 
      5.2511E-1, 5.1496E-1, 5.0250E-1, 4.9358E-1, 4.8628E-1, 4.7896E-1, 0.47672};

      temp = {283.6};
      twt = 1.5419e-2; c = 2.4915; tbeta = 5.4693e-1;
      oscE   = {0.15,    0.305  };
      oscW   = {0.14588, 0.29176};
      dka    = 0.05*2;
      kappa  = {1.4345E-1, 1.4327E-1, 1.4254E-1, 1.4105E-1, 1.3949E-1, 1.4005E-1,
        1.4594E-1, 1.6001E-1, 1.8330E-1, 2.1466E-1, 2.5246E-1, 2.9787E-1, 
        3.5823E-1, 4.4795E-1, 5.8492E-1, 7.8208E-1, 1.0360, 1.3171, 1.5668, 
        1.7109, 1.6919, 1.5145, 1.2664, 1.0697, 9.8035E-1, 9.5375E-1, 9.2609E-1, 
        8.7794E-1, 8.2510E-1, 7.8982E-1, 7.8442E-1, 8.0831E-1, 8.5214E-1, 
        9.0378E-1, 9.5285E-1, 9.9272E-1, 1.0204, 1.0353, 1.0379, 1.0291, 1.0103, 
        9.8324E-1, 9.5046E-1, 9.1501E-1, 8.8008E-1, 8.4863E-1, 8.2315E-1, 
        8.0551E-1, 7.9701E-1, 7.9830E-1, 8.0928E-1, 8.2889E-1, 8.5507E-1, 
        8.8492E-1, 9.1516E-1, 9.4272E-1, 9.6535E-1, 9.8197E-1, 9.9266E-1, 
        9.9838E-1, 1.0005, 1.0002, 9.9847E-1, 9.9560E-1, 9.9174E-1, 9.8699E-1, 
        9.8160E-1, 9.7611E-1, 9.7126E-1, 9.6790E-1, 9.6676E-1, 9.6839E-1, 
        9.7306E-1, 9.8080E-1, 9.9141E-1, 1.0046, 1.0198, 1.0366, 1.0542, 1.0718, 
        1.0888, 1.1044, 1.1179, 1.1287, 1.1365, 1.1410, 1.1421, 1.1399, 1.1349, 
        1.1273, 1.1177, 1.1066, 1.0944, 1.0815, 1.0683, 1.0550, 1.0419, 1.0291, 
        1.0169, 1.0051, 9.9361E-1, 9.8229E-1, 9.7093E-1, 9.5950E-1, 9.4813E-1, 
        9.3709E-1, 9.2672E-1, 9.1730E-1, 9.0903E-1, 9.0201E-1, 8.9629E-1, 
        8.9196E-1, 8.8915E-1, 8.8802E-1, 8.8871E-1, 8.9124E-1, 8.9548E-1, 
        9.0120E-1, 9.0812E-1, 9.1602E-1, 9.2473E-1, 9.3417E-1, 9.4427E-1, 
        9.5488E-1, 9.6578E-1, 9.7672E-1, 9.8744E-1, 9.9776E-1, 1.0076, 1.0169, 
        1.0257, 1.0338, 1.0412, 1.0477, 1.0532, 1.0576, 1.0611, 1.0636, 1.0653, 
        1.0662, 1.0661, 1.0650, 1.0630, 1.0601, 1.0566, 1.0526, 1.0484, 1.0440, 
        1.0394, 1.0345, 1.0295, 1.0244, 1.0194, 1.0147, 1.0105, 1.0068, 1.0036, 
        1.0008, 9.9838E-1, 9.9630E-1, 9.9454E-1, 9.9311E-1, 9.9196E-1, 9.9106E-1, 
        9.9031E-1, 9.8964E-1, 9.8902E-1, 9.8842E-1, 9.8785E-1, 9.8730E-1, 
        9.8673E-1, 9.8605E-1, 9.8521E-1, 9.8413E-1, 9.8282E-1, 9.8132E-1, 
        9.7971E-1, 9.7805E-1, 9.7640E-1, 9.7479E-1, 9.7321E-1, 9.7169E-1, 
        9.7024E-1, 9.6893E-1, 9.6780E-1, 9.6695E-1, 9.6642E-1, 9.6628E-1, 
        9.6654E-1, 9.6724E-1, 9.6840E-1, 9.7003E-1, 9.7215E-1, 9.7476E-1, 
        9.7783E-1, 9.8134E-1, 9.8523E-1, 9.8945E-1, 9.9393E-1, 9.9862E-1, 1.0034, 
        1.0083, 1.0132, 1.0181, 1.0228, 1.0272, 1.0314, 1.0352, 1.0385, 1.0414, 
        1.0437, 1.0454, 1.0465, 1.0470, 1.0469, 1.0460, 1.0446, 1.0425, 1.0398, 
        1.0367, 1.0331, 1.0292, 1.0248, 1.0202, 1.0154, 1.0105, 1.0056, 1.0008, 
        9.9609E-1, 9.9154E-1, 9.8716E-1, 9.8301E-1, 9.7914E-1, 9.7563E-1, 
        9.7256E-1, 9.6998E-1, 9.6790E-1, 9.6633E-1, 9.6522E-1, 9.6456E-1, 
        9.6434E-1, 9.6455E-1, 9.6522E-1, 9.6632E-1, 9.6783E-1, 9.6967E-1, 
        9.7177E-1, 9.7405E-1, 9.7648E-1, 9.7901E-1 };
      cfrac  = 0.73194;

      auto oscEnergiesWeights = ranges::view::zip(oscE,oscW);
      auto secondaryScatterInput = std::make_tuple(nss,b7,std::vector<double>(0));


      auto out = leapr( nphon, awr, iel, npr, ncold, aws, lat, alpha, beta, temp, delta, 
                        rho, twt, c, tbeta, oscEnergiesWeights, dka, kappa, cfrac, secondaryScatterInput );

 
      std::vector<double> sab = std::get<0>(out),
      sab_alpha_0 { 6.9151768E3, 6.0129601E3, 3.0356372E3, 1.080692E-2, 
      2.279549E-3, 1.127651E-3, 6.959806E-6, 5.582740E-6, 4.444004E-6, 
      4.798484E-17, 7.664351E-19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      sab_alpha_9 { 4.436529E-1, 4.436539E-1, 4.436560E-1, 4.475444E-1, 
      4.504088E-1, 4.495092E-1, 7.396225E-2, 6.660250E-2, 5.945740E-2, 
      1.866872E-4, 7.062946E-5, 1.954067E-4, 1.572299E-6, 5.743717E-6, 0, 0, 0, 
      0, 0, 0, 0, 0 },
      sab_alpha_12 { 1.321267E-8, 1.321271E-8, 1.321277E-8, 1.334767E-8, 
      1.355301E-8, 1.390300E-8, 4.897118E-8, 5.181986E-8, 5.481936E-8, 
      1.247930E-5, 1.812182E-5, 2.572030E-5, 2.370183E-4, 3.671722E-4, 
      2.127344E-4, 1.129773E-4, 4.161206E-6, 1.468186E-6, 1.456277E-18, 0, 0, 
      1.017557E-16},
      sab_alpha_15 { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
      0, 0, 0};
  
      
      checkPartOfVec( std::get<0>(out), sab_alpha_0,   0*beta.size() );
      checkPartOfVec( std::get<0>(out), sab_alpha_9,   9*beta.size() );
      checkPartOfVec( std::get<0>(out), sab_alpha_12, 12*beta.size() );
      checkPartOfVec( std::get<0>(out), sab_alpha_15, 15*beta.size() );

    } // WHEN
  } // GIVEN


  GIVEN( "Beryllium oxide" ) {
    WHEN( "a secondary scatterer is considered"){

      nphon = 100;
      awr   = 8.93478;  npr = 1; iel = 3; ncold = 0; 
      nss   = 1; b7 = 0; aws   = 15.858; sps = 3.7481; 
      lat   = 1;
      alpha = { .252, .504, .756, 1.008, 1.260, 1.512, 1.764, 2.016, 2.268, 
      2.520, 2.772, 3.024, 3.282, 3.544, 3.813, 4.087, 4.366, 4.652, 4.943, 
      5.241, 5.545, 5.855, 6.172, 6.495, 6.825, 7.162, 7.507, 7.858, 8.217, 
      8.583, 8.957, 9.339, 9.729, 10.13, 10.53, 10.95, 11.37, 11.81, 12.25, 
      12.69, 13.16, 13.63, 14.11, 14.60, 15.10, 15.61, 16.13, 16.66, 17.21, 
      17.76 };
      beta = { 0., .1513, .3025, .4537, .6049, .7561, .9073, 1.059, 1.210, 
      1.361, 1.512, 1.663, 1.815, 1.966, 2.117, 2.268, 2.419, 2.571, 2.722, 
      2.873, 3.024, 3.176, 3.327, 3.478, 3.629, 3.780, 3.932, 4.083, 4.241, 
      4.408, 4.583, 4.766, 4.958, 5.159, 5.371, 5.592, 5.825, 6.069, 6.325, 
      6.593, 6.875, 7.170, 7.480, 7.805, 8.146, 8.504, 8.879, 9.273, 9.686, 
      10.12, 10.57, 11.05, 11.55, 12.07, 12.62, 13.20, 13.81, 14.44, 15.11, 
      15.81, 16.54, 17.31, 18.12, 18.96, 19.85, 20.78, 21.76, 22.78, 23.86, 
      24.99, 26.17, 27.41, 28.71, 30.08, 31.51, 33.01, 34.59, 36.24, 37.98, 
      39.80 };
      temp  = { 296.0 };                                          

      delta = .0016518;
      std::vector<double> rho_Be = { 0, .3, .7, .9, 1., 1.2, 1.6, 2, 2.2, 3, 
      3.5, 4.5, 5.5, 6.8, 8, 9.2, 10.9, 12.9, 15.5, 18.6, 22, 26, 30.5, 35, 39, 
      40, 34, 28, 26, 24.4, 23, 21.3, 19.8, 17, 14.1, 12, 10, 9, 9, 8.5, 7.5, 6, 
      4.6, 3.1, 1.6, 0.5, 0., 0, 4, 15, 38, 52, 70, 105, 165, 230, 200, 170, 
      145, 136, 134, 112, 96, 89, 84, 75, 87, 81, 66, 59, 68, 105, 95, 97, 135, 
      163, 130, 111, 92, 67, 45, 19, 7, 0 },
      rho_O = {0, 0.4, 0.8, 1, 1.4, 2, 2.5, 3.5, 4.8, 6.2, 8.9, 11, 14, 17.2, 
      21.5, 26.5, 34, 40, 46, 58, 60, 93, 110, 129, 141, 142, 125, 101, 93, 92, 
      91, 95, 95, 98, 108, 93, 78, 98, 112, 115, 145, 160, 190, 190, 120, 43, 
      0, 0, 1, 9, 19, 26, 35, 48, 66, 92, 82, 56, 44, 35, 29, 21, 15, 11.5, 9, 
      8, 7, 6, 5.2, 4.5, 5, 5.9, 6, 5, 4, 2.5, 1.8, 1, 0.5, 0.50, 0.2, 0, 0, 0};



      twt    = 0.0;   c = 0.0;   tbeta = 1.0;
      oscE   = { };                                 
      oscW   = { };                            
      dka    = 0.0;
      cfrac  = 0.0;
      auto secondaryScatterInput = std::make_tuple(nss,b7,rho_O);
      auto oscEnergiesWeights = ranges::view::zip(oscE,oscW);

      auto out = leapr( nphon, awr, iel, npr, ncold, aws, lat, alpha, beta, temp, delta, 
      rho_Be, twt, c, tbeta, oscEnergiesWeights, dka, kappa, cfrac, secondaryScatterInput );

      std::vector<double> effectiveTemps1 {596.67224686},
                          effectiveTemps2 {427.92269393},
                          debyeWaller1    {0.4933530515},
                          debyeWaller2    {0.8818971809},
      sab1_alpha_0 { 6.29193E-2, 3.42835E-2, 1.32566E-2, 1.09786E-2, 1.04764E-2, 
      1.15726E-2, 1.32879E-2, 1.49366E-2, 1.83317E-2, 2.27266E-2, 2.73285E-2, 
      2.50529E-2, 1.59443E-2, 1.24747E-2, 9.38148E-3, 5.98628E-3, 4.12817E-3, 
      3.54511E-3, 2.27881E-3, 1.05735E-3, 5.64969E-4, 3.84247E-3, 1.49046E-2, 
      3.26393E-2, 5.41286E-2, 3.62788E-2, 3.04328E-2, 2.10853E-2, 1.66747E-2, 
      1.56121E-2, 1.54433E-2, 1.90942E-2, 2.47544E-2, 1.25443E-2, 1.79039E-3, 
      7.95924E-4, 6.80115E-4, 5.87810E-4, 5.16439E-4, 4.10645E-4, 3.93152E-4, 
      5.85319E-4, 6.02568E-4, 5.13025E-4, 5.20746E-4, 6.66682E-4, 4.46563E-4, 
      2.40055E-4, 1.78008E-4, 8.33936E-5, 1.46131E-5, 1.18548E-5, 1.00354E-5, 
      1.05085E-5, 8.56971E-6, 6.13765E-6, 3.50533E-6, 1.31623E-6, 4.39957E-7, 
      1.58419E-7, 1.18657E-7, 8.32124E-8, 3.99688E-8, 1.32923E-8, 3.53850E-9, 
      1.28276E-9, 6.92649E-10, 2.61334E-10, 6.18222E-11, 1.38797E-11, 
      5.16923E-12, 1.51297E-12, 2.66579E-13, 5.30734E-14, 1.37056E-14, 
      2.11175E-15, 2.98620E-16, 5.42533E-17, 5.88530E-18, 6.95140E-19},

      sab2_alpha_0 { 4.81133E-2, 2.25962E-2, 1.20589E-2, 1.12154E-2, 1.33425E-2, 
      1.67238E-2, 2.04055E-2, 2.65700E-2, 3.19258E-2, 4.46239E-2, 5.69658E-2, 
      5.12138E-2, 3.23843E-2, 2.82217E-2, 2.70557E-2, 2.50986E-2, 2.36814E-2, 
      2.79966E-2, 3.70647E-2, 2.35842E-2, 1.38590E-3, 2.36073E-3, 5.34433E-3, 
      9.17828E-3, 1.35511E-2, 7.41699E-3, 4.75938E-3, 2.83036E-3, 2.11969E-3, 
      1.63906E-3, 1.34857E-3, 1.21429E-3, 8.39224E-4, 6.48088E-4, 5.20179E-4, 
      3.81904E-4, 2.24828E-4, 2.11820E-4, 2.04270E-4, 1.10165E-4, 6.00981E-5, 
      5.21596E-5, 3.56711E-5, 1.81567E-5, 1.08936E-5, 6.61733E-6, 3.78630E-6, 
      2.08283E-6, 1.12873E-6, 6.68799E-7, 3.04509E-7, 1.61864E-7, 6.97864E-8, 
      3.14532E-8, 1.33687E-8, 5.27542E-9, 2.12871E-9, 7.68895E-10, 2.65565E-10, 
      8.52744E-11, 2.60616E-11, 7.47241E-12, 1.97548E-12, 4.92910E-13, 
      1.11466E-13, 2.34034E-14, 4.48506E-15, 7.91720E-16, 1.24476E-16, 
      1.77399E-17, 2.28907E-18, 2.61793E-19, 2.65729E-20, 2.34741E-21, 
      1.83137E-22, 1.23911E-23, 7.12281E-25, 3.53512E-26, 1.46139E-27, 
      5.01966E-29 }, 

      sab1_alpha_40 {1.63025E-2, 1.34122E-2, 1.07498E-2, 1.04380E-2, 1.07038E-2, 
      1.13794E-2, 1.22731E-2, 1.33229E-2, 1.47410E-2, 1.63851E-2, 1.79150E-2, 
      1.81921E-2, 1.73979E-2, 1.73418E-2, 1.74424E-2, 1.75355E-2, 1.78880E-2, 
      1.84952E-2, 1.91087E-2, 1.98182E-2, 2.08557E-2, 2.28477E-2, 2.63222E-2, 
      3.11368E-2, 3.55129E-2, 3.40037E-2, 3.30059E-2, 3.18897E-2, 3.17875E-2, 
      3.26455E-2, 3.44468E-2, 3.76722E-2, 4.09518E-2, 4.07117E-2, 3.92462E-2, 
      3.89610E-2, 3.91337E-2, 3.97233E-2, 4.06877E-2, 4.14868E-2, 4.32466E-2, 
      4.68919E-2, 4.86620E-2, 4.91484E-2, 5.12188E-2, 5.54778E-2, 5.58055E-2, 
      5.47091E-2, 5.50506E-2, 5.48200E-2, 5.40918E-2, 5.45857E-2, 5.46061E-2, 
      5.59416E-2, 5.58157E-2, 5.47575E-2, 5.30733E-2, 5.04379E-2, 4.80906E-2, 
      4.61523E-2, 4.39258E-2, 4.13817E-2, 3.78084E-2, 3.41029E-2, 3.06037E-2, 
      2.73056E-2, 2.38374E-2, 2.02915E-2, 1.69484E-2, 1.39929E-2, 1.12939E-2, 
      8.83676E-3, 6.76262E-3, 5.04682E-3, 3.64352E-3, 2.54962E-3, 1.72591E-3, 
      1.12413E-3, 7.03259E-4, 4.22101E-4}, 

      sab2_alpha_40 {2.73175E-2, 2.61311E-2, 2.61790E-2, 2.75520E-2, 2.95027E-2, 
      3.17372E-2, 3.41880E-2, 3.71261E-2, 4.01019E-2, 4.39533E-2, 4.77031E-2, 
      4.89360E-2, 4.85091E-2, 4.99536E-2, 5.21296E-2, 5.43386E-2, 5.68473E-2, 
      6.03683E-2, 6.42309E-2, 6.50303E-2, 6.44331E-2, 6.63948E-2, 6.86272E-2, 
      7.11246E-2, 7.37068E-2, 7.49702E-2, 7.65116E-2, 7.82267E-2, 7.99989E-2, 
      8.09320E-2, 8.14420E-2, 8.27510E-2, 8.44695E-2, 8.60277E-2, 8.72207E-2, 
      8.77037E-2, 8.77097E-2, 8.82857E-2, 8.87857E-2, 8.84452E-2, 8.78812E-2, 
      8.69544E-2, 8.56027E-2, 8.41192E-2, 8.19525E-2, 7.92605E-2, 7.63569E-2, 
      7.29113E-2, 6.90145E-2, 6.47682E-2, 6.02597E-2, 5.53468E-2, 5.03166E-2, 
      4.51875E-2, 4.00037E-2, 3.48628E-2, 2.98784E-2, 2.52314E-2, 2.08564E-2, 
      1.69053E-2, 1.34237E-2, 1.03966E-2, 7.84376E-3, 5.77762E-3, 4.11937E-3, 
      2.84986E-3, 1.90235E-3, 1.22853E-3, 7.59864E-4, 4.51201E-4, 2.56830E-4, 
      1.39227E-4, 7.17379E-5, 3.48773E-5, 1.60507E-5, 6.93966E-6, 2.79609E-6, 
      1.05342E-6, 3.65857E-7, 1.17493E-7};

      THEN( "the scattering law is correctly returned on a scaled grid" ){
        checkPartOfVec( std::get<0>(out), sab1_alpha_0,     0  );
        checkPartOfVec( std::get<0>(out), sab1_alpha_40,40*beta.size()  );
        REQUIRE( ranges::equal(std::get<1>(out),effectiveTemps1,equal) );
        REQUIRE( ranges::equal(std::get<2>(out),debyeWaller1   ,equal) );

        checkPartOfVec( std::get<4>(out), sab2_alpha_0,     0  );
        checkPartOfVec( std::get<4>(out), sab2_alpha_40,40*beta.size()  );
        REQUIRE( ranges::equal(std::get<5>(out),effectiveTemps2,equal) );
        REQUIRE( ranges::equal(std::get<6>(out),debyeWaller2,   equal) );

      } // THEN

    } // WHEN
  } // GIVEN
  */
//} // TEST CASE

