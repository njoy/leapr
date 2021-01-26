#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "leapr.hpp"
#include <range/v3/all.hpp>
#include "generalTools/testing.h"
#include <variant>
#include "tsl_h2o_tape24.h"
#include "tsl_sic_tape24.h"
#include "tsl_be_tape24.h"
#include "tsl_ortho_h_tape24.h"
#include "tsl_liq_ch4_tape24.h"
#include "input-tsl-SiO2-alpha.h"
#include "input-tsl-BeO.h"
#include "input-tsl-D2O.h"
#include "tsl-SiO2-alpha.h"

using Tabulated         = section::Type<7,4>::TabulatedFunctions;
using CoherentElastic   = section::Type<7,2>::CoherentElastic;
using IncoherentElastic = section::Type<7,2>::IncoherentElastic;

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

      REQUIRE( my_law.LTHR() == good_law.LTHR() );
      REQUIRE( my_law.temperatureDependenceFlag() == 
             good_law.temperatureDependenceFlag() );
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
      REQUIRE(my_law.thermalScatteringValues().size() == 
            good_law.thermalScatteringValues().size());
      for (size_t j = 0; j < my_law.thermalScatteringValues().size(); ++j){
        checkVec(my_law.thermalScatteringValues()[j],
                 good_law.thermalScatteringValues()[j]);
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
      auto my_value   =   my_table.scatteringFunctions()[ibeta];
      auto good_value = good_table.scatteringFunctions()[ibeta];
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

    auto myTempEff   = my_MT4.principalEffectiveTemperature();
    auto goodTempEff = good_MT4.principalEffectiveTemperature();
    checkVec( myTempEff.moderatorTemperatures(), 
            goodTempEff.moderatorTemperatures() );
    checkVec( myTempEff.effectiveTemperatures(), 
            goodTempEff.effectiveTemperatures() );


  } // if mf7 is present
}








TEST_CASE( "LEAPR" ){

  GIVEN( "ENDF-B/VIII Be in BeO input" ){
    njoy::njoy21::lipservice::iRecordStream<char> 
      iss( (std::istringstream(tslBeO)) );

    njoy::njoy21::lipservice::LEAPR leapr(iss);
    nlohmann::json jsonLEAPR(leapr);

    auto args = nlohmann::json::object();
    njoy::LEAPR::LEAPR leaprInstance;

    leaprInstance( jsonLEAPR, std::cout, std::cerr, args );

    njoy::ENDFtk::tree::Tape<std::string> treeTape(njoy::utility::slurpFileToMemory("tape24"));
    long lineNumber = 1;
    int mat = jsonLEAPR["mat"];
    Tape          tape     = treeTape.parse( lineNumber );
    Material      material = tape.material(mat).front();
    file::Type<7> my_MF7   = material.MF(7_c);

    njoy::ENDFtk::tree::Tape<std::string> treeTape2(njoy::utility::slurpFileToMemory("tsl_BeO_tape24"));
    lineNumber = 1;
    Tape          tape2     = treeTape2.parse( lineNumber );
    Material      material2 = tape2.material(mat).front();
    file::Type<7> good_MF7  = material2.MF(7_c);


    check_MF7(my_MF7,good_MF7);

  } // GIVEN

  GIVEN( "Single temperature D in D2O ENDFB/VIII.0 input (coarse alpha beta)" ){
    njoy::njoy21::lipservice::iRecordStream<char> 
      iss( (std::istringstream(tslD2O)) );

    njoy::njoy21::lipservice::LEAPR leapr(iss);
    nlohmann::json jsonLEAPR(leapr);

    auto args = nlohmann::json::object();
    njoy::LEAPR::LEAPR leaprInstance;
    leaprInstance( jsonLEAPR, std::cout, std::cerr, args );

    njoy::ENDFtk::tree::Tape<std::string> treeTape(njoy::utility::slurpFileToMemory("tape24"));
    long lineNumber = 1;
    int mat = jsonLEAPR["mat"];
    Tape          tape     = treeTape.parse( lineNumber );
    Material      material = tape.material(mat).front();
    file::Type<7> my_MF7   = material.MF(7_c);

    njoy::ENDFtk::tree::Tape<std::string> treeTape2(njoy::utility::slurpFileToMemory("tsl_D2O_tape24"));
    lineNumber = 1;
    Tape          tape2     = treeTape2.parse( lineNumber );
    Material      material2 = tape2.material(mat).front();
    file::Type<7> good_MF7  = material2.MF(7_c);

    check_MF7(my_MF7,good_MF7);

  } // GIVEN



  GIVEN( "test SiO2-alpha ENDF-B/VIII.0 input" ){
    njoy::njoy21::lipservice::iRecordStream<char> 
      iss( (std::istringstream(tslSiO2)) );

    njoy::njoy21::lipservice::LEAPR leapr(iss);
    nlohmann::json jsonLEAPR(leapr);

    auto args = nlohmann::json::object();
    njoy::LEAPR::LEAPR leaprInstance;
    leaprInstance( jsonLEAPR, std::cout, std::cerr, args );

    njoy::ENDFtk::tree::Tape<std::string> treeTape(njoy::utility::slurpFileToMemory("tape24"));
    long lineNumber = 1;
    int mat = jsonLEAPR["mat"];
    Tape          tape     = treeTape.parse( lineNumber );
    Material      material = tape.material(mat).front();
    file::Type<7> my_MF7   = material.MF(7_c);

    // Bring in the good (legacy) inputs
    auto begin = sio2_small_alphaBetaGrid.begin(), 
           end = sio2_small_alphaBetaGrid.end();
    lineNumber = 1;
    StructureDivision division(begin,end,lineNumber);
    file::Type<7>     h2o_MF7(division,begin,end,lineNumber);
  
    check_MF7(my_MF7,h2o_MF7);


  } // GIVEN
  
  GIVEN( "H in H2O ENDF-B/VIII.0 input" ) {
    njoy::njoy21::lipservice::iRecordStream<char> iss( std::istringstream(
        "20/\n"                             // Card1
        "title/\n"                          // Card2
        "2 1 100/\n"                        // Card3
        "1 1001 0 0 2e-38/\n"               // Card4
        "0.9991673 20.43608 2 0 0 0/\n"     // Card5
        "1 1 15.85751 3.7939 1\n"           // Card6
        "2 3 1/\n"                          // Card7
        "1e-5 12.0/\n"                      // Card8
        "0.0 1.0 9.0/\n"                    // Card9
        "283.6/\n"                          // Card10 --
        ".001265 119/\n"                    // Card11  |
        "0.000000E+00 4.005204E-01 1.119462E+00 1.698629E+00 2.010136E+00 \n"
        "2.072919E+00 1.898872E+00 1.580587E+00 1.322490E+00 1.206897E+00 \n"
        "1.171350E+00 1.193990E+00 1.252900E+00 1.308722E+00 1.372691E+00 \n"
        "1.458393E+00 1.544730E+00 1.632884E+00 1.737451E+00 1.839542E+00 \n"
        "1.926015E+00 2.027550E+00 2.141893E+00 2.232498E+00 2.288891E+00 \n"
        "2.340017E+00 2.381759E+00 2.383409E+00 2.343025E+00 2.309382E+00 \n"
        "2.315281E+00 2.363291E+00 2.456813E+00 2.608183E+00 2.803355E+00 \n"
        "3.023636E+00 3.287621E+00 3.596284E+00 3.917294E+00 4.268025E+00 \n"
        "4.680325E+00 5.108988E+00 5.509916E+00 5.935525E+00 6.396799E+00 \n"
        "6.827404E+00 7.218567E+00 7.614806E+00 7.963608E+00 8.212350E+00 \n"
        "8.407257E+00 8.565375E+00 8.611701E+00 8.551322E+00 8.491609E+00 \n"
        "8.456424E+00 8.382308E+00 8.252262E+00 8.082789E+00 7.892225E+00 \n"
        "7.709093E+00 7.552676E+00 7.385949E+00 7.188253E+00 6.997005E+00 \n"
        "6.834275E+00 6.658732E+00 6.452670E+00 6.251437E+00 6.078947E+00 \n"
        "5.906580E+00 5.721030E+00 5.545298E+00 5.382871E+00 5.199881E+00 \n"
        "4.998420E+00 4.809347E+00 4.637318E+00 4.460711E+00 4.285237E+00 \n"
        "4.119970E+00 3.957456E+00 3.802720E+00 3.667038E+00 3.537688E+00 \n"
        "3.392823E+00 3.237844E+00 3.083621E+00 2.938483E+00 2.794379E+00 \n"
        "2.618337E+00 2.398675E+00 2.174743E+00 1.973006E+00 1.767542E+00 \n"
        "1.559136E+00 1.388105E+00 1.245415E+00 1.106893E+00 9.899717E-01 \n"
        "9.134314E-01 8.420908E-01 7.698313E-01 7.170319E-01 6.787253E-01 \n"
        "6.379317E-01 6.056904E-01 5.856277E-01 5.618792E-01 5.377198E-01 \n"
        "5.222762E-01 5.089940E-01 4.943599E-01 4.838341E-01 4.732544E-01 \n"
        "4.588599E-01 4.485080E-01 4.417946E-01 4.331451E-01 /\n" // Card12
        "6.9210e-3 3.5910e+0 5.2308e-1/\n"  // Card13
        "2/\n"                              // Card14
        "0.205 0.415/\n"                    // Card15
        "0.15667 0.31333/\n"                // Card16
        "550/\n"                            // Card10 --
        ".001265 119/\n"                    // Card11  |
        "0.000000E+00 4.324438E-01 1.149558E+00 1.832510E+00 2.240700E+00 \n"
        "2.487719E+00 2.709575E+00 2.874102E+00 2.991449E+00 3.119467E+00 \n"
        "3.233615E+00 3.348222E+00 3.487496E+00 3.640314E+00 3.791807E+00 \n"
        "3.963414E+00 4.135150E+00 4.291491E+00 4.451899E+00 4.642881E+00 \n"
        "4.819324E+00 4.985605E+00 5.123753E+00 5.316881E+00 5.490964E+00 \n"
        "5.629008E+00 5.779357E+00 5.934624E+00 6.069339E+00 6.214651E+00 \n"
        "6.328005E+00 6.433464E+00 6.565717E+00 6.626268E+00 6.671826E+00 \n"
        "6.788980E+00 6.783930E+00 6.845387E+00 6.873965E+00 6.874528E+00 \n"
        "6.823557E+00 6.769115E+00 6.716161E+00 6.661474E+00 6.617604E+00 \n"
        "6.522212E+00 6.395003E+00 6.329291E+00 6.214303E+00 6.122996E+00 \n"
        "6.002230E+00 5.893373E+00 5.767771E+00 5.630944E+00 5.495732E+00 \n"
        "5.382765E+00 5.241958E+00 5.089607E+00 4.982199E+00 4.846551E+00 \n"
        "4.738482E+00 4.579191E+00 4.449045E+00 4.332928E+00 4.195822E+00 \n"
        "4.038285E+00 3.910536E+00 3.780957E+00 3.657666E+00 3.521597E+00 \n"
        "3.388630E+00 3.246092E+00 3.096439E+00 2.962325E+00 2.813099E+00 \n"
        "2.667498E+00 2.523798E+00 2.387101E+00 2.249823E+00 2.108490E+00 \n"
        "1.995867E+00 1.880804E+00 1.767768E+00 1.665855E+00 1.574605E+00 \n"
        "1.478388E+00 1.397494E+00 1.322468E+00 1.263067E+00 1.193514E+00 \n"
        "1.135380E+00 1.081620E+00 1.035438E+00 9.859974E-01 9.494125E-01 \n"
        "9.155681E-01 8.800506E-01 8.444502E-01 8.097923E-01 7.877114E-01 \n"
        "7.617647E-01 7.387365E-01 7.203094E-01 7.016640E-01 6.885593E-01 \n"
        "6.740936E-01 6.574666E-01 6.430402E-01 6.322819E-01 6.192831E-01 \n"
        "6.127931E-01 6.031065E-01 5.967796E-01 5.921670E-01 5.861204E-01 \n"
        "5.801695E-01 5.781248E-01 5.760360E-01 5.787275E-01/\n" // Card12
        "2.129e-2 2.2343e+1 5.0871e-1/\n"   // Card13
        "2/\n"                              // Card14
        "0.205 0.415/\n"                    // Card15
        "0.15667 0.31333/\n"                // Card16
        "/"                                 // Card20
      ) );
    njoy::njoy21::lipservice::LEAPR leapr(iss);
    nlohmann::json jsonLEAPR(leapr);

    auto args = nlohmann::json::object();
    njoy::LEAPR::LEAPR leaprInstance;
    leaprInstance( jsonLEAPR, std::cout, std::cerr, args );

    njoy::ENDFtk::tree::Tape<std::string> treeTape(njoy::utility::slurpFileToMemory("tape20"));
    long lineNumber = 1;
    int mat = jsonLEAPR["mat"];
    Tape          tape     = treeTape.parse( lineNumber );
    Material      material = tape.material(mat).front();
    file::Type<7> my_MF7   = material.MF(7_c);

    // Bring in the good (legacy) inputs
    auto begin = h2o.begin(), end = h2o.end();
    lineNumber = 1;
    StructureDivision division(begin,end,lineNumber);
    file::Type<7>     h2o_MF7(division,begin,end,lineNumber);
  
    check_MF7(my_MF7,h2o_MF7);
  } // GIVEN


  GIVEN( "Si in SiC ENDF-B/VIII.0 input" ) {
    njoy::njoy21::lipservice::iRecordStream<char> iss( std::istringstream(
        "20/\n"                             // Card1
        "title/\n"                          // Card2
        "1 1 100/\n"                        // Card3
        "43 143 0 0 1e-75/\n"               // Card4
        "27.844 2.042 1 0 0 0/\n"     // Card5
        "0\n"           // Card6
        "2 3 1/\n"                          // Card7
        "0.5 2.5/\n"                      // Card8
        "1.0 2.0 6.0/\n"                    // Card9
        "300.0/\n"                          // Card10 --
        "0.001 120/\n"                    // Card11  |
        "0.000000E+00 4.450636E-03 1.456132E-02 2.924358E-02 5.227066E-02\n"
        "8.221577E-02 1.220293E-01 1.759206E-01 2.323762E-01 2.994495E-01\n"
        "3.817621E-01 4.804249E-01 6.073396E-01 7.250425E-01 8.596940E-01\n"
        "1.028602E+00 1.225708E+00 1.483551E+00 1.760605E+00 2.073063E+00\n"
        "2.444496E+00 2.896588E+00 3.413903E+00 3.945487E+00 4.500293E+00\n"
        "5.141146E+00 5.841038E+00 6.654291E+00 7.603319E+00 8.658871E+00\n"
        "9.969959E+00 1.168247E+01 1.437138E+01 1.669329E+01 1.772658E+01\n"
        "1.858584E+01 1.951784E+01 2.057091E+01 2.169448E+01 2.289473E+01\n"
        "2.419614E+01 2.552416E+01 2.692603E+01 2.858531E+01 3.087500E+01\n"
        "2.936516E+01 2.460824E+01 2.145302E+01 1.886116E+01 1.639432E+01\n"
        "1.539083E+01 1.569239E+01 1.599049E+01 1.625384E+01 1.655882E+01\n"
        "1.703456E+01 1.269941E+01 6.601566E+00 5.448713E+00 6.157348E+00\n"
        "7.068402E+00 8.408568E+00 1.031724E+01 1.324826E+01 1.894926E+01\n"
        "2.742393E+01 3.172214E+01 2.669816E+01 1.913398E+01 1.464408E+01\n"
        "1.174114E+01 9.469814E+00 7.628035E+00 5.864667E+00 3.898618E+00\n"
        "2.482198E+00 1.713215E+00 6.918037E-01 3.168440E-02 0.000000E+00\n"
        "0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00 0.000000E+00\n"
        "0.000000E+00 0.000000E+00 0.000000E+00 5.172309E-03 3.464289E-01\n"
        "3.641323E+00 1.485216E+01 2.669359E+01 3.211585E+01 2.519812E+01\n"
        "1.278351E+01 5.693553E+00 1.124512E+00 0.000000E+00 2.044555E-02\n"
        "2.178784E-01 8.314270E-01 2.457587E+00 6.139903E+00 6.256002E+00\n"
        "3.370903E+00 2.669282E+00 2.364179E+00 2.200114E+00 2.076445E+00\n"
        "1.994763E+00 1.954955E+00 1.929053E+00 1.937071E+00 2.001562E+00\n"
        "2.196842E+00 2.097044E+00 1.486082E+00 9.119468E-01 4.707854E-01\n"
        "0 0 1.0/\n"                        // Card13
        "0/\n"                              // Card14
        "/"                                 // Card20
      ) );
    njoy::njoy21::lipservice::LEAPR leapr(iss);
    nlohmann::json jsonLEAPR(leapr);

    auto args = nlohmann::json::object();
    njoy::LEAPR::LEAPR leaprInstance;
    leaprInstance( jsonLEAPR, std::cout, std::cerr, args );

    njoy::ENDFtk::tree::Tape<std::string> treeTape(njoy::utility::slurpFileToMemory("tape20"));
    long lineNumber = 1;
    int mat = jsonLEAPR["mat"];
    Tape          tape     = treeTape.parse( lineNumber );
    Material      material = tape.material(mat).front();
    file::Type<7> my_MF7   = material.MF(7_c);

    auto begin = sic.begin(), end = sic.end();
    lineNumber = 1;
    StructureDivision division(begin,end,lineNumber);
    njoy::ENDFtk::file::Type<7> sic_MF7(division,begin,end,lineNumber);
  
    check_MF7(my_MF7,sic_MF7);

      
  } // GIVEN

    


  GIVEN( "Be-metal ENDF-B/VIII.0 input" ) {
    WHEN( "Multiple temperatures considered" ){
      njoy::njoy21::lipservice::iRecordStream<char> iss( std::istringstream(
        "20/\n"                             // Card1
        "title/\n"                          // Card2
        "3 1 100/\n"                        // Card3
        "26 126 0 0 1e-75/\n"               // Card4
        "8.93478 6.153875 1 2 0 0/\n"     // Card5
        "0\n"           // Card6
        "2 2 1/\n"                          // Card7
        "10 20/\n"                      // Card8
        "15 25/\n"                    // Card9
        "296.0/\n"                          // Card10 --
        "0.00069552 127/\n"                    // Card11  |
        "0.0000E+00 7.2477E-04 3.7084E-03 8.0087E-03 1.0642E-02 1.5897E-02 \n"
        "2.7372E-02 4.1843E-02 5.0214E-02 6.5036E-02 8.3674E-02 9.9329E-02 \n"
        "1.1977E-01 1.4296E-01 1.6484E-01 1.8945E-01 2.1887E-01 2.3537E-01 \n"
        "2.6166E-01 3.0003E-01 3.4054E-01 3.8728E-01 4.2481E-01 4.7598E-01 \n"
        "5.1890E-01 5.7400E-01 6.2970E-01 6.5754E-01 7.2042E-01 7.9118E-01 \n"
        "8.6756E-01 9.2948E-01 1.0030E+00 1.1163E+00 1.2048E+00 1.2870E+00 \n"
        "1.4139E+00 1.5249E+00 1.6221E+00 1.7638E+00 1.8924E+00 2.0388E+00 \n"
        "2.2056E+00 2.3709E+00 2.5558E+00 2.7595E+00 3.0108E+00 3.2603E+00 \n"
        "3.5066E+00 3.7442E+00 4.0067E+00 4.3677E+00 4.7164E+00 5.0820E+00 \n"
        "5.5881E+00 6.0898E+00 6.5510E+00 7.0877E+00 7.5931E+00 8.0736E+00 \n"
        "8.6232E+00 9.2283E+00 9.9334E+00 1.0613E+01 1.1278E+01 1.1973E+01 \n"
        "1.2784E+01 1.3744E+01 1.4739E+01 1.5918E+01 1.7654E+01 1.9834E+01 \n"
        "2.1455E+01 2.2574E+01 2.3744E+01 2.4900E+01 2.6227E+01 2.7931E+01 \n"
        "2.9747E+01 2.9884E+01 2.7358E+01 2.4817E+01 2.3690E+01 2.3242E+01 \n"
        "2.3624E+01 2.3473E+01 2.2368E+01 2.1447E+01 2.0724E+01 2.1121E+01 \n"
        "2.4240E+01 2.7607E+01 2.7643E+01 2.5431E+01 2.3755E+01 2.3377E+01 \n"
        "2.3410E+01 2.3504E+01 2.3647E+01 2.3681E+01 2.3805E+01 2.3714E+01 \n"
        "2.3385E+01 2.3050E+01 2.2244E+01 2.1008E+01 1.9536E+01 1.8341E+01 \n"
        "1.8075E+01 1.8606E+01 1.9599E+01 2.1037E+01 2.3193E+01 2.4016E+01 \n"
        "2.3573E+01 2.5664E+01 3.0187E+01 3.1256E+01 2.7257E+01 2.2765E+01 \n"
        "1.4893E+01 6.8192E+00 3.8444E+00 2.4718E+00 1.3358E+00 3.5968E-01 0\n"
        "0 0 1.0/\n"                        // Card13
        "0/\n"                              // Card14
        "-400\n" 
        "-500\n" 
        "/"                                 // Card20
      ) );
  
      njoy::njoy21::lipservice::LEAPR leapr(iss);
      nlohmann::json jsonLEAPR(leapr);
  
      auto args = nlohmann::json::object();
      njoy::LEAPR::LEAPR leaprInstance;
      leaprInstance( jsonLEAPR, std::cout, std::cerr, args );
  
      njoy::ENDFtk::tree::Tape<std::string> treeTape(njoy::utility::slurpFileToMemory("tape20"));
      long lineNumber = 1;
      int mat = jsonLEAPR["mat"];
      Tape          tape     = treeTape.parse( lineNumber );
      Material      material = tape.material(mat).front();
      file::Type<7> my_MF7   = material.MF(7_c);

      auto begin = be.begin(), end = be.end();
      lineNumber = 1;
      StructureDivision division(begin,end,lineNumber);
      njoy::ENDFtk::file::Type<7> be_MF7(division,begin,end,lineNumber);
  
      check_MF7(my_MF7,be_MF7);
      
    } // WHEN

  } // GIVEN


  GIVEN( "Liquid CH4 ENDF-B/VIII.0 input" ) {
    WHEN( "" ) {
    njoy::njoy21::lipservice::iRecordStream<char> iss( std::istringstream(
        "20/\n"                             // Card1
        "title/\n"                          // Card2
        "3 1 100/\n"                        // Card3
        "33 1001 0 0 1e-75/\n"               // Card4
        "0.99917 20.478 4 0 0 0/\n"     // Card5
        "1 1 11.898 4.7392 1\n"           // Card6
        "4 2 0/\n"                          // Card7
        "2.0 4.0 6.0 8.0 /\n"                      // Card8
        "0.0 10.0 /\n"                    // Card9
        "90.0/\n"                          // Card10 --
        "0.0004 45/\n"                    // Card11  |
        "0. .004 .018 .034 .050 .068 .087 .109\n"
        ".133 .156 .178 .203 .223 .243 .261 .277 .291 .299 .298 .288 \n"
        ".278 .267 .253 .237 .219 .202 .186 .173 .161 .150 .138 .118 \n"
        ".097 .078 .062 .047 .034 .023 .017 .013 .01 .008 .006 .003 0/\n"
        "0.015 200.0 0.305/\n"                        // Card13
        "4/\n"                              // Card14
        "0.162 0.190 0.361 0.374/\n"
        "0.308 0.186 0.042 0.144/\n"
        "-100\n" 
        "-110\n" 
        "/"                                 // Card20
      ) );

      njoy::njoy21::lipservice::LEAPR leapr(iss);
      nlohmann::json jsonLEAPR(leapr);

      auto args = nlohmann::json::object();
      njoy::LEAPR::LEAPR leaprInstance;
      leaprInstance( jsonLEAPR, std::cout, std::cerr, args );


      njoy::ENDFtk::tree::Tape<std::string> treeTape(njoy::utility::slurpFileToMemory("tape20"));
      long lineNumber = 1;
      int mat = jsonLEAPR["mat"];
      Tape          tape     = treeTape.parse( lineNumber );
      Material      material = tape.material(mat).front();
      file::Type<7> my_MF7   = material.MF(7_c);

      // Bring in the good (legacy) inputs
      auto begin = ch4.begin(), end = ch4.end();
      lineNumber = 1;
      StructureDivision division(begin,end,lineNumber);
      file::Type<7>     h2o_MF7(division,begin,end,lineNumber);
  
      check_MF7(my_MF7,h2o_MF7);

      
    } // WHEN

  } // GIVEN

} // TEST CASE

