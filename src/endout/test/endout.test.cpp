#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "endout/endout.h"
#include "generalTools/testing.h"
#include "coher/coher.h"


TEST_CASE( "finalize the debye-waller coefficient output" ){
  GIVEN( "No secondary scatterers considered" ){
  } // GIVEN
  GIVEN( "Secondary scatterers considered" ){
    WHEN( "SCT approximation (BeO test 23)" ){
        std::cout.precision(15);

        int numSecondaryScatterers = 1, secondaryScatterType = 0;
        double awr = 8.9347799999999999, aws = 15.858000000000001;
        std::vector<double> 

        temps { 296.0, 400.0, 500.0, 600.0, 700.0, 800.0, 1000.0, 1200.0 },
        correctDWPIX {2.18024395, 2.72608260, 3.28130480, 3.85303452, 4.43480683, 
              5.02310000, 6.21177434, 7.41006538},
        correctDWP1 {2.16475875, 2.59570940, 3.04447306, 3.51490437, 3.99970405, 
              4.49435669, 5.50273307, 6.5268008},
        dwpix { 0.88189718, 1.49011626, 2.24201096, 3.15918678, 4.24222697, 
                5.49139855, 8.48861475, 12.1513474 },
        dwp1  { 0.49335305, 0.79941568, 1.17203004, 1.62375823, 2.15567111, 
                2.76830621, 4.23677157, 6.03029012 };

        scaleDebyeWallerCoefficients(numSecondaryScatterers, secondaryScatterType, 
                                     dwpix, dwp1, temps, awr, aws );
        REQUIRE( ranges::equal(dwpix,correctDWPIX,equal) );
        REQUIRE( ranges::equal(dwp1,correctDWP1,equal) );


    } // WHEN
    WHEN( "Secondary scatterers use free gas approximation" ){
      THEN( "DWP1 is not changed and DWPIX is scaled by awr" ){
        int numSecondaryScatterers = 1, secondaryScatterType = 1;
        std::vector<double> 
        dwp1 (5,0.0),
        dwpix { 0.27366867, 0.27913484, 0.43494593, 0.86192089, 1.44790731 },
        temps { 296.0, 300.0, 400.0, 600.0, 800.0 },
        correctDWPIX{10.73794694, 10.80639075, 12.62883113, 16.68414797, 21.02028737 },
        correctDWP1 (5,0.0);
        
        double awr = 0.99917, aws = 15.85316; 
        scaleDebyeWallerCoefficients(numSecondaryScatterers, secondaryScatterType, 
                                     dwpix, dwp1, temps, awr, aws );
        REQUIRE( ranges::equal(dwpix,correctDWPIX,equal) );
        REQUIRE( ranges::equal(dwp1,correctDWP1,equal) );

      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE


  /*

TEST_CASE( "processing coherent elastic scattering data" ){
  GIVEN( "high tolerance and eight temperatures" ){
    std::vector<double> bragg (10000,0.0);
    double maxEnergy = 5.0;
    int iel = 3, npr = 1;
    auto out = coher(iel,npr,bragg,maxEnergy);
    int numEdges = int(0.5*std::get<1>(out));

    double tol = 1e-2;
    std::vector<double> 
    temps { 100.0, 900.0, 2000.0 },
    dwpix { 1.3891141273070633, 5.6158501929886118, 12.242592321730232 },
    dwp1  { 1.5294981884837682, 4.9959938487181246, 10.689571188951792 };

    int numSecondaryScatterers = 1, secondaryScatterType = 0;

    auto totalSCR = processCoherentElastic(bragg,dwpix,dwp1,numSecondaryScatterers,secondaryScatterType,numEdges,tol,temps);
    std::vector<std::vector<double>> correctSCR { 
    {1.061174E-3, 0.000000, 3.754368E-3, 0.01864404, 4.244696E-3, 0.03053734, 
     4.815542E-3, 0.04661845, 7.999064E-3, 0.05933238, 9.550566E-3, 0.05933238, 
     0.01126310, 0.1005425, 0.01232428, 0.1005425, 0.01330493, 0.1486804, 0.01501747, 
     0.1574092, 0.01550780, 0.1923674, 0.01607865, 0.2006081, 0.01697878, 0.2008283, 
     0.01926217, 0.2085001, 0.02073315, 0.2087926, 0.02081367, 0.2087926, 0.02456804, 
     0.2419633, 0.02628058, 0.2543204, 0.02652935, 0.2543204, 0.02734175, 0.2661549, 
     0.02824189, 0.2671144, 0.03028372, 0.2960110, 0.03052527, 0.3074238, 0.03199626, 
     0.3076443, 0.03378931, 0.3285054, 0.03485049, 0.3285054, 0.03583114, 0.3799438, 
     0.03779246, 0.3799438, 0.03803401, 0.3995157, 0.03820227, 0.4027672, 0.04154682, 
     0.4258681, 0.04195663, 0.4304213, 0.04325936, 0.4307764, 0.04333988, 0.4307764, 
     0.04505242, 0.4476931, 0.04611359, 0.4476931, 0.04880679, 0.4556435, 0.04929711, 
     0.4717409, 0.04946537, 0.4877951, 0.04986796, 0.4954783, 0.05076810, 0.4961058, 
     0.05199753, 0.5461104, 0.05460298, 0.5461104, 0.05575190, 0.5855024, 0.06006989, 
     0.5888576, 0.06031866, 0.5926588, 0.06431459, 0.6025879, 0.06701500, 0.6055862, 
     0.06791514, 0.6249048, 0.07133300, 0.6306709, 0.07158177, 0.6306709, 0.07166951, 
     0.6537167, 0.07533614, 0.6874693, 0.07704868, 0.7621183, 0.08293261, 0.7976725, 
     0.08595510, 0.8174247, 0.09385920, 0.8556731, 0.09827216, 0.8618886, 0.1009726, 
     0.9107558, 0.1061174, 0.9256253, 0.1120674, 0.9997089, 0.1183467, 0.9999859, 
     0.1197444, 1.011994, 0.1211349, 1.033715, 0.1278973, 1.033715, 0.1279850, 
     1.036807, 0.1284021, 1.067964, 0.1351573, 1.093864, 0.1399067, 1.121800, 
     0.1434195, 1.176430, 0.1506651, 1.186621, 0.1511698, 1.191693, 0.1521360, 
     1.191894, 0.1528091, 1.214293, 0.1614378, 1.216558, 0.1616866, 1.216558, 
     0.1617743, 1.218816, 0.1621914, 1.236427, 0.1678265, 1.266505, 0.1771139, 
     1.268519, 0.1772088, 1.280743, 0.1790896, 1.280822, 0.1793384, 1.332283, 
     0.1884720, 1.354000, 0.1943559, 1.383420, 0.2047777, 1.391446, 0.2056190, 
     1.405659, 0.2079901, 1.445747, 0.2184119, 1.448597, 0.2192532, 1.451427, 
     0.2211124, 1.460524, 0.2222613, 1.465383, 0.2230076, 1.489431, 0.2342707, 
     1.490685, 0.2359616, 1.503516, 0.2387642, 1.541698, 0.2515427, 1.542770, 
     0.2518792, 1.558103, 0.2537816, 1.574731, 0.2668089, 1.584011, 0.2680600, 
     1.585955, 0.2716606, 1.625938, 0.2853320, 1.628450, 0.2863271, 1.630931, 
     0.2866780, 1.662116, 0.3013445, 1.662864, 0.3018493, 1.671250, 0.3054499, 
     1.677739, 0.3066793, 1.712747, 0.3224803, 1.716599, 0.3236292, 1.719789, 
     0.3243755, 1.721024, 0.3262347, 1.721323, 0.3317304, 1.725274, 0.3404686, 
     1.736913, 0.3438204, 1.767960, 0.3633025, 1.773390, 0.3655198, 1.773922, 
     0.3667492, 1.791408, 0.3776097, 1.797173, 0.3830838, 1.817541, 0.4023760, 
     1.817900, 0.4024709, 1.820436, 0.4038903, 1.836881, 0.4168731, 1.850971, 
     0.4244696, 1.868511, 0.4470186, 1.870094, 0.4477649, 1.870616, 0.4482697, 
     1.877715, 0.4556246, 1.878501, 0.4582589, 1.882708, 0.4679778, 1.904490, 
     0.4914269, 1.904490, 0.4920710, 1.905458, 0.4924808, 1.910180, 0.5017671, 
     1.918583, 0.5136082, 1.932640, 0.5393108, 1.933347, 0.5398888, 1.933486, 
     0.5402336, 1.933829, 0.5428680, 1.936544, 0.5469589, 1.936608, 0.5473976, 
     1.944574, 0.5613611, 1.956466, 0.5894358, 1.958357, 0.5924500, 1.958747, 
     0.5926265, 1.960019, 0.5951504, 1.963225, 0.6112363, 1.973420, 0.6419970, 
     1.973695, 0.6440388, 1.974566, 0.6450256, 1.981802, 0.6632338, 1.989518, 
     0.6965183, 1.989766, 0.6970231, 1.992757, 0.7173537, 1.998876, 0.7533097, 
     1.999046, 0.7538878, 1.999924, 0.7570930, 2.002407, 0.7735959, 2.006826, 
     0.8122896, 2.007340, 0.8177121, 2.009315, 0.8319605, 2.012973, 0.8735889, 
     2.012989, 0.8749638, 2.013005, 0.8755213, 2.013082, 0.8767280, 2.013083, 
     0.8770129, 2.013160, 0.8772328, 2.013267, 0.8787482, 2.013808, 0.8890451, 
     2.014066, 0.8920303, 2.014080, 0.8924474, 2.017160, 0.9370828, 2.017161, 
     0.9373027, 2.017223, 0.9411675, 2.017768, 0.9525173, 2.018008, 0.9550566, 
     2.020188, 1.003863, 2.020221, 1.007517, 2.020241, 1.008607, 2.020332, 
     1.013875, 2.020404, 1.015127, 2.020590, 1.019788, 2.022263, 1.071442, 
     2.022289, 1.072240, 2.022307, 1.076411, 2.022367, 1.079858, 2.022474, 5.0, 2.026704 },
    {0.0, 0.01759756, 0.02873892, 0.04367158, 0.05491318, 0.05491318, 0.08956622, 
     0.08956622, 0.1287927, 0.1357207, 0.1632581, 0.1696926, 0.1698622, 0.1755662, 
     0.1757789, 0.1757789, 0.1985081, 0.2067552, 0.2067552, 0.2145257, 0.2151470, 
     0.2332805, 0.2404159, 0.2405507, 0.2529543, 0.2529543, 0.2825927, 0.2825927, 
     0.2934940, 0.2953004, 0.3074904, 0.3098779, 0.3100604, 0.3100604, 0.3185183, 
     0.3185183, 0.3222702, 0.3298097, 0.3373094, 0.3408765, 0.3411638, 0.3636310, 
     0.3636310, 0.3803365, 0.3816680, 0.3831706, 0.3868616, 0.3879308, 0.3947252, 
     0.3966492, 0.3966492, 0.4042995, 0.4148894, 0.4377015, 0.4476261, 0.4528892, 
     0.4619137, 0.4632839, 0.4736185, 0.4765238, 0.4897325, 0.4897773, 0.4916798, 
     0.4950482, 0.4950482, 0.4954798, 0.4998003, 0.5030374, 0.5062828, 0.5122954, 
     0.5132987, 0.5137942, 0.5138135, 0.5159471, 0.5161360, 0.5161360, 0.5163234, 
     0.5177755, 0.5200495, 0.5201815, 0.5209814, 0.5209865, 0.5242455, 0.5254405, 
     0.5269193, 0.5272630, 0.5278637, 0.5294974, 0.5295963, 0.5296933, 0.5299963, 
     0.5301552, 0.5309331, 0.5309672, 0.5313072, 0.5322763, 0.5322987, 0.5326167, 
     0.5329517, 0.5331047, 0.5331361, 0.5337479, 0.5337790, 0.5338093, 0.5341880, 
     0.5341953, 0.5342759, 0.5343349, 0.5346474, 0.5346744, 0.5346963, 0.5347047, 
     0.5347067, 0.5347307, 0.5347925, 0.5349490, 0.5349693, 0.5349712, 0.5350331, 
     0.5350504, 0.5351065, 0.5351072, 0.5351124, 0.5351453, 0.5351684, 0.5351940, 
     0.5351956, 0.5351961, 0.5352033, 0.5352040, 0.5352077, 0.5352239, 0.5352239, 
     0.5352244, 0.5352268, 0.5352306, 0.5352358, 0.5352359, 0.5352360, 0.5352360, 
     0.5352367, 0.5352367, 0.5352385, 0.5352406, 0.5352408, 0.5352408, 0.5352410, 
     0.5352413, 0.5352421, 0.5352422, 0.5352422, 0.5352426, 0.5352428, 0.5352428, 
     0.5352429, 0.5352430, 0.5352430, 0.5352430, 0.5352430, 0.5352431, 0.5352431, 
     0.5352431, 0.5352431, 0.5352431, 0.5352431, 0.5352431, 0.5352431, 0.5352431, 
     0.5352431, 0.5352431, 0.5352431, 0.5352431, 0.5352431, 0.5352431, 0.5352431, 
     0.5352431, 0.5352431, 0.5352431, 0.5352431, 0.5352431, 0.5352431, 0.5352431, 
     0.5352431, 0.5352431, 0.5352431, 0.5352431, 0.5352431, 0.5352431, 0.5352431},
    {0.0, 0.01604264, 0.02607757, 0.03933944, 0.04856999, 0.04856999, 0.07482497, 
     0.07482497, 0.1030868, 0.1078720, 0.1266638, 0.1309935, 0.1311051, 0.1346536, 
     0.1347812, 0.1347812, 0.1471884, 0.1515043, 0.1515043, 0.1554657, 0.1557755, 
     0.1643737, 0.1677370, 0.1677982, 0.1731928, 0.1731928, 0.1854507, 0.1854507, 
     0.1897211, 0.1904258, 0.1948050, 0.1956541, 0.1957170, 0.1957170, 0.1985040, 
     0.1985040, 0.1996311, 0.2018688, 0.2040854, 0.2051294, 0.2052116, 0.2114505, 
     0.2114505, 0.2156795, 0.2159826, 0.2163225, 0.2170791, 0.2172842, 0.2185588, 
     0.2188906, 0.2188906, 0.2201989, 0.2218535, 0.2252705, 0.2265565, 0.2271895, 
     0.2280828, 0.2282045, 0.2290630, 0.2292756, 0.2301105, 0.2301129, 0.2302124, 
     0.2303827, 0.2303827, 0.2304011, 0.2305837, 0.2306995, 0.2308028, 0.2309783, 
     0.2310028, 0.2310148, 0.2310152, 0.2310646, 0.2310682, 0.2310682, 0.2310716, 
     0.2310983, 0.2311347, 0.2311364, 0.2311465, 0.2311466, 0.2311859, 0.2311974, 
     0.2312097, 0.2312119, 0.2312157, 0.2312254, 0.2312258, 0.2312263, 0.2312276, 
     0.2312282, 0.2312314, 0.2312315, 0.2312325, 0.2312352, 0.2312353, 0.2312359, 
     0.2312366, 0.2312368, 0.2312368, 0.2312376, 0.2312376, 0.2312376, 0.2312380, 
     0.2312380, 0.2312380, 0.2312381, 0.2312382, 0.2312382, 0.2312382, 0.2312382, 
     0.2312382, 0.2312382, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 
     0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 
     0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 
     0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 
     0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 
     0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 
     0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 
     0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 
     0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 
     0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 
     0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383, 0.2312383 }};

    for ( size_t t = 0; t < temps.size(); ++t ){
      REQUIRE( ranges::equal( totalSCR[t], correctSCR[t], equal ) );
    }

  } // GIVEN

  GIVEN( "high tolerance and eight temperatures" ){
    std::vector<double> bragg (10000,0.0);
    double maxEnergy = 5.0;
    int iel = 3, npr = 1;
    auto out = coher(iel,npr,bragg,maxEnergy);
    int numEdges = int(0.5*std::get<1>(out));

    double tol = 0.9;
    std::vector<double> 
    temps { 296.0, 400.0, 500.0, 600.0, 700.0, 800.0, 1000.0, 1200.0 },
    dwpix {2.18024395, 2.72608260, 3.28130480, 3.85303452, 4.43480683, 5.02310000, 6.21177434, 7.41006538},
    dwp1  {2.16475875, 2.59570940, 3.04447306, 3.51490437, 3.99970405, 4.49435669, 5.50273307, 6.5268008};

    int numSecondaryScatterers = 1, secondaryScatterType = 0;

    auto totalSCR = processCoherentElastic(bragg,dwpix,dwp1,numSecondaryScatterers,secondaryScatterType,numEdges,tol,temps);
    std::vector<std::vector<double>> correctSCR { 
    {1.061174E-3, 0, 3.754368E-3, 0.01844542, 4.244696E-3, 0.03019557, 
     4.815542E-3, 0.04605728, 7.999064E-3, 0.05848436, 9.550566E-3, 0.05848436, 
     0.01126310, 0.09839143, 0.01232428, 0.09839143, 0.01330493, 0.1447364, 
     0.01501747, 0.1530992, 5.0, 1.351837 },
    {0, 0.01831063, 0.02996375, 0.04567693, 0.05791133, 0.05791133, 0.09694990, 
        0.09694990, 0.1421058, 0.1502268, 9.754622},
    {0, 0.01817311, 0.02972733, 0.04528931, 0.05732877, 0.05732877, 0.09549435, 
        0.09549435, 0.1394598, 0.1473396, 9.457270},
    {0, 0.01803146, 0.02948390, 0.04489046, 0.05673086, 0.05673086, 0.09401090, 
        0.09401090, 0.1367739, 0.1444109, 9.158226},
    {0, 0.01788763, 0.02923684, 0.04448595, 0.05612603, 0.05612603, 0.09252107, 
        0.09252107, 0.1340874, 0.1414836, 8.861996},
    {0, 0.01774276, 0.02898811, 0.04407901, 0.05551915, 0.05551915, 0.09103709, 
        0.09103709, 0.1314227, 0.1385822, 8.571053},
    {0, 0.01745246, 0.02849001, 0.04326493, 0.05430996, 0.05430996, 0.08811286, 
        0.08811286, 0.1262051, 0.1329074, 8.009906},
    {0, 0.01716365, 0.02799491, 0.04245695, 0.05311619, 0.05311619, 0.08526853, 
        0.08526853, 0.1211735, 0.1274430, 7.479787}};

    for ( size_t t = 0; t < temps.size(); ++t ){
      REQUIRE( ranges::equal( totalSCR[t], correctSCR[t], equal ) );
    }

  } // GIVEN
  GIVEN( "low (normal) tolerance and a single temperature" ){
      THEN( "the SCR vector is returned correctly" ){
        std::vector<double> bragg (10000,0.0);
        double maxEnergy = 5.0;
        int iel = 3, npr = 1;
        auto out = coher(iel,npr,bragg,maxEnergy);
        int numEdges = int(0.5*std::get<1>(out));

        double tol = 9e-8;
        std::vector<double> 
        temps { 296.0 },//, 400.0, 500.0, 600.0, 700.0, 800.0, 1000.0, 1200.0 },
        dwpix {2.18024395 },//, 2.72608260, 3.28130480, 3.85303452, 4.43480683, 5.02310000, 6.21177434, 7.41006538},
        dwp1  {2.16475875 };//, 2.59570940, 3.04447306, 3.51490437, 3.99970405, 4.49435669, 5.50273307, 6.5268008};

        int numSecondaryScatterers = 1, secondaryScatterType = 0;
        std::vector<std::vector<double>> correctSCR { { 
        1.061174E-3, 0.000000, 3.754368E-3, 0.01844542, 4.244696E-3, 0.03019557, 
        4.815542E-3, 0.04605728, 7.999064E-3, 0.05848436, 9.550566E-3, 0.05848436, 
        0.01126310, 0.09839143, 0.01232428, 0.09839143, 0.01330493, 0.1447364, 
        0.01501747, 0.1530992, 0.01550780, 0.1865445, 0.01607865, 0.1944157, 
        0.01697878, 0.1946256, 0.01926217, 0.2018871, 0.02073315, 0.2021628, 
        0.02081367, 0.2021628, 0.02456804, 0.2330883, 0.02628058, 0.2445528, 
        0.02652935, 0.2445528, 0.02734175, 0.2554993, 0.02824189, 0.2563846, 
        0.03028372, 0.2828894, 0.03052527, 0.2933505, 0.03199626, 0.2935517, 
        0.03378931, 0.3124959, 0.03485049, 0.3124959, 0.03583114, 0.3589361, 
        0.03779246, 0.3589361, 0.03803401, 0.3764955, 0.03820227, 0.3794114, 
        0.04154682, 0.3999303, 0.04195663, 0.4039698, 0.04325936, 0.4042837, 
        0.04333988, 0.4042837, 0.04505242, 0.4191601, 0.04611359, 0.4191601, 
        0.04880679, 0.4260771, 0.04929711, 0.4400626, 0.04946537, 0.4540039, 
        0.04986796, 0.4606684, 0.05076810, 0.4612112, 0.05199753, 0.5043222, 
        0.05460298, 0.5043222, 0.05575190, 0.5379218, 0.06006989, 0.5407486, 
        0.06031866, 0.5439489, 0.06431459, 0.5522136, 0.06701500, 0.5546901, 
        0.06791514, 0.5706061, 0.07133300, 0.5753105, 0.07158177, 0.5753105, 
        0.07166951, 0.5940949, 0.07533614, 0.6213200, 0.07704868, 0.6812390, 
        0.08293261, 0.7093026, 0.08595510, 0.7247594, 0.09385920, 0.7540229, 
        0.09827216, 0.7587189, 0.1009726, 0.7953557, 0.1061174, 0.8063413, 
        0.1120674, 0.8601529, 0.1183467, 0.8603506, 0.1197444, 0.8688838, 0.1211349, 
        0.8842583, 0.1278973, 0.8842583, 0.1279850, 0.8864049, 0.1284021, 0.9080055, 
        0.1351573, 0.9256194, 0.1399067, 0.9443617, 0.1434195, 0.9806480, 0.1506651, 
        0.9872785, 0.1511698, 0.9905737, 0.1521360, 0.9907037, 0.1528091, 1.005188, 
        0.1614378, 1.006617, 0.1616866, 1.006617, 0.1617743, 1.008041, 0.1621914, 
        1.019128, 0.1678265, 1.037763, 0.1771139, 1.038978, 0.1772088, 1.046352, 
        0.1790896, 1.046399, 0.1793384, 1.077251, 0.1884720, 1.089936, 0.1943559, 
        1.106835, 0.2047777, 1.111309, 0.2056190, 1.119215, 0.2079901, 1.141363, 
        0.2184119, 1.142891, 0.2192532, 1.144405, 0.2211124, 1.149246, 0.2222613, 
        1.151824, 0.2230076, 1.164552, 0.2342707, 1.165195, 0.2359616, 1.171740, 
        0.2387642, 1.191062, 0.2515427, 1.191585, 0.2518792, 1.199059, 0.2537816, 
        1.207120, 0.2668089, 1.211455, 0.2680600, 1.212360, 0.2716606, 1.230781, 
        0.2853320, 1.231894, 0.2863271, 1.232990, 0.2866780, 1.246755, 0.3013445, 
        1.247072, 0.3018493, 1.250616, 0.3054499, 1.253331, 0.3066793, 1.267926, 
        0.3224803, 1.269461, 0.3236292, 1.270729, 0.3243755, 1.271218, 0.3262347, 
        1.271336, 0.3317304, 1.272870, 0.3404686, 1.277276, 0.3438204, 1.288919, 
        0.3633025, 1.290845, 0.3655198, 1.291032, 0.3667492, 1.297174, 0.3776097, 
        1.299137, 0.3830838, 1.305966, 0.4023760, 1.306080, 0.4024709, 1.306884, 
        0.4038903, 1.312080, 0.4168731, 1.316370, 0.4244696, 1.321595, 0.4470186, 
        1.322037, 0.4477649, 1.322183, 0.4482697, 1.324159, 0.4556246, 1.324373, 
        0.4582589, 1.325511, 0.4679778, 1.331243, 0.4914269, 1.331243, 0.4920710, 
        1.331481, 0.4924808, 1.332640, 0.5017671, 1.334648, 0.5136082, 1.337895, 
        0.5393108, 1.338047, 0.5398888, 1.338077, 0.5402336, 1.338151, 0.5428680, 
        1.338728, 0.5469589, 1.338741, 0.5473976, 1.340412, 0.5613611, 1.342810, 
        0.5894358, 1.343162, 0.5924500, 1.343234, 0.5926265, 1.343468, 0.5951504, 
        1.344055, 0.6112363, 1.345838, 0.6419970, 1.345882, 0.6440388, 1.346021, 
        0.6450256, 1.347170, 0.6632338, 1.348333, 0.6965183, 1.348367, 0.6970231, 
        1.348776, 0.7173537, 1.349567, 0.7533097, 1.349587, 0.7538878, 1.349689, 
        0.7570930, 1.349975, 0.7735959, 1.350462, 0.8122896, 1.350512, 0.8177121, 
        1.350704, 0.8319605, 1.351045, 0.8735889, 1.351046, 0.8749638, 1.351047, 
        0.8755213, 1.351054, 0.8767280, 1.351054, 0.8770129, 1.351060, 0.8772328, 
        1.351069, 0.8787482, 1.351113, 0.8890451, 1.351133, 0.8920303, 1.351134, 
        0.8924474, 1.351376, 0.9370828, 1.351376, 0.9373027, 1.351380, 0.9411675, 
        1.351417, 0.9525173, 1.351433, 0.9550566, 1.351576, 1.003863, 1.351578, 
        1.007517, 1.351579, 1.008607, 1.351584, 1.013875, 1.351588, 1.015127, 
        1.351599, 1.019788, 1.351690, 1.071442, 1.351691, 1.072240, 1.351692, 
        1.076411, 1.351695, 1.079858, 1.351700, 1.086642, 1.351756, 1.141328, 
        1.351756, 1.143567, 1.351758, 1.146712, 1.351763, 1.155619, 1.351793, 
        1.213414, 1.351793, 1.215688, 1.351796, 1.226717, 1.351813, 1.288085, 
        1.351813, 1.288809, 1.351813, 1.290776, 1.351815, 1.299938, 1.351824, 
        1.365075, 1.351824, 1.365629, 1.351825, 1.366590, 1.351825, 1.367936, 
        1.351825, 1.368763, 1.351825, 1.369224, 1.351825, 1.372589, 1.351826, 
        1.375282, 1.351831, 1.444142, 1.351831, 1.444339, 1.351831, 1.446359, 
        1.351831, 1.446615, 1.351831, 1.447393, 1.351831, 1.449412, 1.351831, 
        1.452747, 1.351834, 1.525773, 1.351834, 1.525903, 1.351834, 1.526512, 
        1.351834, 1.526936, 1.351834, 1.530821, 1.351834, 1.532335, 1.351836, 
        1.609042, 1.351836, 1.609224, 1.351836, 1.609897, 1.351836, 1.614046, 
        1.351836, 1.695459, 1.351836, 1.695656, 1.351836, 1.696248, 1.351836, 
        1.697878, 1.351837, 1.783325, 1.351837, 1.783463, 1.351837, 1.783834, 
        1.351837, 1.873237, 1.351837, 1.873983, 1.351837, 1.875665, 1.351837, 
        5.000000, 1.351837 }};
      auto totalSCR = processCoherentElastic(bragg,dwpix,dwp1,numSecondaryScatterers,secondaryScatterType,numEdges,tol,temps);
      REQUIRE( ranges::equal( totalSCR[0], correctSCR[0], equal ) );

    } // THEN
  } // GIVEN
} // TEST CASE
*/







TEST_CASE( "endout" ){
  REQUIRE( true );
  std::vector<double> alphas, betas, sab, temps, secondaryScatterVecThing (10,0.0), dwpix, dwp1;
  alphas = {1.1, 2.2, 3.3, 4.5, 5.8};
  betas = {0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7};
  sab = { 0.8812416231, 0.8226588624, 0.0371931173, 0.0370414517, 0.0453434977, 
          0.0449511778, 0.0177215328, 0.4944382555, 0.4922000032, 0.0833173024, 
          0.0740790498, 0.0698111325, 0.0693835223, 0.0339337297, 0.3223316052, 
          0.3256420103, 0.1180584841, 0.1055070878, 0.0815560280, 0.0810404205, 
          0.0458832427, 0.2182146674, 0.2226695868, 0.1282393202, 0.1182525679, 
          0.0860463732, 0.0851264533, 0.0537751752, 0.1498324933, 0.1534650775, 
          0.1195514020, 0.1135384304, 0.0852010403, 0.0844404310, 0.0572611451};
  temps = { 296.0 };
  dwpix = { 0.27366867553080776 };
  dwp1  = { 0.0 };





  std::cout.precision(15);
//std::cout << sab[0+4*betas.size()] << std::endl;
  double awr = 0.99917, spr = 20.449, aws = 15.85316, sps = 3.8883;

  int numSecondaryScatterers = 1, secondaryScatterType = 1;

  //endout(sab,awr,aws,spr,sps,temps,numSecondaryScatterers,secondaryScatterType,secondaryScatterVecThing,alphas,betas,dwpix,dwp1);

} // TEST CASE
