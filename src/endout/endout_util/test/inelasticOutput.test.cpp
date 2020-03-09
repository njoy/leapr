#include "catch.hpp"
#include "endout/endout_util/inelasticOutput.h"
#include "generalTools/testing.h"


TEST_CASE( "writing" ){
  using namespace njoy::ENDFtk;
  using ScatteringLawConstants = section::Type< 7, 4 >::ScatteringLawConstants;

  double za  = 127.0,
         awr = 8.934780e+0;
  int    lasym = 0,
           lat = 1;
  ScatteringLawConstants 
    constants( 0, 1.976285e+2, 5.000001e+0, 6.153875e+0, 8.934780e+0, 1 );

  std::vector<double> alphas { 0.1, 0.2, 0.3 }, 
                      betas  { 0.0, 0.2, 0.4, 0.6 },
  temps { 296.0, 400.0, 1200.0 },
  sab_temp_1 {3.803356e-2, 1.186118e-2, 5.35523e-3,  5.494297e-3,   // a0b0 a0b1 a0b2 a0b3
              7.283326e-2, 2.289232e-2, 1.153210e-2, 1.067055e-2,   // a1b0 a1b1 a1b2 a1b3
              1.046095e-1, 3.313806e-2, 1.680395e-2, 1.554271e-2 }, // a2b0 a2b1 a2b2 a2b3
  sab_temp_2 {6.888776e-2, 2.157749e-2, 1.085073e-2, 1.007578e-2,    
              1.308459e-1, 4.144943e-2, 2.104045e-2, 1.952161e-2,
              1.864122e-1, 5.972006e-2, 3.059757e-2, 2.836719e-2 },
  sab_temp_3 {5.744355e-1, 1.848110e-1, 9.506814e-2, 8.835550e-2, 
              1.011715e+0, 3.376009e-1, 1.787335e-1, 1.656053e-1,
              1.337622e+0, 4.627526e-1, 2.518537e-1, 2.328031e-1 };


  std::vector<std::vector<double>> sabFull { sab_temp_1, sab_temp_2, sab_temp_3 };




  writeToENDF(/*za, awr, lasym, lat, constants,*/ sabFull, alphas, betas, temps);
} // TEST CASE




TEST_CASE( "processing inelastic scattering data" ){
  GIVEN( "" ){  
    std::cout.precision(15);
    int lat = 1;

    std::vector<double>
      alphas {0.1, 0.2, 0.3},
      betas  {0.0, 0.2, 0.4, 0.6},
      temps  {296.0, 400.0, 1200.0};
    std::vector<std::vector<double>> fullSAB {
     { 3.80335716E-2, 1.30979776E-2, 7.237890430E-3, 7.39846449E-3, 7.283325728E-2, 
       2.52793721E-2, 1.40624592E-2, 1.436865581E-2, 0.104609458,   3.659346961E-2, 
       2.04910576E-2, 2.092936745E-2 },
     { 6.88877596E-2, 2.32208158E-2, 1.256643491E-2, 1.255765314E-2, 0.1308459204, 
       4.460619358E-2, 2.436734637E-2, 2.433019342E-2, 0.186412246, 6.426830031E-2, 
       3.543562590E-2, 3.535462320E-2 },
     { 0.5744354712, 0.1893883947,  9.98357465889E-2, 9.508460039E-2, 1.011714767, 
       0.3459625751, 0.1876968818, 0.178217691,  1.3376224712, 0.474214044, 
       0.264483953, 0.2505332640 } };
    std::vector<std::vector<double>> correctOut_beta0 (betas.size());
    std::vector<std::vector<double>> correctOut_beta1 (betas.size());
    
    WHEN( "isym = 0. No cold option used and symmetric S(a,b) requested" ){
      int isym = 0;
      AND_WHEN( "No log option is used" ){
        int ilog = 0;

        auto out0 = getSABreadyToWrite(fullSAB,temps,alphas,betas,isym,ilog,lat,0);
        auto out1 = getSABreadyToWrite(fullSAB,temps,alphas,betas,isym,ilog,lat,1);

        correctOut_beta0 = { {0.03803357, 0.07283326, 0.1046095}, 
                             {0.06888776, 0.1308459,  0.1864122}, 
                             {0.5744355,  1.011715,   1.337622 } };
        correctOut_beta1 = { {0.01186118, 0.02289232, 0.03313806}, 
                             {0.02157749, 0.04144943, 0.05972006}, 
                             {0.1848110,  0.3376009,  0.4627526}};
        auto toWrite0 = std::get<1>(out0),
             toWrite1 = std::get<1>(out1);

        REQUIRE( 0.0 == Approx(std::get<0>(out0)).epsilon(1e-6) );
        REQUIRE( 0.2 == Approx(std::get<0>(out1)).epsilon(1e-6) );

        for ( size_t a = 0; a < alphas.size(); ++a ){
          REQUIRE( ranges::equal(toWrite0[a],correctOut_beta0[a],equal) );
          REQUIRE( ranges::equal(toWrite1[a],correctOut_beta1[a],equal) );
        }
      } // AND WHEN
      AND_WHEN( "Log option is used - log10(S(a,b)) returned" ){
        int ilog = 1;
        auto out0 = getSABreadyToWrite(fullSAB,temps,alphas,betas,isym,ilog,lat,0);
        auto out1 = getSABreadyToWrite(fullSAB,temps,alphas,betas,isym,ilog,lat,1);
        auto toWrite0 = std::get<1>(out0);
        auto toWrite1 = std::get<1>(out1);
    
        correctOut_beta0 = {{-3.269286,  -2.619583, -2.257521},
                            {-2.675277,  -2.033735, -1.679795}, 
                            {-0.5543675, 0.01164668, 0.2908938}};

        correctOut_beta1 = {{-4.434485, -3.776954, -3.407073 }, 
                            {-3.836105, -3.183281, -2.818087 }, 
                            {-1.688422, -1.085891, -0.7705627 }};

        REQUIRE( 0.0 == Approx(std::get<0>(out0)).epsilon(1e-6) );
        REQUIRE( 0.2 == Approx(std::get<0>(out1)).epsilon(1e-6) );
        for ( size_t a = 0; a < alphas.size(); ++a ){
          REQUIRE( ranges::equal(toWrite0[a],correctOut_beta0[a],equal) );
          REQUIRE( ranges::equal(toWrite1[a],correctOut_beta1[a],equal) );
        }
      } // AND WHEN
    } // AND WHEN


    WHEN( "isym = 1. Ncold option is used and symmetric S(a,b) requested" ){
      int isym = 1;

      std::vector<std::vector<double>> fullSAB_1 {{0.04952457206, 0.02137508058, 
      0.01418638218, 0.01066516070, 0.09102180354, 0.03892860779, 0.02592329276, 0.02082273114, 
      0.1247384714, 0.05227166194, 0.03452343287, 0.02930870218}, {0.09008423197, 0.03788786274, 
      0.02508671669, 0.01814987511, 0.1644316393, 0.06838177727, 0.04497266379, 0.03526222649, 
      0.2235133181, 0.09144674298, 0.05954545052, 0.04951281781}, {0.7234388542, 0.2700462858, 
      0.1634772382, 0.1318445112, 1.219476539, 0.4629902995, 0.2819640002, 0.2432061438, 
      1.537234344, 0.5986513806, 0.3704009522, 0.3293786514}};

      std::vector<std::vector<double>> fullSAB_2{{0.02664995522, 6.70965907E-3, 
        3.90001531E-3, 4.03631178E-3, 0.05008646793, 0.01298147618, 7.66290745E-3, 
        8.00561935E-3, 0.07026869416, 0.01879841914, 0.01131220822, 0.01182831703},
        {0.04833979060, 0.01258480466, 7.54491499E-3, 8.04413321E-3, 0.09022655754, 
        0.02454353198, 0.01507582038, 0.01608328746, 0.12546381684, 0.03554293119, 
        0.02236444746, 0.02370712345}, {0.38751953235, 0.11412210887, 0.07529282489, 
        0.08297012476, 0.66382656949, 0.21325391570, 0.14676749529, 0.15729955905, 
        0.84980927694, 0.29180756145, 0.20664996179, 0.21832926074}};
      AND_WHEN( "No log option is used" ){
        int ilog = 0;

        std::vector<std::vector<std::vector<double>>> correctOut {
        {{7.920232E-3, 0.01546351, 0.02176542}, {0.01456276, 0.02829306, 0.03972719},
         {0.1225139, 0.2259945, 0.3060686}},
        {{0.01163372, 0.02125872, 0.02831138},{0.02166160, 0.03883250, 0.05141565},
         {0.1556705, 0.2684989, 0.3527126}},
        {{0.01935670, 0.03525270, 0.04733581},{0.03520655, 0.06354242, 0.08497509},
         {0.2635194, 0.4518001, 0.5841824}},
        {{0.02664996, 0.05008647, 0.07026869},{0.04833979, 0.09022656, 0.1254638},
         {0.3875195, 0.6638266, 0.8498093}},
        {{7.409296E-3, 0.01433510, 0.02075859},{0.01354326, 0.02641275, 0.03824986},
         {0.1169487, 0.2185358, 0.2990350}},
        {{4.755753E-3, 9.344296E-3, 0.01379432},{8.737910E-3, 0.01745960, 0.02590069},
         {0.07906871, 0.1541278, 0.2170133}},
        {{5.435183E-3, 0.01078014, 0.01592768}, {0.01002557, 0.02004493, 0.02954668},
         {0.08928908, 0.1692794, 0.2349571}}};

        std::vector<double> betaVals {-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6};


        for ( int b = 0; b < 7; ++b ){
          auto out0 = getSABreadyToWrite(fullSAB_1,temps,alphas,betas,isym,ilog,lat,b,fullSAB_2);
          auto toWrite0 = std::get<1>(out0);
          REQUIRE( betaVals[b] == Approx(std::get<0>(out0)).epsilon(1e-6) );
          for ( size_t a = 0; a < alphas.size(); ++a ){
            REQUIRE( ranges::equal(toWrite0[a],correctOut[b][a],equal) );
          }
        }

      } // AND WHEN

      AND_WHEN( "Log option is used - log10(S(a,b)) returned" ){
        int ilog = 1;

        std::vector<std::vector<std::vector<double>>> correctOut 
        {{{-4.838335, -4.169272, -3.827433}, {-4.229287, -3.565139, -3.225720}, 
          {-2.099531, -1.487244, -1.183946}},
         {{-4.453847, -3.850988, -3.564492}, {-3.832214, -3.248498, -2.967813}, 
          {-1.860014, -1.314908, -1.042102}},
         {{-3.944717, -3.345213, -3.050488}, {-3.346523, -2.756048, -2.465397}, 
          {-1.333628, -0.7945154, -0.5375421}},
         {{-3.624968, -2.994004, -2.655429}, {-3.029500, -2.405431, -2.075738}, 
          {-0.9479890, -0.4097344, -0.1627433}},
         {{-4.905020, -4.245045, -3.874795}, {-4.301867, -3.633908, -3.263615}, 
          {-2.146020, -1.520806, -1.207195}},
         {{-5.348400, -4.672989, -4.283498}, {-4.740084, -4.047866, -3.653486}, 
          {-2.537438, -1.869973, -1.527797}},
         {{-5.214862, -4.530050, -4.139697}, {-4.602616, -3.909779, -3.521784}, 
          {-2.415876, -1.776205, -1.448352}}};

        std::vector<double> betaVals {-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6};


        for ( int b = 0; b < 7; ++b ){
          auto out0 = getSABreadyToWrite(fullSAB_1,temps,alphas,betas,isym,ilog,
                                         lat,b,fullSAB_2);
          auto toWrite0 = std::get<1>(out0);
          REQUIRE( betaVals[b] == Approx(std::get<0>(out0)).epsilon(1e-6) );
          for ( size_t a = 0; a < alphas.size(); ++a ){
            REQUIRE( ranges::equal(toWrite0[a],correctOut[b][a],equal) );
          }
        }
      } // AND WHEN

    } // WHEN







    WHEN( "isym = 2. No cold option used and non-symmetric S(a,b) requested" ){
      int isym = 2;

      AND_WHEN( "No log option is used" ){
        int ilog = 0;
        auto out0 = getSABreadyToWrite(fullSAB,temps,alphas,betas,isym,ilog,lat,0);
        auto out1 = getSABreadyToWrite(fullSAB,temps,alphas,betas,isym,ilog,lat,1);
        auto toWrite0 = std::get<1>(out0),
             toWrite1 = std::get<1>(out1);
  
        correctOut_beta0 =
        {{0.03803357, 0.07283326, 0.1046095}, 
         {0.06888776, 0.1308459,  0.1864122}, 
         {0.5744355,  1.011715,   1.337622}};
        correctOut_beta1 =
        {{0.01309798, 0.02527937, 0.03659347},
         {0.02322082, 0.04460619, 0.06426830},
         {0.1893884, 0.3459626, 0.4742140}};

        REQUIRE( 0.0 == Approx(std::get<0>(out0)).epsilon(1e-6) );
        REQUIRE( 0.2 == Approx(std::get<0>(out1)).epsilon(1e-6) );
        for ( size_t a = 0; a < alphas.size(); ++a ){
          REQUIRE( ranges::equal(toWrite0[a],correctOut_beta0[a],equal) );
          REQUIRE( ranges::equal(toWrite1[a],correctOut_beta1[a],equal) );
        }
      } // AND WHEN

      AND_WHEN( "Log option is used - log10(S(a,b)) returned" ){
        int ilog = 1;
        auto out0 = getSABreadyToWrite(fullSAB,temps,alphas,betas,isym,ilog,lat,0);
        auto out1 = getSABreadyToWrite(fullSAB,temps,alphas,betas,isym,ilog,lat,1);
        auto toWrite0 = std::get<1>(out0);
        auto toWrite1 = std::get<1>(out1);
    
        correctOut_beta0 = {{-3.269286, -2.619583, -2.257521},
                            {-2.675277, -2.033735, -1.679795}, 
                            {-0.5543675, 0.01164668, 0.2908938}};
        correctOut_beta1 = {{-4.335297, -3.677767, -3.307885},
                            {-3.762706, -3.109883, -2.744689},
                            {-1.663955, -1.061425, -0.7460965}};

        REQUIRE( 0.0 == Approx(std::get<0>(out0)).epsilon(1e-6) );
        REQUIRE( 0.2 == Approx(std::get<0>(out1)).epsilon(1e-6) );
        for ( size_t a = 0; a < alphas.size(); ++a ){
          REQUIRE( ranges::equal(toWrite0[a],correctOut_beta0[a],equal) );
          REQUIRE( ranges::equal(toWrite1[a],correctOut_beta1[a],equal) );
        }

      } // AND WHEN

    } // WHEN




    WHEN( "isym = 3. Ncold option is used and non-symmetric S(a,b) requested" ){
      int isym = 3;

      std::vector<std::vector<double>> fullSAB_1 {{0.04952457206, 0.02137508058, 
      0.01418638218, 0.01066516070, 0.09102180354, 0.03892860779, 0.02592329276, 
      0.02082273114, 0.1247384714, 0.05227166194, 0.03452343287, 0.02930870218}, 
     {0.09008423197, 0.03788786274, 0.02508671669, 0.01814987511, 0.1644316393, 
      0.06838177727, 0.04497266379, 0.03526222649, 0.2235133181, 0.09144674298, 
      0.05954545052, 0.04951281781}, {0.7234388542, 0.2700462858, 0.1634772382, 
      0.1318445112, 1.219476539, 0.4629902995, 0.2819640002, 0.2432061438, 
      1.537234344, 0.5986513806, 0.3704009522, 0.3293786514}};

      std::vector<std::vector<double>> fullSAB_2{{0.02664995522, 6.70965907E-3, 
        3.90001531E-3, 4.03631178E-3, 0.05008646793, 0.01298147618, 7.66290745E-3, 
        8.00561935E-3, 0.07026869416, 0.01879841914, 0.01131220822, 0.01182831703},
        {0.04833979060, 0.01258480466, 7.54491499E-3, 8.04413321E-3, 0.09022655754, 
        0.02454353198, 0.01507582038, 0.01608328746, 0.12546381684, 0.03554293119, 
        0.02236444746, 0.02370712345}, {0.38751953235, 0.11412210887, 0.07529282489, 
        0.08297012476, 0.66382656949, 0.21325391570, 0.14676749529, 0.15729955905, 
        0.84980927694, 0.29180756145, 0.20664996179, 0.21832926074}};
      AND_WHEN( "No log option is used" ){
        int ilog = 0;

        std::vector<std::vector<std::vector<double>>> correctOut 
            {{{1.066516E-2, 2.082273E-2, 2.930870E-2},
{1.814988E-2, 3.526223E-2, 4.951282E-2},
{1.318445E-1, 2.432061E-1, 3.293787E-1}},
{{1.418638E-2, 2.592329E-2, 3.452343E-2},
{2.508672E-2, 4.497266E-2, 5.954545E-2},
{1.634772E-1, 2.819640E-1, 3.704010E-1}},
{{2.137508E-2, 3.892861E-2, 5.227166E-2},
{3.788786E-2, 6.838178E-2, 9.144674E-2},
{2.700463E-1, 4.629903E-1, 5.986514E-1}},
{{2.664996E-2, 5.008647E-2, 7.026869E-2},
{4.833979E-2, 9.022656E-2, 1.254638E-1},
{3.875195E-1, 6.638266E-1, 8.498093E-1}},
{{6.709659E-3, 1.298148E-2, 1.879842E-2},
{1.258480E-2, 2.454353E-2, 3.554293E-2},
{1.141221E-1, 2.132539E-1, 2.918076E-1}},
{{3.900015E-3, 7.662907E-3, 1.131221E-2},
{7.544915E-3, 1.507582E-2, 2.236445E-2},
{7.529282E-2, 1.467675E-1, 2.066500E-1}},
{{4.036312E-3, 8.005619E-3, 1.182832E-2},
{8.044133E-3, 1.608329E-2, 2.370712E-2},
{8.297012E-2, 1.572996E-1, 2.183293E-1}}};

        std::vector<double> betaVals {-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6};


        for ( int b = 0; b < 7; ++b ){
          auto out0 = getSABreadyToWrite(fullSAB_1,temps,alphas,betas,isym,ilog,lat,b,fullSAB_2);
          auto toWrite0 = std::get<1>(out0);
          REQUIRE( betaVals[b] == Approx(std::get<0>(out0)).epsilon(1e-6) );
          for ( size_t a = 0; a < alphas.size(); ++a ){
            REQUIRE( ranges::equal(toWrite0[a],correctOut[b][a],equal) );
          }
        }

      } // AND WHEN

      AND_WHEN( "Log option is used - log10(S(a,b)) returned" ){
        int ilog = 1;

        std::vector<std::vector<std::vector<double>>> correctOut {
        {{-4.540773, -3.871710, -3.529871},
         {-4.009092, -3.344943, -3.005524},
         {-2.026132, -1.413846, -1.110547}},
        {{-4.255473, -3.652613, -3.366117},
         {-3.685417, -3.101700, -2.821015},
         {-1.811082, -1.265976, -0.9931692}},
        {{-3.845529, -3.246026, -2.951301},
         {-3.273124, -2.682649, -2.391999},
         {-1.309162, -0.7700492, -0.5130759}},
        {{-3.624968, -2.994004, -2.655429},
         {-3.029500, -2.405431, -2.075738},
         {-0.9479890, -0.4097344, -0.1627433}},
        {{-5.004207, -4.344232, -3.973983},
         {-4.375265, -3.707307, -3.337014},
         {-2.170486, -1.545272, -1.231661}},
        {{-5.546775, -4.871364, -4.481873},
         {-4.886881, -4.194663, -3.800283},
         {-2.586370, -1.918906, -1.576729}},
        {{-5.512424, -4.827612, -4.437259},
         {-4.822812, -4.129975, -3.741980},
         {-2.489275, -1.849603, -1.521751}}};
        


        std::vector<double> betaVals {-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6};

        for ( int b = 0; b < 7; ++b ){
          auto out0 = getSABreadyToWrite(fullSAB_1,temps,alphas,betas,isym,ilog,
                                         lat,b,fullSAB_2);
          auto toWrite0 = std::get<1>(out0);
          REQUIRE( betaVals[b] == Approx(std::get<0>(out0)).epsilon(1e-6) );
          for ( size_t a = 0; a < alphas.size(); ++a ){
            REQUIRE( ranges::equal(toWrite0[a],correctOut[b][a],equal) );
          }
        }
      } // AND WHEN
    } // WHEN
  } // GIVEN
} // TEST CASE


TEST_CASE( "processing inelastic scattering data ( all betas )" ){
  GIVEN( "" ){  
    std::cout.precision(15);
    int lat = 1;

    std::vector<double>
      alphas {0.1, 0.2, 0.3},
      betas  {0.0, 0.2},//, 0.4, 0.6},
      temps  {296.0, 400.0, 1200.0};
    std::vector<std::vector<double>> fullSAB {
     { 3.80335716E-2, 1.30979776E-2, 7.237890430E-3, 7.39846449E-3, 7.283325728E-2, 
       2.52793721E-2, 1.40624592E-2, 1.436865581E-2, 0.104609458,   3.659346961E-2, 
       2.04910576E-2, 2.092936745E-2 },
     { 6.88877596E-2, 2.32208158E-2, 1.256643491E-2, 1.255765314E-2, 0.1308459204, 
       4.460619358E-2, 2.436734637E-2, 2.433019342E-2, 0.186412246, 6.426830031E-2, 
       3.543562590E-2, 3.535462320E-2 },
     { 0.5744354712, 0.1893883947,  9.98357465889E-2, 9.508460039E-2, 1.011714767, 
       0.3459625751, 0.1876968818, 0.178217691,  1.3376224712, 0.474214044, 
       0.264483953, 0.2505332640 } };
    std::vector<std::vector<double>> correctOut_beta0 (betas.size());
    std::vector<std::vector<double>> correctOut_beta1 (betas.size());
    
    WHEN( "isym = 0. No cold option used and symmetric S(a,b) requested" ){
      int isym = 0;
      AND_WHEN( "No log option is used" ){
        int ilog = 0;
        inelasticOutput(alphas,betas,fullSAB,temps,isym,ilog,lat);
      } // AND WHEN
    } // WHEN
  } // GIVEN
} // TEST CASE




/*


*/







