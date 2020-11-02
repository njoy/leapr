#define CATCH_CONFIG_MAIN
#include <iostream>
#include "catch.hpp"
#include "coherentElastic/coherentElastic.h"
#include "generalTools/testing.h"
#include <range/v3/all.hpp>


TEST_CASE( "coherentElastic" ){


  int iel, npr = 1;
  std::vector<double> b ( 60000, 0.0 );
  double emax = 5.0;
  std::vector<double> braggVals_0_to_99 (100);
  std::vector<double> braggVals_500_to_599 (100);
  std::vector<double> braggVals_600_to_699 (100);
  std::vector<double> braggVals_2000_to_2099 (100);
  GIVEN( "graphite is the requested material" ){
    iel = 1;
    WHEN( "1 principal scattering atoms in compound" ){
      npr = 1;
      auto output = coherentElastic( iel, npr, b, emax );
      int nw   = std::get<0>(output);
      int maxb = std::get<1>(output);

      braggVals_0_to_99 = {4.555814E-4, 0.000000000, 1.822325E-3, 1.369340E-2, 
        4.100233E-3, 0.000000000, 4.515834E-3, 1.631010E-3, 4.971416E-3, 
        9.326892E-3, 6.338160E-3, 2.753430E-3, 7.289303E-3, 6.846702E-3, 
        8.616068E-3, 7.084716E-3, 1.138953E-2, 0.000000000, 1.180513E-2, 
        2.017531E-3, 1.354750E-2, 1.506663E-2, 1.400308E-2, 2.12213E-28, 
        1.536983E-2, 2.829053E-2, 1.590537E-2, 5.214409E-3, 1.640093E-2, 
        4.564468E-3, 1.764773E-2, 1.89034E-28, 1.806333E-2, 8.155052E-4, 
        1.851892E-2, 4.832470E-3, 1.988566E-2, 1.554482E-3, 2.083680E-2, 
        2.429742E-2, 2.091676E-2, 1.515683E-3, 2.216357E-2, 4.417302E-3, 
        2.232349E-2, 0.000000000, 2.493704E-2, 1.59024E-28, 2.535264E-2, 
        1.376715E-3, 2.683932E-2, 4.014128E-3, 2.915721E-2, 3.423351E-3, 
        2.945287E-2, 3.831891E-3, 2.994843E-2, 2.026696E-2, 3.161084E-2, 
        1.232928E-3, 3.206642E-2, 7.344830E-3, 3.343316E-2, 2.397711E-3, 
        3.367304E-2, 1.194578E-3, 3.446427E-2, 1.180786E-3, 3.571107E-2, 
        6.959940E-3, 3.587099E-2, 1.32590E-28, 3.690209E-2, 0.000000000, 
        3.890014E-2, 2.222850E-3, 4.038683E-2, 3.272330E-3, 4.064251E-2, 
        8.698722E-3, 4.109809E-2, 3.71617E-28, 4.141793E-2, 3.231341E-3, 
        4.246483E-2, 1.702005E-2, 4.270471E-2, 1.697218E-2, 4.300038E-2, 
        6.342651E-3, 4.474274E-2, 3.56160E-28, 4.555814E-2, 2.738681E-3, 
        4.722055E-2, 1.008765E-3, 4.793181E-2, 1.602004E-2, 4.801177E-2, 
        2.000837E-3 };
      braggVals_500_to_599={1.899752263, 3.180806E-4, 1.900156839, 5.406474E-3, 
        1.900944241, 8.425241E-3, 1.902129752, 2.557797E-2, 1.906704233, 
        9.898097E-2, 1.924831731, 5.173886E-1, 2.021217751, 1.541571E-2, 
        2.022442286, 2.927698E-3, 2.024057289, 3.851379E-3, 2.025155416, 
        1.272094E-1, 2.045105240, 5.731394E-1, 2.147446376, 2.973445E-2, 
        2.155016495, 2.389188E-3, 2.155075141, 1.253819E-2, 2.158001108, 
        6.499047E-2, 2.169023401, 5.746446E-1, 2.277761344, 1.42420E-27, 
        2.277861663, 1.655197E-2, 2.281316703, 2.902396E-3, 2.281919269, 
        7.361348E-2, 2.296586214, 6.828792E-1, 2.411696567, 5.04745E-27, 
        2.412260108, 8.468275E-4, 2.412276864, 5.363164E-3, 2.412477720, 
        6.440965E-2, 2.427793679, 6.608077E-1, 2.549570411, 2.195730E-2, 
        2.552768814, 2.030029E-2, 2.556981565, 2.193049E-3, 2.558192209, 
        2.999576E-2, 2.562645795, 7.141287E-1, 2.690785822, 1.603541E-3, 
        2.691390298, 5.610789E-3, 2.693306475, 7.879202E-3, 2.695385437, 
        3.188995E-2, 2.701142564, 7.190485E-1, 2.836369129, 4.165089E-3, 
        2.836399776, 1.301366E-2, 2.837652310, 2.627506E-2, 2.843002213, 
        4.810053E-3, 2.843283984, 7.930234E-1, 2.985463256, 4.820958E-3, 
        2.985640150, 2.283299E-3, 2.986861602, 3.170728E-3, 2.987974793, 
        9.256786E-3, 2.989070056, 8.526749E-1, 3.139001636, 7.423543E-4, 
        3.139178530, 6.680131E-3, 3.140491186, 2.473927E-4, 3.140566588, 
        6.183789E-3 };
      THEN( "bragg edges vector is correctly output" ){
        checkPartOfVec(b,braggVals_500_to_599,500);
        restAreZero(nw,b);
        REQUIRE( nw   == 852 );
        REQUIRE( maxb == 690 );
      } // THEN
    } // WHEN
    WHEN( "3 principal scattering atoms in compound" ){
      npr = 3;
      auto output = coherentElastic( iel, npr, b, emax );
      int nw   = std::get<0>(output);
      int maxb = std::get<1>(output);

      braggVals_0_to_99 = { 4.555814E-4, 0.000000000, 1.822325E-3, 4.564468E-3, 
        4.100233E-3, 0.000000000, 4.515834E-3, 5.436701E-4, 4.971416E-3, 
        3.108964E-3, 6.338160E-3, 9.178102E-4, 7.289303E-3, 2.282234E-3, 
        8.616068E-3, 2.361572E-3, 1.138953E-2, 0.000000000, 1.180513E-2, 
        6.725103E-4, 1.354750E-2, 5.022209E-3, 1.400308E-2, 7.07379E-29, 
        1.536983E-2, 9.430179E-3, 1.590537E-2, 1.738136E-3, 1.640093E-2, 
        1.521489E-3, 1.764773E-2, 6.30115E-29, 1.806333E-2, 2.718350E-4, 
        1.851892E-2, 1.610823E-3, 1.988566E-2, 5.181606E-4, 2.083680E-2, 
        8.099140E-3, 2.091676E-2, 5.052277E-4, 2.216357E-2, 1.472434E-3, 
        2.232349E-2, 0.000000000, 2.493704E-2, 5.30080E-29, 2.535264E-2, 
        4.589051E-4, 2.683932E-2, 1.338042E-3, 2.915721E-2, 1.141117E-3, 
        2.945287E-2, 1.277297E-3, 2.994843E-2, 6.755654E-3, 3.161084E-2, 
        4.109760E-4, 3.206642E-2, 2.448276E-3, 3.343316E-2, 7.992373E-4, 
        3.367304E-2, 3.981927E-4, 3.446427E-2, 3.935953E-4, 3.571107E-2, 
        2.319980E-3, 3.587099E-2, 4.41969E-29, 3.690209E-2, 0.000000000, 
        3.890014E-2, 7.409499E-4, 4.038683E-2, 1.090776E-3, 4.064251E-2, 
        2.899574E-3, 4.109809E-2, 1.23872E-28, 4.141793E-2, 1.077113E-3, 
        4.246483E-2, 5.673352E-3, 4.270471E-2, 5.657395E-3, 4.300038E-2, 
        2.114217E-3, 4.474274E-2, 1.18720E-28, 4.555814E-2, 9.128937E-4, 
        4.722055E-2, 3.362551E-4, 4.793181E-2, 5.340016E-3, 4.801177E-2, 
        6.669459E-4 };
      THEN( "bragg edges vector is correctly output" ){
        checkVec(b,braggVals_0_to_99,braggVals_0_to_99.size());
        restAreZero(nw,b);
        REQUIRE( nw   == 852 );
        REQUIRE( maxb == 690 );
      } // THEN
    } // WHEN
  } // GIVEN





  GIVEN( "beryllium metal is the requested material" ){
    iel = 2;
    WHEN( "1 principal scattering atom in compound" ){
      npr = 1;
      auto output = coherentElastic( iel, npr, b, emax );
      int nw   = std::get<0>(output);
      int maxb = std::get<1>(output);

      braggVals_0_to_99 = { 1.5928452E-3, 0.0000000000, 5.2198010E-3, 
        8.9780237E-3, 6.3713807E-3, 1.0835019E-2, 6.8126461E-3, 4.7152102E-2, 
        1.1591182E-2, 1.2049634E-2, 1.4335607E-2, 0.0000000000, 1.5659403E-2, 
        2.0733858E-2, 1.7252248E-2, 0.0000000000, 1.9555407E-2, 2.7830790E-2, 
        2.0879204E-2, 4.4890118E-3, 2.2030784E-2, 3.4960910E-2, 2.2472049E-2, 
        2.5961969E-2, 2.5485523E-2, 5.4175095E-3, 2.7250585E-2, 7.8586836E-3, 
        2.9995009E-2, 0.0000000000, 3.0705324E-2, 7.4033942E-3, 3.5214810E-2, 
        2.0739405E-2, 3.6538607E-2, 6.7867480E-3, 3.8131452E-2, 3.9860918E-2, 
        3.9821129E-2, 0.0000000000, 4.1144926E-2, 2.5582299E-2, 4.2909987E-2, 
        1.2525311E-2, 4.5040930E-2, 1.8338140E-2, 4.6364727E-2, 6.0248170E-3, 
        4.6978209E-2, 1.1970698E-2, 4.8571054E-2, 3.4509615E-2, 5.3349589E-2, 
        2.2466332E-2, 5.7342426E-2, 3.6116730E-3, 6.0700333E-2, 2.6214666E-2, 
        6.2562227E-2, 1.5553515E-2, 6.7857413E-2, 5.4269672E-2, 7.2463731E-2, 
        4.8005807E-2, 7.6359736E-2, 2.8168048E-2, 7.8049413E-2, 4.6384650E-3, 
        8.2193019E-2, 2.7150104E-2, 8.3269214E-2, 2.9071970E-2, 8.8123134E-2, 
        1.7480455E-2, 8.9888196E-2, 2.1287287E-2, 9.7852422E-2, 1.2441512E-2, 
        9.8928617E-2, 4.1013331E-2, 1.0194209E-1, 2.6761195E-2, 1.0716189E-1, 
        4.7286135E-2, 1.1351182E-1, 2.3103002E-2, 1.1458802E-1, 5.3467779E-2, 
        1.1760149E-1, 3.0112907E-2, 1.2282129E-1, 2.9464778E-2, 1.2902046E-1, 
        5.1362008E-2, 1.3686640E-1, 3.1356724E-2, 1.4085924E-1, 7.8135812E-2, 
        1.4892030E-1, 1.3446856E-2 };

      braggVals_500_to_599 = { 4.7065676E00, 1.6604899E-1, 4.7282533E00, 
        6.4950631E-1, 4.8183566E00, 1.3988329E00, 5.0000000000, 1.3988329E00, 
        5.1687479E19, 1.0430078E-6, 5.4286323E19, 2.1980656E-8, 5.4340330E19, 
        1.0335175E-6, 5.6853152E19, 1.1973872E-6, 5.9708819E19, 4.6588986E-9, 
        5.9715242E19, 3.8813193E-8, 5.9780690E19, 1.0862852E-8, 5.9824736E19, 
        9.3087658E-9, 5.9828462E19, 1.5514126E-9, 5.9830520E19, 6.2052733E-9, 
        5.9849224E19, 3.0237884E-8, 5.9931282E19, 9.2269901E-7, 6.2264810E19, 
        1.3321733E-6, 6.5388441E19, 2.3737156E-8, 6.5465932E19, 9.5279738E-7, 
        6.7922452E19, 1.4115651E-6, 7.1327857E19, 7.1035314E-9, 7.1362922E19, 
        5.1111418E-8, 7.1511669E19, 1.9153661E-8, 7.1549894E19, 1.0561234E-7, 
        7.1750377E19, 8.8028590E-7, 7.3826078E19, 1.4825885E-6, 7.7532303E19, 
        4.0884732E-9, 7.7534373E19, 2.7255794E-9, 7.7550167E19, 0.0000000000, 
        7.7554283E19, 4.0878938E-9, 7.7557406E19, 6.8128358E-9, 7.7575471E19, 
        9.5290865E-8, 7.7776996E19, 9.1111482E-7, 7.9975689E19, 1.6443680E-6, 
        8.3998954E19, 3.9279465E-9, 8.3999368E19, 1.9637110E-8, 8.4049599E19, 
        9.3266209E-7, 8.6371284E19, 1.8045704E-6, 9.0699910E19, 6.3001005E-9, 
        9.0704063E19, 1.0079690E-8, 9.0716807E19, 1.5117120E-8, 9.0767702E19, 
        3.2110741E-8, 9.0854905E19, 1.7479449E-7, 9.1273523E19, 6.8371463E-7, 
        9.3012864E19, 1.4725069E-6, 9.6519281E19, 1.4725069E-6, 0.0000000000, 
        0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 0.0000000000, 
        0.0000000000, 0.0000000000 ,0.00};
      THEN( "bragg edges vector is correctly output" ){
        checkPartOfVec(b,braggVals_0_to_99);
        checkPartOfVec(b,braggVals_500_to_599,500);
        restAreZero(nw,b);
        REQUIRE( nw   == 592 );
        REQUIRE( maxb == 508 );
      } // THEN
    } // WHEN
  } // GIVEN





  GIVEN( "beryllium oxide is the requested material" ){
    iel = 3;
    WHEN( "1 principal scattering atom in compound" ){
      npr = 1;
      auto output = coherentElastic( iel, npr, b, emax );
      int nw   = std::get<0>(output);
      int maxb = std::get<1>(output);

      braggVals_0_to_99 = { 1.0611740E-3, 0.0000000000, 3.7543682E-3, 
        1.9057139E-2, 4.2446962E-3, 1.2191662E-2, 4.8155422E-3, 1.6539555E-2, 
        7.9990644E-3, 1.3321646E-2, 9.5505664E-3, 0.0000000000, 1.1263105E-2, 
        4.4010576E-2, 1.2324279E-2, 0.0000000000, 1.3304935E-2, 5.2025441E-2, 
        1.5017473E-2, 9.5285693E-3, 1.5507801E-2, 3.8270343E-2, 1.6078647E-2, 
        9.0515232E-3, 1.6978785E-2, 2.4321227E-4, 1.9262169E-2, 8.5846952E-3, 
        2.0733153E-2, 3.3013941E-4, 2.0813671E-2, 0.0000000000, 2.4568039E-2, 
        3.8285727E-2, 2.6280577E-2, 1.4405843E-2, 2.6529351E-2, 0.0000000000, 
        2.7341751E-2, 1.3882358E-2, 2.8241889E-2, 1.1314703E-3, 3.0283719E-2, 
        3.4483992E-2, 3.0525273E-2, 1.3638848E-2, 3.1996257E-2, 2.6575439E-4, 
        3.3789314E-2, 2.5409518E-2, 3.4850488E-2, 0.0000000000, 3.5831144E-2, 
        6.3404752E-2, 3.7792456E-2, 0.0000000000, 3.8034010E-2, 2.4437198E-2, 
        3.8202266E-2, 4.0638873E-3, 4.1546824E-2, 2.9441046E-2, 4.1956634E-2, 
        5.8167065E-3, 4.3259362E-2, 4.5710922E-4, 4.3339880E-2, 0.0000000000, 
        4.5052418E-2, 2.2005288E-2, 4.6113592E-2, 0.0000000000, 4.8806786E-2, 
        1.0570998E-2, 4.9297114E-2, 2.1464794E-2, 4.9465370E-2, 2.1428256E-2, 
        4.9867960E-2, 1.0279349E-2, 5.0768098E-2, 8.4390699E-4, 5.1997528E-2, 
        6.7737237E-2, 5.4602984E-2, 0.0000000000, 5.5751897E-2, 5.4543503E-2, 
        6.0069891E-2, 4.7642846E-3, 6.0318665E-2, 5.4055677E-3, 6.4314587E-2, 
        1.4452725E-2, 6.7015001E-2, 4.4336398E-3, 6.7915139E-2, 2.8717537E-2, 
        7.1332995E-2, 8.7440148E-3 };
      braggVals_500_to_599 = { 2.469195365, 1.041153E-1, 2.478734324, 
        8.075297E-1, 2.547878893, 1.518366455, 2.675464624, 8.454533E-2, 
        2.683036147, 4.668016E-2, 2.686724437, 8.637980E-1, 2.760113703, 
        1.664943516, 2.898372226, 1.371762E-3, 2.898925341, 9.848081E-3, 
        2.899165798, 1.116743E-4, 2.899223559, 3.926320E-2, 2.903203943, 
        8.294982E-1, 2.980837905, 1.849293301, 3.130023021, 1.074771E-4, 
        3.130181863, 1.761500E-2, 3.131994385, 5.802377E-2, 3.135550520, 
        1.877500E-1, 3.154453419, 6.694512E-1, 3.210051500, 1.897611779, 
        3.371090032, 1.307338E-2, 3.371489331, 2.740183E-2, 3.373658291, 
        1.364833E-1, 3.387911710, 6.695892E-1, 3.447754487, 2.018569562, 
        3.620272088, 2.359544E-2, 3.623718825, 4.394513E-2, 3.624595109, 
        5.794462E-2, 3.629859393, 8.094373E-1, 3.693946866, 2.150195585, 
        3.878815433, 2.110670E-2, 3.880296469, 8.034692E-1, 3.948628638, 
        2.364601146, 4.146201316, 1.183514E-2, 4.147705568, 1.178631E-2, 
        4.148225415, 3.264936E-1, 4.173012251, 4.097022E-1, 4.211799803, 
        2.490940271, 4.422841742, 0.000000000, 4.422894477, 5.706881E-3, 
        4.423154400, 1.908630E-2, 4.425410638, 2.129948E-2, 4.426004701, 
        5.455574E-4, 4.426466786, 4.428739E-2, 4.429803027, 1.014795E-1, 
        4.440428111, 5.105310E-1, 4.483460359, 2.580070534, 4.707784476, 
        0.000000000, 4.708044399, 8.763193E-5, 4.708247021, 2.841540E-2, 
        4.710216650, 4.380657E-5 };
      braggVals_600_to_699 = { 4.710394958, 3.916851E-2, 4.711798767, 
        3.214987E-2, 4.713641267, 3.276534E-2, 4.716333364, 5.387837E-1, 
        4.763610308, 2.801887250, 5.000000000, 2.801887250, 4.148162E19, 
        1.683046E-5, 4.355674E19, 1.288690E-7, 4.356315E19, 8.487957E-8, 
        4.358272E19, 4.352578E-7, 4.363562E19, 1.335506E-5, 4.525081E19, 
        1.910801E-5, 4.751348E19, 4.026319E-8, 4.751531E19, 7.428161E-7, 
        4.766499E19, 1.405762E-6, 4.784913E19, 1.090323E-5, 4.918388E19, 
        2.050093E-5, 5.164678E19, 1.141528E-6, 5.179294E19, 6.302739E-7, 
        5.186414E19, 1.166297E-5, 5.328083E19, 2.248001E-5, 5.594976E19, 
        1.852149E-8, 5.596043E19, 1.329684E-7, 5.596507E19, 1.507823E-9, 
        5.596619E19, 5.301304E-7, 5.604303E19, 1.119985E-5, 5.754166E19, 
        2.496909E-5, 6.042151E19, 1.451153E-9, 6.042458E19, 2.378372E-7, 
        6.045956E19, 7.834350E-7, 6.052821E19, 2.534995E-6, 6.089311E19, 
        9.038908E-6, 6.196637E19, 2.562149E-5, 6.507503E19, 1.765163E-7, 
        6.508274E19, 3.699785E-7, 6.512461E19, 1.842793E-6, 6.539976E19, 
        9.040771E-6, 6.655495E19, 2.725465E-5, 6.988521E19, 3.185849E-7, 
        6.995174E19, 5.933457E-7, 6.996866E19, 7.823664E-7, 7.007028E19, 
        1.092899E-5, 7.130741E19, 2.903186E-5, 7.487609E19, 2.849820E-7, 
        7.490468E19, 1.084841E-5, 7.622375E19, 3.192676E-5, 8.003767E19, 
        1.597977E-7, 8.006671E19, 1.591384E-7, 8.007674E19, 4.408306E-6, 
        8.055522E19, 5.531786E-6 };
      THEN( "bragg edges vector is correctly output" ){
        checkPartOfVec(b,braggVals_0_to_99);
        checkPartOfVec(b,braggVals_500_to_599,500);
        checkPartOfVec(b,braggVals_600_to_699,600);
        restAreZero(nw,b);
        REQUIRE( nw   == 740 );
        REQUIRE( maxb == 612 );

      } // THEN
    } // WHEN

    WHEN( "3 principal scattering atoms in compound" ){
      npr = 3;
      auto output = coherentElastic( iel, npr, b, emax );
      int nw   = std::get<0>(output);
      int maxb = std::get<1>(output);

      braggVals_0_to_99 = { 1.06117E-3, 0.00000000, 3.75436E-3, 6.35237E-3, 
      4.24469E-3, 4.06388E-3, 4.81554E-3, 5.51318E-3, 7.99906E-3, 4.44054E-3, 
      9.55056E-3, 0.00000000, 1.12631E-2, 1.46701E-2, 1.23242E-2, 0.00000000, 
      1.33049E-2, 1.73418E-2, 1.50174E-2, 3.17618E-3, 1.55078E-2, 1.27567E-2, 
      1.60786E-2, 3.01717E-3, 1.69787E-2, 8.10707E-5, 1.92621E-2, 2.86156E-3, 
      2.07331E-2, 1.10046E-4, 2.08136E-2, 0.00000000, 2.45680E-2, 1.27619E-2, 
      2.62805E-2, 4.80194E-3, 2.65293E-2, 0.00000000, 2.73417E-2, 4.62745E-3, 
      2.82418E-2, 3.77156E-4, 3.02837E-2, 1.14946E-2, 3.05252E-2, 4.54628E-3, 
      3.19962E-2, 8.85847E-5, 3.37893E-2, 8.46983E-3, 3.48504E-2, 0.00000000, 
      3.58311E-2, 2.11349E-2, 3.77924E-2, 0.00000000, 3.80340E-2, 8.14573E-3, 
      3.82022E-2, 1.35462E-3, 4.15468E-2, 9.81368E-3, 4.19566E-2, 1.93890E-3, 
      4.32593E-2, 1.52369E-4, 4.33398E-2, 0.00000000, 4.50524E-2, 7.33509E-3, 
      4.61135E-2, 0.00000000, 4.88067E-2, 3.52366E-3, 4.92971E-2, 7.15493E-3, 
      4.94653E-2, 7.14275E-3, 4.98679E-2, 3.42644E-3, 5.07680E-2, 2.81302E-4, 
      5.19975E-2, 2.25790E-2, 5.46029E-2, 0.00000000, 5.57518E-2, 1.81811E-2, 
      6.00698E-2, 1.58809E-3, 6.03186E-2, 1.80185E-3, 6.43145E-2, 4.81757E-3, 
      6.70150E-2, 1.47787E-3, 6.79151E-2, 9.57251E-3, 7.13329E-2, 2.91467E-3 };
      THEN( "bragg edges vector is correctly output" ){
        checkVec(b,braggVals_0_to_99,braggVals_0_to_99.size());
        restAreZero(nw,b);
        REQUIRE( nw   == 740 );
        REQUIRE( maxb == 612 );
      } // THEN
    } // WHEN
  } // GIVEN



  GIVEN( "aluminum is the requested material" ){
    iel = 4; npr = 1;
    auto output = coherentElastic( iel, npr, b, emax );
    int nw   = std::get<0>(output);
    int maxb = std::get<1>(output);

    braggVals_0_to_99 = { 3.759016E-3, 5.508122E-3, 5.012021E-3, 3.57763E-3, 
    1.00240E-2, 5.05953E-3, 1.37830E-2, 8.62956E-3, 1.50360E-2, 2.75406E-3, 
    2.00480E-2, 1.78881E-3, 2.38071E-2, 6.56611E-3, 2.50601E-2, 6.39986E-3, 
    3.00721E-2, 5.84224E-3, 3.38311E-2, 7.34416E-3, 4.00961E-2, 2.52976E-3, 
    4.38551E-2, 9.67567E-3, 4.51081E-2, 5.96271E-3, 5.01202E-2, 4.52538E-3, 
    5.38792E-2, 4.36466E-3, 5.51322E-2, 4.31478E-3, 6.01442E-2, 1.37703E-3, 
    6.39032E-2, 8.01549E-3, 6.51562E-2, 3.96902E-3, 7.01683E-2, 7.64929E-3, 
    7.39273E-2, 1.11784E-2, 8.01923E-2, 8.94407E-4, 8.39513E-2, 3.49661E-3, 
    8.52043E-2, 6.94162E-3, 9.02163E-2, 5.05953E-3, 9.39754E-2, 7.71137E-3, 
    9.52284E-2, 3.28305E-3, 1.00240E-1, 3.19993E-3, 1.03999E-1, 9.42470E-3, 
    1.05252E-1, 6.24562E-3, 1.10264E-1, 3.05101E-3, 1.14023E-1, 6.00060E-3, 
    1.20288E-1, 2.92112E-3, 1.24047E-1, 8.62956E-3, 1.25300E-1, 3.57763E-3, 
    1.30312E-1, 8.41957E-3, 1.34071E-1, 8.30070E-3, 1.35324E-1, 3.67208E-3, 
    1.44095E-1, 5.33785E-3, 1.45348E-1, 7.97219E-3, 1.50360E-1, 5.22546E-3, 
    1.54119E-1, 5.16134E-3, 1.60384E-1, 1.26488E-3, 1.64143E-1, 1.25031E-2, 
    1.65396E-1, 4.98228E-3, 1.70408E-1, 4.90846E-3, 1.74167E-1, 7.28281E-3, 
    1.75420E-1, 4.83783E-3, 1.80432E-1, 2.98135E-3, 1.84191E-1, 5.50812E-3 };

    braggVals_500_to_599 = { 9.35995E-1, 5.235946E-3, 9.372480E-1, 7.848667E-4, 
    9.46019E-1, 5.20813E-3, 9.47272E-1, 4.16374E-3, 9.52284E-1, 1.03819E-3, 
    9.56043E-1, 2.07230E-3, 9.66067E-1, 2.83459E-3, 9.67320E-1, 7.72570E-4, 
    9.72332E-1, 3.85288E-3, 9.76091E-1, 4.10182E-3, 9.77344E-1, 1.53719E-3, 
    9.82356E-1, 7.66635E-4, 9.86115E-1, 1.78540E-3, 9.87368E-1, 3.05874E-3, 
    9.92380E-1, 2.79676E-3, 9.96139E-1, 1.52262E-3, 1.00240429, 1.26488E-3, 
    1.00616331, 4.04005E-3, 1.00741631, 2.77581E-3, 1.01242833, 1.51032E-3, 
    1.01618735, 3.26631E-3, 1.01744035, 1.50660E-3, 1.02245238, 7.51452E-4, 
    1.02621139, 4.00040E-3, 1.02746440, 1.49923E-3, 1.03247642, 4.48677E-3, 
    1.03623543, 2.73694E-3, 1.04250046, 4.96128E-4, 1.04625948, 1.48570E-3, 
    1.04751248, 3.46457E-3, 1.05252450, 9.87519E-4, 1.05628352, 2.95728E-3, 
    1.05753653, 9.85177E-4, 1.06254855, 1.47427E-3, 1.06630756, 4.41502E-3, 
    1.06756057, 1.47081E-3, 1.07257259, 7.33685E-4, 1.07633161, 2.68547E-3, 
    1.08259663, 1.21713E-3, 1.08635565, 2.43005E-3, 1.08760865, 1.45719E-3, 
    1.09262068, 1.93846E-3, 1.09637969, 2.90270E-3, 1.09763270, 1.45052E-3, 
    1.10640373, 1.44476E-3, 1.10765674, 2.40657E-3, 1.11266876, 2.88137E-3, 
    1.11642778, 3.11623E-3, 1.12269280, 4.78081E-4, 1.12645182, 4.29554E-3 };

    braggVals_2000_to_2099={ 2.394595E18, 1.03396E-8, 2.394595E18, 1.03396E-8,
    2.394595E18, 1.03396E-8, 2.394595E18, 1.03396E-8, 2.394595E18, 1.03396E-8, 
    2.394595E18, 1.03396E-8, 2.394595E18, 1.03396E-8, 2.394595E18, 1.03396E-8, 
    2.394595E18, 1.03396E-8, 2.394595E18, 1.03396E-8, 2.394595E18, 1.03396E-8, 
    2.394595E18, 1.03396E-8, 2.394595E18, 1.03396E-8, 2.394595E18, 1.03396E-8, 
    2.394595E18, 1.03396E-8, 2.394595E18, 1.03396E-8, 2.394595E18, 1.03396E-8, 
    2.394595E18, 1.03396E-8, 2.394595E18, 1.03396E-8, 2.394595E18, 1.03396E-8, 
    2.394595E18, 1.03396E-8, 2.394595E18, 1.03396E-8, 2.394595E18, 1.03396E-8, 
    2.394595E18, 1.03396E-8, 2.394595E18, 1.03396E-8, 2.394595E18, 1.03396E-8, 
    2.394595E18, 1.03396E-8, 2.394595E18, 1.03396E-8, 2.394595E18, 1.03396E-8, 
    2.394595E18, 1.03396E-8, 2.394595E18, 1.03396E-8, 2.394595E18, 1.03396E-8, 
    2.394595E18, 1.03396E-8, 2.394595E18, 1.03396E-8, 2.394595E18, 1.03396E-8, 
    2.394595E18, 1.03396E-8, 2.418783E18, 1.028777E-8,2.418783E18, 1.028777E-8, 
    2.418783E18, 1.028777E-8,2.418783E18, 1.028777E-8,2.418783E18, 1.028777E-8, 
    2.418783E18, 1.028777E-8,2.418783E18, 1.028777E-8,2.418783E18, 1.028777E-8, 
    2.418783E18, 1.028777E-8,2.418783E18, 1.028777E-8,2.418783E18, 1.028777E-8, 
    2.418783E18, 1.028777E-8,2.418783E18, 1.028777E-8,2.418783E18, 1.028777E-8};

    THEN( "bragg edges vector is correctly output" ){
      checkPartOfVec(b,braggVals_0_to_99);
      checkPartOfVec(b,braggVals_500_to_599,500);
      checkPartOfVec(b,braggVals_2000_to_2099,2000);
      restAreZero(nw,b);
      REQUIRE( nw   == 59582 );
      REQUIRE( maxb == 1234  );

    } // THEN
  } // GIVEN




  GIVEN( "lead is the requested material" ){
    iel = 5; npr = 1;
    auto output = coherentElastic( iel, npr, b, emax );
    int nw   = std::get<0>(output);
    int maxb = std::get<1>(output);


    braggVals_0_to_99 = { 2.51410E-3, 2.46417E-3, 3.35213E-3, 1.60052E-3, 
      6.70427E-3, 2.26348E-3, 9.21837E-3, 3.86061E-3, 1.00564E-2, 1.23208E-3, 
      1.34085E-2, 8.00263E-4, 1.59226E-2, 2.93748E-3, 1.67606E-2, 2.86311E-3, 
      2.01128E-2, 2.61365E-3, 2.26269E-2, 3.28556E-3, 2.68170E-2, 1.13174E-3, 
      2.93311E-2, 4.32861E-3, 3.01692E-2, 2.66754E-3, 3.35213E-2, 2.02452E-3, 
      3.60354E-2, 1.95262E-3, 3.68735E-2, 1.93030E-3, 4.02256E-2, 6.16043E-4, 
      4.27397E-2, 3.58589E-3, 4.35777E-2, 1.77562E-3, 4.69299E-2, 3.42207E-3, 
      4.94440E-2, 5.00090E-3, 5.36341E-2, 4.00131E-4, 5.61482E-2, 1.56428E-3, 
      5.69863E-2, 3.10547E-3, 6.03384E-2, 2.26348E-3, 6.28525E-2, 3.44984E-3, 
      6.36906E-2, 1.46874E-3, 6.70427E-2, 1.43155E-3, 6.95568E-2, 4.21633E-3, 
      7.03948E-2, 2.79410E-3, 7.37470E-2, 1.36493E-3, 7.62611E-2, 2.68449E-3, 
      8.04512E-2, 1.30682E-3, 8.29653E-2, 3.86061E-3, 8.38034E-2, 1.60052E-3, 
      8.71555E-2, 3.76667E-3, 8.96696E-2, 3.71349E-3, 9.05076E-2, 1.64278E-3, 
      9.63739E-2, 2.38799E-3, 9.72119E-2, 3.56652E-3, 1.00564E-1, 2.33772E-3, 
      1.03078E-1, 2.30903E-3, 1.07268E-1, 5.65871E-4, 1.09782E-1, 5.59354E-3, 
      1.10620E-1, 2.22892E-3, 1.13972E-1, 2.19590E-3, 1.16486E-1, 3.25811E-3, 
      1.17324E-1, 2.16430E-3, 1.20676E-1, 1.33377E-3, 1.23191E-1, 2.46417E-3 };
    braggVals_500_to_599 = { 6.26011E-1, 2.34240E-3, 6.26849E-1, 3.51126E-4, 
      6.32715E-1, 2.32996E-3, 6.33553E-1, 1.86273E-3, 6.36906E-1, 4.64457E-4, 
      6.39420E-1, 9.27087E-4, 6.46124E-1, 1.26811E-3, 6.46962E-1, 3.45625E-4, 
      6.50314E-1, 1.72366E-3, 6.52828E-1, 1.83503E-3, 6.53666E-1, 6.87696E-4, 
      6.57018E-1, 3.42970E-4, 6.59532E-1, 7.98737E-4, 6.60370E-1, 1.36839E-3, 
      6.63723E-1, 1.25118E-3, 6.66237E-1, 6.81178E-4, 6.70427E-1, 5.65871E-4, 
      6.72941E-1, 1.80740E-3, 6.73779E-1, 1.24181E-3, 6.77131E-1, 6.75676E-4, 
      6.79645E-1, 1.46125E-3, 6.80483E-1, 6.74010E-4, 6.83835E-1, 3.36178E-4, 
      6.86350E-1, 1.78966E-3, 6.87188E-1, 6.70714E-4, 6.90540E-1, 2.00725E-3, 
      6.93054E-1, 1.22442E-3, 6.97244E-1, 2.21953E-4, 6.99758E-1, 6.64662E-4, 
      7.00596E-1, 1.54995E-3, 7.03948E-1, 4.41787E-4, 7.06462E-1, 1.32300E-3, 
      7.07300E-1, 4.40739E-4, 7.10653E-1, 6.59548E-4, 7.13167E-1, 1.97515E-3, 
      7.14005E-1, 6.57998E-4, 7.17357E-1, 3.28229E-4, 7.19871E-1, 1.20140E-3, 
      7.24061E-1, 5.44510E-4, 7.26575E-1, 1.08713E-3, 7.27413E-1, 6.51905E-4, 
      7.30765E-1, 8.67211E-4, 7.33279E-1, 1.29858E-3, 7.34117E-1, 6.48921E-4, 
      7.39984E-1, 6.46344E-4, 7.40822E-1, 1.07663E-3, 7.44174E-1, 1.28904E-3, 
      7.46688E-1, 1.39411E-3, 7.50878E-1, 2.13879E-4, 7.53392E-1, 1.92170E-3 };
    braggVals_2000_to_2099 = { 1.601551E18, 1.264298E-8, 1.601551E18, 
      1.264298E-8, 1.601551E18, 1.264298E-8, 1.601551E18, 1.264298E-8, 
      1.601551E18, 1.264298E-8, 1.601551E18, 1.264298E-8, 1.601551E18, 
      1.264298E-8, 1.601551E18, 1.264298E-8, 1.601551E18, 1.264298E-8, 
      1.601551E18, 1.264298E-8, 1.601551E18, 1.264298E-8, 1.601551E18, 
      1.264298E-8, 1.601551E18, 1.264298E-8, 1.601551E18, 1.264298E-8, 
      1.601551E18, 1.264298E-8, 1.601551E18, 1.264298E-8, 1.601551E18, 
      1.264298E-8, 1.601551E18, 1.264298E-8, 1.601551E18, 1.264298E-8, 
      1.601551E18, 1.264298E-8, 1.601551E18, 1.264298E-8, 1.601551E18, 
      1.264298E-8, 1.601551E18, 1.264298E-8, 1.601551E18, 1.264298E-8, 
      1.601551E18, 1.264298E-8, 1.601551E18, 1.264298E-8, 1.601551E18, 
      1.264298E-8, 1.601551E18, 1.264298E-8, 1.601551E18, 1.264298E-8, 
      1.601551E18, 1.264298E-8, 1.601551E18, 1.264298E-8, 1.601551E18, 
      1.264298E-8, 1.601551E18, 1.264298E-8, 1.601551E18, 1.264298E-8, 
      1.601551E18, 1.264298E-8, 1.601551E18, 1.264298E-8, 1.617729E18, 
      1.257960E-8, 1.617729E18, 1.257960E-8, 1.617729E18, 1.257960E-8, 
      1.617729E18, 1.257960E-8, 1.617729E18, 1.257960E-8, 1.617729E18, 
      1.257960E-8, 1.617729E18, 1.257960E-8, 1.617729E18, 1.257960E-8, 
      1.617729E18, 1.257960E-8, 1.617729E18, 1.257960E-8, 1.617729E18, 
      1.257960E-8, 1.617729E18, 1.257960E-8, 1.617729E18, 1.257960E-8, 
      1.617729E18, 1.257960E-8 };

    THEN( "bragg edges vector is correctly output" ){
      checkPartOfVec(b,braggVals_0_to_99);
      checkPartOfVec(b,braggVals_500_to_599,500);
      checkPartOfVec(b,braggVals_2000_to_2099,2000);
      restAreZero(nw,b);
      REQUIRE( nw   == 59582 );
      REQUIRE( maxb == 1234  );

    } // THEN
  } // GIVEN




  GIVEN( "iron is the requested material" ){
    iel = 6; npr = 1;
    auto output = coherentElastic( iel, npr, b, emax );
    int nw   = std::get<0>(output);
    int maxb = std::get<1>(output);
    braggVals_0_to_99 = { 5.00050E-3, 8.71143E-2, 1.00010E-2, 3.07995E-2, 
      1.50015E-2, 1.00590E-1, 2.00020E-2, 4.35571E-2, 2.50025E-2, 7.79174E-2, 
      3.00030E-2, 2.37095E-2, 3.50035E-2, 1.31704E-1, 4.00040E-2, 1.53997E-2, 
      4.50045E-2, 8.71143E-2, 5.00050E-2, 5.50959E-2, 5.50055E-2, 5.25319E-2, 
      6.00060E-2, 5.02954E-2, 6.50065E-2, 1.44967E-1, 7.50075E-2, 8.99713E-2, 
      8.00080E-2, 2.17785E-2, 8.50085E-2, 8.45133E-2, 9.00090E-2, 5.13326E-2, 
      9.50095E-2, 1.19912E-1, 1.00010E-1, 3.89587E-2, 1.05010E-1, 7.60396E-2, 
      1.10011E-1, 3.71456E-2, 1.15011E-1, 7.26583E-2, 1.20012E-1, 1.18547E-2, 
      1.25012E-1, 1.21960E-1, 1.30013E-1, 3.41690E-2, 1.35013E-1, 1.34121E-1, 
      1.40014E-1, 6.58522E-2, 1.45014E-1, 3.23534E-2, 1.55015E-1, 1.25169E-1, 
      1.60016E-1, 7.69989E-3, 1.65016E-1, 1.21317E-1, 1.70017E-1, 5.97599E-2, 
      1.75017E-1, 5.89000E-2, 1.80018E-1, 4.35571E-2, 1.85018E-1, 1.43215E-1, 
      1.90019E-1, 2.82636E-2, 1.95019E-1, 5.57978E-2, 2.00020E-1, 2.75479E-2, 
      2.05020E-1, 5.44199E-2, 2.10021E-1, 5.37681E-2, 2.15021E-1, 1.32848E-1, 
      2.20022E-1, 2.62659E-2, 2.25022E-1, 1.29862E-1, 2.35023E-1, 1.01655E-1, 
      2.40024E-1, 2.51477E-2, 2.45024E-1, 1.12004E-1, 2.50025E-1, 3.07995E-2, 
      2.55025E-1, 4.87937E-2, 2.60026E-1, 7.24835E-2, 2.65026E-1, 7.17964E-2 };
    braggVals_500_to_599 = { 1.37013750, 1.84196E-2, 1.37513800, 1.31329E-2, 
      1.38013850, 1.04873E-2, 1.38513900, 5.23419E-2, 1.39013951, 7.83715E-3, 
      1.39514001, 4.17231E-2, 1.40014051, 5.20607E-3, 1.40514101, 7.79520E-3, 
      1.41014151, 1.03751E-2, 1.41514201, 5.43733E-2, 1.42014252, 5.16928E-3, 
      1.42514302, 1.54806E-2, 1.43514402, 2.05687E-2, 1.44014452, 5.13326E-3, 
      1.44514503, 3.58706E-2, 1.45014553, 1.53465E-2, 1.45514603, 1.78735E-2, 
      1.46014653, 2.03919E-2, 1.46514703, 3.05356E-2, 1.47014753, 5.92738E-3, 
      1.47514804, 2.02879E-2, 1.48014854, 1.01268E-2, 1.48514904, 4.80214E-2, 
      1.49014954, 3.02783E-2, 1.49515004, 1.00759E-2, 1.50015055, 7.54432E-3, 
      1.50515105, 4.51906E-2, 1.51515205, 1.00091E-2, 1.52015255, 7.49452E-3, 
      1.52515305, 9.97630E-3, 1.53015356, 2.24099E-2, 1.53515406, 3.23172E-2, 
      1.54015456, 9.92760E-3, 1.54515506, 9.91152E-3, 1.55015556, 9.89552E-3, 
      1.55515607, 9.87960E-3, 1.56515707, 4.18539E-2, 1.57015757, 4.91614E-3, 
      1.57515807, 3.43583E-2, 1.58015857, 4.90056E-3, 1.58515908, 1.95713E-2, 
      1.59516008, 3.41422E-2, 1.60516108, 3.40357E-2, 1.61016159, 1.94187E-2, 
      1.61516209, 1.45415E-2, 1.62016259, 9.67937E-3, 1.62516309, 3.86578E-2, 
      1.63516410, 1.92697E-2, 1.64016460, 7.21512E-3, 1.64516510, 4.80276E-3 };
    braggVals_2000_to_2099 = { 2.992398E18, 2.312332E-9, 2.992398E18, 
      2.312332E-9, 2.992398E18, 2.312332E-9, 2.992398E18, 2.312332E-9, 
      2.992398E18, 2.312332E-9, 2.992398E18, 2.312332E-9, 2.992398E18, 
      2.312332E-9, 2.992398E18, 2.312332E-9, 2.992398E18, 2.312332E-9, 
      2.992398E18, 2.312332E-9, 2.992398E18, 2.312332E-9, 2.992398E18, 
      2.312332E-9, 2.992398E18, 2.312332E-9, 2.992398E18, 2.312332E-9, 
      2.992398E18, 2.312332E-9, 2.992398E18, 2.312332E-9, 2.992398E18, 
      2.312332E-9, 2.992398E18, 2.312332E-9, 2.992398E18, 2.312332E-9, 
      2.992398E18, 2.312332E-9, 2.992398E18, 2.312332E-9, 2.992398E18, 
      2.312332E-9, 2.992398E18, 2.312332E-9, 2.992398E18, 2.312332E-9, 
      2.992398E18, 2.312332E-9, 2.992398E18, 2.312332E-9, 2.992398E18, 
      2.312332E-9, 2.992398E18, 2.312332E-9, 2.992398E18, 2.312332E-9, 
      2.992398E18, 2.312332E-9, 2.992398E18, 2.312332E-9, 2.992398E18, 
      2.312332E-9, 2.992398E18, 2.312332E-9, 2.992398E18, 2.312332E-9, 
      2.992398E18, 2.312332E-9, 2.992398E18, 2.312332E-9, 2.992398E18, 
      2.312332E-9, 2.992398E18, 2.312332E-9, 2.992398E18, 2.312332E-9, 
      2.992398E18, 2.312332E-9, 2.992398E18, 2.312332E-9, 2.992398E18, 
      2.312332E-9, 2.992398E18, 2.312332E-9, 2.992398E18, 2.312332E-9, 
      2.992398E18, 2.312332E-9, 2.992398E18, 2.312332E-9, 2.992398E18, 
      2.312332E-9, 2.992398E18, 2.312332E-9, 2.992398E18, 2.312332E-9, 
      2.992398E18, 2.312332E-9 };

    THEN( "bragg edges vector is correctly output" ){
      checkPartOfVec(b,braggVals_0_to_99);
      checkPartOfVec(b,braggVals_500_to_599,500);
      checkPartOfVec(b,braggVals_2000_to_2099,2000);

      restAreZero(nw,b);
      REQUIRE( nw   == 59246 );
      REQUIRE( maxb == 1282  );


    } // THEN
  } // GIVEN
  

} // TEST CASE