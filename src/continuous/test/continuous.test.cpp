#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "continuous/continuous.h"
#include "generalTools/testing.h"

TEST_CASE( "continuous treatment" ){
  int nphon;
  double delta, tbeta, tev, sc, scaling, lambda_s, t_bar;
  std::vector<double> alpha, beta, rho, sabCorrect;
  std::tuple<double,double> output;

  GIVEN( "simplified water example" ){
    nphon = 100; 
    tbeta = 0.444444; tev = 2.5507297688E-2; 
    delta = 0.00255/tev; 
    sc = 0.0253/tev; scaling = sc;

    rho   = { 0.0, 0.0005, 0.001, 0.002, 0.0035, 0.005, 0.0075, 0.01, 0.013, 
    0.0165, 0.02, 0.0245, 0.029, 0.034, 0.0395, 0.045, 0.0506, 0.0562, 0.0622, 
    0.0686, 0.075, 0.083, 0.091, 0.099, 0.107, 0.115, 0.1197, 0.1214, 0.1218, 
    0.1195, 0.1125, 0.1065, 0.1005, 0.09542, 0.09126, 0.0871, 0.0839, 0.0807, 
    0.07798, 0.07574, 0.0735, 0.07162, 0.06974, 0.06804, 0.06652, 0.065, 0.0634, 
    0.0618, 0.06022, 0.05866, 0.0571, 0.05586, 0.05462, 0.0535, 0.0525, 0.0515, 
    0.05042, 0.04934, 0.04822, 0.04706, 0.0459, 0.04478, 0.04366, 0.04288, 
    0.04244, 0.042, 0.};


    WHEN( "alpha and beta grid is relatively small" ){
      alpha = {0.01, 0.08, 4.00, 6.00, 20, 50};
      beta = {0.0, 0.4, 1.0, 5.0, 9.0, 20};
      for ( auto& a : alpha ){ a *= scaling; }
      for ( auto& b : beta  ){ b *= sc;      }
      std::vector<double> sab( alpha.size()*beta.size(), 0.0 );

      output = continuous( nphon, delta, tbeta, rho, alpha, beta, sab );
      THEN( "continuous output matches expected value" ){
        lambda_s =   0.2352041964494244, t_bar = 1.9344942054735312;
        sabCorrect = {5.6553448E-4, 2.9998239E-4, 3.5728611E-4, 1.3283598E-4, 
        7.2155584E-8, 7.038367E-16, 4.4601537E-3, 2.3709915E-3, 2.8246954E-3, 
        1.0657726E-3, 4.6129213E-6, 2.919488E-12, 1.0134916E-1, 6.0785511E-2, 
        7.3518778E-2, 5.1338891E-2, 9.9527483E-3, 2.6923598E-5, 1.0257124E-1, 
        6.5337140E-2, 7.9600320E-2, 6.9281786E-2, 1.9855721E-2, 1.4530853E-4, 
        2.6874094E-2, 2.4511956E-2, 3.1123242E-2, 6.6847691E-2, 6.6284139E-2, 
        1.1948925E-2, 1.0067245E-3, 1.2017244E-3, 1.5855351E-3, 6.8388574E-3, 
        1.7465674E-2, 4.3725435E-2};
        
        REQUIRE(ranges::equal(sabCorrect, sab, equal));
        REQUIRE(lambda_s == Approx(std::get<0>(output)).epsilon(1e-6));
        REQUIRE(t_bar    == Approx(std::get<1>(output)).epsilon(1e-6));
      } // THEN
    } // WHEN

  
    WHEN( "alpha and beta grid is larger" ){
      sabCorrect = { 5.6553448E-4, 5.6838546E-4, 5.7123645E-4, 5.7408744E-4, 
      5.7693842E-4, 5.7978941E-4, 5.8264040E-4, 5.8549139E-4, 5.8834237E-4, 
      5.9119336E-4, 5.9404435E-4, 3.1648544E-4, 2.9150645E-4, 2.9998239E-4, 
      2.8804518E-4, 3.1222006E-4, 3.2063479E-4, 3.3305120E-4, 3.4840621E-4, 
      3.5728611E-4, 4.8899523E-4, 2.1589873E-4, 8.8731671E-5, 1.2984180E-7, 
      3.0467685E-8, 5.463617E-12, 7.038367E-16, 8.790798E-19, 1.1287644E-3, 
      1.1344539E-3, 1.1401434E-3, 1.1458329E-3, 1.1515224E-3, 1.1572119E-3, 
      1.1629014E-3, 1.1685909E-3, 1.1742804E-3, 1.1799699E-3, 1.1856594E-3, 
      6.3184092E-4, 5.8200398E-4, 5.9892822E-4, 5.7512420E-4, 6.2337635E-4, 
      6.4018338E-4, 6.6497494E-4, 6.9563024E-4, 7.1336611E-4, 9.7634697E-4, 
      4.3157066E-4, 1.7762681E-4, 5.1895921E-7, 1.2203981E-7, 4.377827E-11, 
      1.128190E-14, 1.523015E-17, 2.2483397E-3, 2.2596689E-3, 2.2709982E-3, 
      2.2823275E-3, 2.2936568E-3, 2.3049861E-3, 2.3163153E-3, 2.3276446E-3, 
      2.3389739E-3, 2.3503032E-3, 2.3616325E-3, 1.2591745E-3, 1.1599831E-3, 
      1.1937214E-3, 1.1463938E-3, 1.2425089E-3, 1.2760334E-3, 1.3254520E-3, 
      1.3865443E-3, 1.4219205E-3, 1.9461373E-3, 8.6223168E-4, 3.5590311E-4, 
      2.0725741E-6, 4.8950451E-7, 3.513314E-10, 1.811652E-13, 2.811790E-16, 
      3.3587836E-3, 3.3757033E-3, 3.3926229E-3, 3.4095426E-3, 3.4264623E-3, 
      3.4433820E-3, 3.4603016E-3, 3.4772213E-3, 3.4941410E-3, 3.5110607E-3, 
      3.5279803E-3, 1.8820258E-3, 1.7339591E-3, 1.7844018E-3, 1.7138290E-3, 
      1.8574203E-3, 1.9075730E-3, 1.9814548E-3, 2.0727671E-3, 2.1256885E-3, 
      2.9094057E-3, 1.2919782E-3, 5.3482022E-4, 4.6559536E-6, 1.1043903E-6, 
      1.1894553E-9, 9.204526E-13, 1.614776E-15, 4.4601537E-3, 4.4826147E-3, 
      4.5050757E-3, 4.5275368E-3, 4.5499978E-3, 4.5724588E-3, 4.5949198E-3, 
      4.6173809E-3, 4.6398419E-3, 4.6623029E-3, 4.6847639E-3, 2.5004198E-3, 
      2.3039535E-3, 2.3709915E-3, 2.2774497E-3, 2.4681329E-3, 2.5348248E-3, 
      2.6330071E-3, 2.7543235E-3, 2.8246954E-3, 3.8661866E-3, 1.7208054E-3, 
      7.1436949E-4, 8.2642113E-6, 1.9686645E-6, 2.8282035E-9, 2.919488E-12, 
      5.712788E-15, 1.0881003E-2, 1.0935701E-2, 1.0990399E-2, 1.1045096E-2, 
      1.1099794E-2, 1.1154492E-2, 1.1209190E-2, 1.1263887E-2, 1.1318585E-2, 
      1.1373283E-2, 1.1427981E-2, 6.1185582E-3, 5.6415024E-3, 5.8058553E-3, 
      5.5801755E-3, 6.0454756E-3, 6.2095539E-3, 6.4501663E-3, 6.7470430E-3, 
      6.9201474E-3, 9.4725406E-3, 4.2742013E-3, 1.8044626E-3, 5.1163131E-5, 
      1.2498055E-5, 4.4997042E-8, 1.164510E-10, 3.697511E-13, 2.0893706E-2, 
      2.0998425E-2, 2.1103145E-2, 2.1207864E-2, 2.1312583E-2, 2.1417302E-2, 
      2.1522021E-2, 2.1626740E-2, 2.1731459E-2, 2.1836178E-2, 2.1940897E-2, 
      1.1808368E-2, 1.0899575E-2, 1.1217778E-2, 1.0792543E-2, 1.1686378E-2, 
      1.2005869E-2, 1.2471393E-2, 1.3044400E-2, 1.3381381E-2, 1.8319834E-2, 
      8.4533348E-3, 3.6662678E-3, 2.0141089E-4, 5.1210439E-5, 3.7029727E-7, 
      1.9256303E-9, 1.005148E-11, 3.0091653E-2, 3.0242029E-2, 3.0392405E-2, 
      3.0542780E-2, 3.0693156E-2, 3.0843532E-2, 3.0993907E-2, 3.1144283E-2, 
      3.1294659E-2, 3.1445034E-2, 3.1595410E-2, 1.7092811E-2, 1.5794511E-2, 
      1.6256612E-2, 1.5655926E-2, 1.6943809E-2, 1.7410363E-2, 1.8085905E-2, 
      1.8915446E-2, 1.9407426E-2, 2.6574283E-2, 1.2533371E-2, 5.5773754E-3, 
      4.4592095E-4, 1.1776220E-4, 1.2828083E-6, 1.0053103E-8, 7.320599E-11, 
      3.8525392E-2, 3.8717352E-2, 3.8909312E-2, 3.9101272E-2, 3.9293232E-2, 
      3.9485192E-2, 3.9677152E-2, 3.9869112E-2, 4.0061072E-2, 4.0253032E-2, 
      4.0444992E-2, 2.1994104E-2, 2.0345619E-2, 2.0942189E-2, 2.0188268E-2, 
      2.1837860E-2, 2.2443438E-2, 2.3314872E-2, 2.4382442E-2, 2.5020883E-2, 
      3.4266742E-2, 1.6510859E-2, 7.5302079E-3, 7.7993141E-4, 2.1352291E-4, 
      3.1149725E-6, 3.2699901E-8, 3.061132E-10, 7.5571073E-2, 7.5941293E-2, 
      7.6311513E-2, 7.6681733E-2, 7.7051953E-2, 7.7422173E-2, 7.7792393E-2, 
      7.8162613E-2, 7.8532833E-2, 7.8903053E-2, 7.9273272E-2, 4.4466407E-2, 
      4.1400373E-2, 4.2631833E-2, 4.1334009E-2, 4.4580444E-2, 4.5867650E-2, 
      4.7656406E-2, 4.9817834E-2, 5.1173555E-2, 7.0169484E-2, 3.8091714E-2, 
      1.9752726E-2, 4.4062956E-3, 1.4501293E-3, 5.4372146E-5, 1.4664119E-6, 
      3.2543351E-8, 1.0134916E-1, 1.0183298E-1, 1.0231680E-1, 1.0280061E-1, 
      1.0328443E-1, 1.0376825E-1, 1.0425206E-1, 1.0473588E-1, 1.0521970E-1, 
      1.0570351E-1, 1.0618733E-1, 6.2686865E-2, 5.8982414E-2, 6.0785511E-2, 
      5.9465001E-2, 6.3851421E-2, 6.5809924E-2, 6.8397123E-2, 7.1457267E-2, 
      7.3518778E-2, 1.0107725E-1, 6.5059818E-2, 3.9786012E-2, 1.4739589E-2, 
      6.0815621E-3, 4.7811040E-4, 2.6923598E-5, 1.2081963E-6, 1.0257124E-1, 
      1.0304994E-1, 1.0352864E-1, 1.0400735E-1, 1.0448605E-1, 1.0496476E-1, 
      1.0544346E-1, 1.0592217E-1, 1.0640087E-1, 1.0687958E-1, 1.0735828E-1, 
      6.6637858E-2, 6.3340739E-2, 6.5337140E-2, 6.4443403E-2, 6.8927079E-2, 
      7.1156412E-2, 7.3978649E-2, 7.7252714E-2, 7.9600320E-2, 1.0980598E-1, 
      8.1575092E-2, 5.6715063E-2, 2.7445910E-2, 1.3387084E-2, 1.6527252E-3, 
      1.4530853E-4, 1.0039192E-5, 9.2892662E-2, 9.3318039E-2, 9.3743416E-2, 
      9.4168793E-2, 9.4594170E-2, 9.5019547E-2, 9.5444924E-2, 9.5870301E-2, 
      9.6295678E-2, 9.6721055E-2, 9.7146432E-2, 6.3319048E-2, 6.0775306E-2, 
      6.2754227E-2, 6.2359164E-2, 6.6470184E-2, 6.8722884E-2, 7.1474553E-2, 
      7.4611684E-2, 7.6986816E-2, 1.0662339E-1, 8.9569333E-2, 6.9212387E-2, 
      4.0050446E-2, 2.2292146E-2, 3.8326351E-3, 4.6638109E-4, 4.4189172E-5, 
      7.9437473E-2, 7.9795769E-2, 8.0154065E-2, 8.0512361E-2, 8.0870657E-2, 
      8.1228953E-2, 8.1587249E-2, 8.1945545E-2, 8.2303841E-2, 8.2662137E-2, 
      8.3020433E-2, 5.6729258E-2, 5.4956670E-2, 5.6807106E-2, 5.6832227E-2, 
      6.0399364E-2, 6.2532610E-2, 6.5061439E-2, 6.7900303E-2, 7.0153630E-2, 
      9.7598494E-2, 9.1219651E-2, 7.7135468E-2, 5.1033608E-2, 3.1664744E-2, 
      7.0896583E-3, 1.1167203E-3, 1.3600086E-4, 4.7732716E-2, 4.7943623E-2, 
      4.8154529E-2, 4.8365436E-2, 4.8576343E-2, 4.8787250E-2, 4.8998156E-2, 
      4.9209063E-2, 4.9419970E-2, 4.9630877E-2, 4.9841783E-2, 3.7972086E-2, 
      3.7556330E-2, 3.8929808E-2, 3.9498415E-2, 4.1739003E-2, 4.3342075E-2, 
      4.5140816E-2, 4.7099160E-2, 4.8803325E-2, 6.8785386E-2, 7.9602002E-2, 
      8.0639220E-2, 6.7554384E-2, 5.1605256E-2, 1.8999509E-2, 4.8550702E-3, 
      9.4687407E-4, 2.6874094E-2, 2.6993039E-2, 2.7111984E-2, 2.7230929E-2, 
      2.7349874E-2, 2.7468819E-2, 2.7587764E-2, 2.7706709E-2, 2.7825655E-2, 
      2.7944600E-2, 2.8063545E-2, 2.3442661E-2, 2.3582995E-2, 2.4511956E-2, 
      2.5132254E-2, 2.6459615E-2, 2.7541130E-2, 2.8713907E-2, 2.9964567E-2, 
      3.1123242E-2, 4.4504338E-2, 6.0366777E-2, 7.0259410E-2, 6.9296765E-2, 
      6.1699292E-2, 3.2765153E-2, 1.1948925E-2, 3.2949346E-3, 1.4974727E-2, 
      1.5042154E-2, 1.5109581E-2, 1.5177009E-2, 1.5244436E-2, 1.5311863E-2, 
      1.5379290E-2, 1.5446718E-2, 1.5514145E-2, 1.5581572E-2, 1.5649000E-2, 
      1.4062135E-2, 1.4334238E-2, 1.4934864E-2, 1.5429147E-2, 1.6207524E-2, 
      1.6901398E-2, 1.7638605E-2, 1.8414586E-2, 1.9164032E-2, 2.7808184E-2, 
      4.2543911E-2, 5.5321907E-2, 6.1777780E-2, 6.2037752E-2, 4.4054173E-2, 
      2.1316072E-2, 7.7446092E-3, 8.4353212E-3, 8.4742359E-3, 8.5131506E-3, 
      8.5520653E-3, 8.5909800E-3, 8.6298948E-3, 8.6688095E-3, 8.7077242E-3, 
      8.7466389E-3, 8.7855536E-3, 8.8244683E-3, 8.3710316E-3, 8.6167465E-3, 
      8.9957769E-3, 9.3428054E-3, 9.8019165E-3, 1.0236265E-2, 1.0692408E-2, 
      1.1168927E-2, 1.1641770E-2, 1.7130898E-2, 2.8743291E-2, 4.0910896E-2, 
      5.0447835E-2, 5.5887987E-2, 5.0464590E-2, 3.0885530E-2, 1.4119899E-2, 
      4.8372362E-3, 4.8601058E-3, 4.8829754E-3, 4.9058449E-3, 4.9287145E-3, 
      4.9515841E-3, 4.9744537E-3, 4.9973233E-3, 5.0201929E-3, 5.0430625E-3, 
      5.0659321E-3, 4.9908487E-3, 5.1729762E-3, 5.4090161E-3, 5.6379362E-3, 
      5.9114564E-3, 6.1802168E-3, 6.4607323E-3, 6.7526287E-3, 7.0473395E-3, 
      1.0503270E-2, 1.8929055E-2, 2.9026523E-2, 3.8832222E-2, 4.6703863E-2, 
      5.1708108E-2, 3.8685646E-2, 2.1538167E-2, 2.8247806E-3, 2.8384183E-3, 
      2.8520559E-3, 2.8656935E-3, 2.8793312E-3, 2.8929688E-3, 2.9066064E-3, 
      2.9202441E-3, 2.9338817E-3, 2.9475193E-3, 2.9611570E-3, 2.9912324E-3, 
      3.1150193E-3, 3.2610145E-3, 3.4071855E-3, 3.5717432E-3, 3.7372755E-3, 
      3.9095484E-3, 4.0885092E-3, 4.2712017E-3, 6.4377777E-3, 1.2269204E-2, 
      2.0016752E-2, 2.8669314E-2, 3.6977980E-2, 4.8822771E-2, 4.3513146E-2, 
      2.8785467E-2, 1.6757407E-3, 1.6839622E-3, 1.6921837E-3, 1.7004052E-3, 
      1.7086267E-3, 1.7168482E-3, 1.7250698E-3, 1.7332913E-3, 1.7415128E-3, 
      1.7497343E-3, 1.7579558E-3, 1.8042052E-3, 1.8847145E-3, 1.9747558E-3, 
      2.0665242E-3, 2.1663515E-3, 2.2682125E-3, 2.3741105E-3, 2.4840772E-3, 
      2.5971324E-3, 3.9527222E-3, 7.8733232E-3, 1.3528940E-2, 2.0530845E-2, 
      2.8124822E-2, 4.3318170E-2, 4.5064473E-2, 3.4745272E-2, 1.0067245E-3, 
      1.0117210E-3, 1.0167174E-3, 1.0217139E-3, 1.0267104E-3, 1.0317068E-3, 
      1.0367033E-3, 1.0416997E-3, 1.0466962E-3, 1.0516927E-3, 1.0566891E-3, 
      1.0950643E-3, 1.1462190E-3, 1.2017244E-3, 1.2588651E-3, 1.3198301E-3, 
      1.3825751E-3, 1.4478029E-3, 1.5155576E-3, 1.5855351E-3, 2.4331310E-3, 
      5.0206427E-3, 9.0124483E-3, 1.4371127E-2, 2.0745127E-2, 3.6609014E-2, 
      4.3725435E-2, 3.8682725E-2 };

      alpha = { 1E-2, 2E-2, 4E-2, 6E-2, 8E-2, 2E-1, 4E-1, 6E-1, 8E-1, 2, 4, 6, 
      8, 10, 15, 20, 25, 30, 35, 40, 45, 50};

      beta = { 0, 1E-2, 2E-2, 3E-2, 4E-2, 5E-2, 6E-2, 7E-2, 8E-2, 9E-2, 0.1, 
      0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 4, 6, 8, 10, 15, 20, 25 };
      for ( auto& a : alpha ){ a *= scaling; }
      for ( auto& b : beta  ){ b *= sc;      }



      std::vector<double> sab( alpha.size()*beta.size(), 0.0 );

      output = continuous( nphon, delta, tbeta, rho, alpha, beta, sab );
      THEN( "continuous output matches expected value" ){
        lambda_s =   0.2352041964494244, t_bar = 1.9344942054735312;
        REQUIRE(ranges::equal(sabCorrect, sab, equal));
        REQUIRE(lambda_s == Approx(std::get<0>(output)).epsilon(1e-6));
        REQUIRE(t_bar    == Approx(std::get<1>(output)).epsilon(1e-6));
      } // THEN
    } // WHEN




  } // GIVEN 

  GIVEN( "Beryllium metal example" ){

    alpha = { 3.052968E-3, 4.022632E-3, 6.084044E-3, 8.016418E-3, 1.212446E-2, 
    1.597535E-2, 2.416197E-2, 3.183614E-2, 4.815068E-2, 6.344398E-2,9.595606E-2,
    1.264330E-1, 1.912240E-1, 2.519593E-1, 3.810768E-1, 5.021119E-1,7.594208E-1, 
    1.000623E+0, 1.513396E+0, 1.994071E+0, 3.015939E+0, 3.973842E+0, 6.010250E0, 
    7.919187E+0, 1.197740E+1, 1.578158E+1, 2.386891E+1, 3.145000E+1, 4.756666E1, 
    6.267447E+1, 6.714871E+1, 8.847604E+1 };

    beta = { 0.000000E+0, 1.091104E-1, 3.000536E-1, 3.818864E-1, 5.728296E-1, 
    6.546624E-1, 8.456056E-1, 9.274384E-1, 1.118382E+0, 1.200214E+0, 1.391158, 
    1.472990E+0, 1.663934E+0, 1.745766E+0, 1.936710E+0, 2.018542E+0, 2.209486, 
    2.291318E+0, 2.482262E+0, 2.564094E+0, 2.755038E+0, 2.836870E+0, 3.027814, 
    3.109647E+0, 3.300590E+0, 3.382423E+0, 4.536305E+0, 5.358183E+0, 7.902297, 
    9.334018E+0, 1.376589E+1, 1.625997E+1, 2.398035E+1, 2.832506E+1, 4.177405E1, 
    4.934258E+1, 7.277088E+1, 8.595535E+1, 1.267677E+2, 1.497352E+2, 1.673147E2, 
    1.976285E+2 }; 


    rho = { 0.0000E+0, 7.2477E-4, 3.7084E-3, 8.0087E-3, 1.0642E-2, 1.5897E-2, 
    2.7372E-2, 4.1843E-2, 5.0214E-2, 6.5036E-2, 8.3674E-2, 9.9329E-2, 1.1977E-1, 
    1.4296E-1, 1.6484E-1, 1.8945E-1, 2.1887E-1, 2.3537E-1, 2.6166E-1, 3.0003E-1, 
    3.4054E-1, 3.8728E-1, 4.2481E-1, 4.7598E-1, 5.1890E-1, 5.7400E-1, 6.2970E-1, 
    6.5754E-1, 7.2042E-1, 7.9118E-1, 8.6756E-1, 9.2948E-1, 1.0030E+0, 1.1163E+0, 
    1.2048E+0, 1.2870E+0, 1.4139E+0, 1.5249E+0, 1.6221E+0, 1.7638E+0, 1.8924E+0, 
    2.0388E+0, 2.2056E+0, 2.3709E+0, 2.5558E+0, 2.7595E+0, 3.0108E+0, 3.2603E+0, 
    3.5066E+0, 3.7442E+0, 4.0067E+0, 4.3677E+0, 4.7164E+0, 5.0820E+0, 5.5881E+0, 
    6.0898E+0, 6.5510E+0, 7.0877E+0, 7.5931E+0, 8.0736E+0, 8.6232E+0, 9.2283E+0, 
    9.9334E+0, 1.0613E+1, 1.1278E+1, 1.1973E+1, 1.2784E+1, 1.3744E+1, 1.4739E+1, 
    1.5918E+1, 1.7654E+1, 1.9834E+1, 2.1455E+1, 2.2574E+1, 2.3744E+1, 2.4900E+1, 
    2.6227E+1, 2.7931E+1, 2.9747E+1, 2.9884E+1, 2.7358E+1, 2.4817E+1, 2.3690E+1, 
    2.3242E+1, 2.3624E+1, 2.3473E+1, 2.2368E+1, 2.1447E+1, 2.0724E+1, 2.1121E+1, 
    2.4240E+1, 2.7607E+1, 2.7643E+1, 2.5431E+1, 2.3755E+1, 2.3377E+1, 2.3410E+1, 
    2.3504E+1, 2.3647E+1, 2.3681E+1, 2.3805E+1, 2.3714E+1, 2.3385E+1, 2.3050E+1, 
    2.2244E+1, 2.1008E+1, 1.9536E+1, 1.8341E+1, 1.8075E+1, 1.8606E+1, 1.9599E+1, 
    2.1037E+1, 2.3193E+1, 2.4016E+1, 2.3573E+1, 2.5664E+1, 3.0187E+1, 3.1256E+1, 
    2.7257E+1, 2.2765E+1, 1.4893E+1, 6.8192E+0, 3.8444E+0, 2.4718E+0, 1.3358E+0, 
    3.5968E-1, 0.0000E+0 };

std::vector<double> 
    everyFifthValue = { 7.528341E-5, 1.271158E-4, 3.154207E-4, 1.027164E-3, 
      7.179574E-4, 5.810458E-5, 3.58833E-18, 0.000000000, 0.000000000, 
      1.382893E-4, 2.721986E-4, 7.799903E-4, 1.137213E-3, 8.270636E-4, 
      3.00455E-10, 0.000000000, 0.000000000, 1.467789E-4, 2.950356E-4, 
      7.349432E-4, 2.082407E-3, 1.328274E-3, 2.029612E-6, 7.16969E-20, 
      0.000000000, 0.000000000, 3.119314E-4, 6.079728E-4, 2.293577E-3, 
      2.152916E-3, 9.246918E-4, 1.44555E-10, 0.000000000, 0.000000000, 
      3.915177E-4, 6.521952E-4, 2.043209E-3, 3.493972E-3, 2.115376E-3, 
      5.818367E-6, 8.41156E-27, 0.000000000, 3.937206E-4, 6.630886E-4, 
      1.640403E-3, 5.336090E-3, 3.731555E-3, 3.074530E-4, 1.40518E-14, 
      0.000000000, 0.000000000, 8.278304E-4, 1.622089E-3, 4.634833E-3, 
      6.756191E-3, 4.919080E-3, 6.453637E-8, 2.75069E-29, 0.000000000, 
      7.677845E-4, 1.531924E-3, 3.797075E-3, 1.074101E-2, 6.859254E-3, 
      5.481408E-5, 1.12817E-15, 0.000000000, 0.000000000, 1.857386E-3, 
      3.590771E-3, 1.347333E-2, 1.265474E-2, 5.471909E-3, 3.421227E-8, 
      0.000000000, 0.000000000, 2.032999E-3, 3.354916E-3, 1.040733E-2, 
      1.777452E-2, 1.080150E-2, 1.554518E-4, 6.37456E-21, 0.000000000, 
      2.354737E-3, 3.905568E-3, 9.486316E-3, 3.064815E-2, 2.149673E-2, 
      1.966517E-3, 1.08563E-10, 0.000000000, 0.000000000, 4.254641E-3, 
      8.154179E-3, 2.296504E-2, 3.344275E-2, 2.448795E-2, 8.840780E-6, 
      4.01926E-22, 0.000000000, 4.584144E-3, 8.766216E-3, 2.108500E-2, 
      5.903327E-2, 3.797827E-2, 1.816492E-3, 4.82295E-11, 0.000000000, 
      0.000000000, 9.279147E-3, 1.725367E-2, 6.299824E-2, 5.936082E-2, 
      2.654920E-2, 6.726930E-6, 8.51879E-30, 0.000000000, 1.155290E-2, 
      1.810867E-2, 5.297783E-2, 8.976278E-2, 5.583998E-2, 4.821793E-3, 
      2.64788E-14, 0.000000000, 1.181477E-2, 1.838391E-2, 4.102864E-2, 
      1.279588E-1, 9.120142E-2, 1.263264E-2, 3.965660E-7, 0.000000000, 
      0.000000000, 2.231310E-2, 3.824733E-2, 9.889100E-2, 1.431105E-1, 
      1.085762E-1, 1.452652E-3, 3.13758E-14, 0.000000000, 2.159523E-2, 
      3.546525E-2, 7.470613E-2, 1.980557E-1, 1.326887E-1, 3.248983E-2, 
      8.047784E-7, 0.000000000, 0.000000000, 3.947250E-2, 6.178125E-2, 
      1.916340E-1, 1.849184E-1, 1.006280E-1, 2.408368E-3, 3.99776E-18, 
      0.000000000, 3.997470E-2, 5.397322E-2, 1.237828E-1, 2.008207E-1, 
      1.415231E-1, 6.259761E-2, 3.884935E-8, 0.000000000, 3.772447E-2, 
      4.998176E-2, 8.148058E-2, 2.020080E-1, 1.629478E-1, 7.241268E-2, 
      1.573214E-3, 4.64100E-17, 0.000000000, 4.309686E-2, 5.873885E-2, 
      1.088398E-1, 1.526378E-1, 1.381639E-1, 5.266990E-2, 4.137640E-7, 
      0.000000000, 2.795324E-2, 3.699844E-2, 5.317895E-2, 9.957868E-2, 
      9.008784E-2, 1.039835E-1, 7.116744E-3, 7.06706E-21, 0.000000000, 
      2.271228E-2, 2.917271E-2, 5.064990E-2, 5.758423E-2, 5.868845E-2, 
      7.424323E-2, 7.353243E-8, 0.000000000, 7.652941E-3, 9.903402E-3, 
      1.409992E-2, 1.885109E-2, 2.276816E-2, 4.160121E-2, 9.569618E-3, 
      1.06221E-17, 2.548565E-3, 3.438908E-3, 4.770335E-3, 6.495148E-3, 
      8.434437E-3, 1.029217E-2, 5.891244E-2, 9.149304E-6, 0.000000000, 
      4.145125E-4, 5.865565E-4, 7.817987E-4, 1.078613E-3, 1.397314E-3, 
      7.203282E-3, 3.912367E-2, 7.62975E-20, 5.896026E-5, 8.428778E-5, 
      1.132740E-4, 1.586193E-4, 2.092513E-4, 4.254586E-4, 1.234288E-2, 
      2.423062E-5, 0.000000000, 1.705554E-6, 2.310225E-6, 3.276033E-6, 
      4.385315E-6, 6.132996E-6, 7.253713E-5, 3.137750E-2, 2.05286E-13, 
      4.616867E-8, 6.276863E-8, 8.955085E-8, 1.206710E-7, 1.703729E-7, 
      4.890455E-7, 3.219968E-4, 6.789171E-3, 1.431288E-8, 1.976221E-8, 
      2.826851E-8, 3.818909E-8, 5.409857E-8, 7.248202E-8, 5.326591E-6, 
      1.343825E-2, 1.11808E-10, 1.35589E-10, 1.94586E-10, 2.63874E-10, 
      3.75923E-10, 5.06612E-10, 4.504120E-9, 8.101733E-6, 1.689654E-3
};

    nphon = 100; 
    tbeta = 1.0; tev = 2.5507297688E-2; 
    delta = 0.00069552/tev; 
    sc = 0.0253/tev; scaling = sc;
    for ( auto& a : alpha ){ a *= scaling; }
    for ( auto& b : beta  ){ b *= sc;      }



      std::vector<double> sab( alpha.size()*beta.size(), 0.0 );

      output = continuous( nphon, delta, tbeta, rho, alpha, beta, sab );


    THEN( "continuous output matches expected value" ){
      lambda_s = 0.60333890798030132; t_bar = 1.4641327610045163;
      for (size_t i = 0; i < everyFifthValue.size(); ++i ){ 
        REQUIRE( sab[5*i] == Approx(everyFifthValue[i]).epsilon(1e-6) );
      }
      REQUIRE(lambda_s == Approx(std::get<0>(output)).epsilon(1e-6));
      REQUIRE(t_bar    == Approx(std::get<1>(output)).epsilon(1e-6));
    } // THEN


  } // GIVEN 


  GIVEN( "arbitrary phonon distribution" ){
    alpha = {0.01, 0.1,  1,   5, 10};
    beta  = {0.0,  0.05, 0.5, 5, 10, 20};
    std::vector<double> sab( alpha.size()*beta.size(), 0.0 );
    rho = {0, 0.0005, 0.001, 0.002, 0.0035, 0.005, 0.0075, 0.01, 0.013, 0.0165, 
    0.02, 0.0245, 0.029, 0.034, 0.0395, 0.045, 0.0506, 0.0562, 0.0622, 0.0686, 
    0.075, 0.083, 0.091, 0.09126, 0.0871, 0.0839, 0.0807, 0.07798, 0.07574, 
    0.0735, 0.07162, 0.06974, 0.0 };
    tbeta = 0.8;

    WHEN( "large number of phonon expansion terms" ){
      nphon = 300; 
      AND_WHEN( "temperature is very low" ){
        tev = 100*8.6173303e-5;   
        sc = 0.0253/tev; scaling = sc;
        delta = 0.005/tev;
        for ( auto& a : alpha ){ a *= scaling; }
        for ( auto& b : beta  ){ b *= sc;      }


 
        output = continuous( nphon, delta, tbeta, rho, alpha, beta, sab );
        THEN( "continuous output matches expected value" ){
          lambda_s = 7.6444701E-2; t_bar = 6.2888965;
          sabCorrect = {4.1395658E-5, 4.4726194E-5, 3.7025816E-5, 1.5765995E-4, 
          8.4143007E-8, 9.152651E-15, 4.0587963E-4, 4.3853189E-4, 3.6318734E-4, 
          1.5528410E-3, 8.3197978E-6, 9.109286E-11, 3.3329337E-3, 3.6007961E-3, 
          2.9949842E-3, 1.3332332E-2, 7.4136325E-4, 8.6242441E-7, 6.9436755E-3, 
          7.4993342E-3, 6.3573122E-3, 3.3430808E-2, 1.0643445E-2, 3.7753740E-4, 
          4.6516788E-3, 5.0220262E-3, 4.3586613E-3, 2.7573268E-2, 1.9902691E-2, 
          3.3329649E-3 };
          REQUIRE(ranges::equal(sabCorrect, sab, equal));
          REQUIRE(lambda_s == Approx(std::get<0>(output)).epsilon(1e-6));
          REQUIRE(t_bar    == Approx(std::get<1>(output)).epsilon(1e-6));
        } // THEN
      } // AND WHEN
      AND_WHEN( "temperature is very high" ){
        tev = 900*8.6173303e-5;   
        sc = 0.0253/tev; scaling = sc;
        delta = 0.005/tev;
        for ( auto& a : alpha ){ a *= scaling; }
        for ( auto& b : beta  ){ b *= sc;      }


        output = continuous( nphon, delta, tbeta, rho, alpha, beta, sab );
        THEN( "continuous output matches expected value" ){
          lambda_s =  1.729454237; t_bar = 1.16861142;
          sabCorrect = {3.3452442E-3, 3.3727931E-3, 1.7063857E-3, 1.7623769E-3, 
          1.1883155E-6, 2.080211E-13, 3.2121369E-2, 3.2383671E-2, 1.6544030E-2, 
          1.7133899E-2, 1.1666113E-4, 2.0741014E-9, 2.1565380E-1, 2.1728142E-1, 
          1.2210329E-1, 1.3033373E-1, 9.6744738E-3, 1.9554632E-5, 2.2530912E-1, 
          2.2675945E-1, 1.8153698E-1, 2.2861657E-1, 1.0322905E-1, 6.8448510E-3, 
          1.0895402E-1, 1.0974843E-1, 1.0983673E-1, 1.6582818E-1, 1.4975245E-1, 
          4.2614699E-2};
          REQUIRE(ranges::equal(sabCorrect, sab, equal));
          REQUIRE(lambda_s == Approx(std::get<0>(output)).epsilon(1e-6));
          REQUIRE(t_bar    == Approx(std::get<1>(output)).epsilon(1e-6));
        } // THEN
      } // AND WHEN
    } // WHEN
  } // GIVEN 
} // TEST CASE 




