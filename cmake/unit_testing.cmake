#######################################################################
# Setup
#######################################################################

message( STATUS "Adding leapr unit testing" )
enable_testing()


#######################################################################
# Unit testing directories
#######################################################################

<<<<<<< HEAD
add_subdirectory( src/coher/coher_util/hexLatticeFactors_util/test )
add_subdirectory( src/coher/coher_util/test )
add_subdirectory( src/coher/test )
add_subdirectory( src/coldh/coldh_util/betaLoop_util/jprimeLoop_util/sumh_util/test )
add_subdirectory( src/coldh/coldh_util/betaLoop_util/jprimeLoop_util/test )
add_subdirectory( src/coldh/coldh_util/betaLoop_util/test )
add_subdirectory( src/coldh/coldh_util/test )
add_subdirectory( src/coldh/test )
#add_subdirectory( src/contin/contin_util/start_util/test )
add_subdirectory( src/continuous/continuous_util/test )
add_subdirectory( src/continuous/test )
add_subdirectory( src/discre/discre_util/oscLoopFuncs_util/test )
add_subdirectory( src/discre/discre_util/test )
add_subdirectory( src/discre/test )
=======
add_subdirectory( src/coherentElastic/coherentElastic_util/hexLatticeFactors_util/test )
add_subdirectory( src/coherentElastic/coherentElastic_util/test )
add_subdirectory( src/coherentElastic/test )
add_subdirectory( src/coldHydrogen/coldHydrogen_util/betaLoop_util/jprimeLoop_util/test )
add_subdirectory( src/coldHydrogen/coldHydrogen_util/betaLoop_util/test )
add_subdirectory( src/coldHydrogen/coldHydrogen_util/test )
add_subdirectory( src/coldHydrogen/test )
add_subdirectory( src/continuous/continuous_util/test )
add_subdirectory( src/continuous/test )
add_subdirectory( src/discreteOscillators/discreteOscillators_util/oscLoopFuncs_util/test )
add_subdirectory( src/discreteOscillators/discreteOscillators_util/test )
add_subdirectory( src/discreteOscillators/test )
>>>>>>> 8fb05255c06fae5c8be6f84083f1c7cfc26f4003
add_subdirectory( src/endout/endout_util/test )
add_subdirectory( src/endout/test )
add_subdirectory( src/generalTools/test )
add_subdirectory( src/skold/test )
add_subdirectory( src/test )
add_subdirectory( src/translational/test )
add_subdirectory( src/translational/translational_util/test )
