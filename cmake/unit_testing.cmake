#######################################################################
# Setup
#######################################################################

message( STATUS "Adding leapr unit testing" )
enable_testing()


#######################################################################
# Unit testing directories
#######################################################################

#add_subdirectory( src/coherentElastic/coherentElastic_util/hexLatticeFactors_util/test )
add_subdirectory( src/coherentElastic/coherentElastic_util/test )
add_subdirectory( src/coherentElastic/test )
#add_subdirectory( src/coldHydrogen/coldHydrogen_util/betaLoop_util/jprimeLoop_util/sumh_util/test )
add_subdirectory( src/coldHydrogen/coldHydrogen_util/betaLoop_util/jprimeLoop_util/test )
add_subdirectory( src/coldHydrogen/coldHydrogen_util/betaLoop_util/test )
add_subdirectory( src/coldHydrogen/coldHydrogen_util/test )
add_subdirectory( src/coldHydrogen/test )
#add_subdirectory( src/continuous/continuous_util/start_util/test )
add_subdirectory( src/continuous/continuous_util/test )
add_subdirectory( src/continuous/test )
add_subdirectory( src/discreteOscillators/discreteOscillators_util/oscLoopFuncs_util/test )
add_subdirectory( src/discreteOscillators/discreteOscillators_util/test )
add_subdirectory( src/discreteOscillators/test )
add_subdirectory( src/endout/endout_util/test )
add_subdirectory( src/endout/test )
add_subdirectory( src/generalTools/test )
add_subdirectory( src/skold/test )
add_subdirectory( src/test )
add_subdirectory( src/translational/test )
add_subdirectory( src/translational/translational_util/test )
