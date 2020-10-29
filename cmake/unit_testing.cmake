#######################################################################
# Setup
#######################################################################

message( STATUS "Adding leapr unit testing" )
enable_testing()


#######################################################################
# Unit testing directories
#######################################################################

add_subdirectory( src/coher/coher_util/hexLatticeFactors_util/test )
add_subdirectory( src/coher/coher_util/test )
add_subdirectory( src/coher/test )
add_subdirectory( src/coldh/coldh_util/betaLoop_util/jprimeLoop_util/sumh_util/test )
add_subdirectory( src/coldh/coldh_util/betaLoop_util/jprimeLoop_util/test )
add_subdirectory( src/coldh/coldh_util/betaLoop_util/test )
add_subdirectory( src/coldh/coldh_util/test )
add_subdirectory( src/coldh/test )
#add_subdirectory( src/contin/contin_util/start_util/test )
add_subdirectory( src/contin/contin_util/test )
add_subdirectory( src/contin/test )
add_subdirectory( src/discre/discre_util/oscLoopFuncs_util/test )
add_subdirectory( src/discre/discre_util/test )
add_subdirectory( src/discre/test )
add_subdirectory( src/endout/endout_util/test )
add_subdirectory( src/endout/test )
add_subdirectory( src/generalTools/test )
add_subdirectory( src/skold/test )
add_subdirectory( src/test )
add_subdirectory( src/trans/test )
add_subdirectory( src/trans/trans_util/test )
