
add_executable( contin.contin_util.start_util.test
                getDebyeWaller.test.cpp
                normalize.test.cpp
                getEffectiveTemp.test.cpp )
target_compile_options( contin.contin_util.start_util.test PRIVATE ${${PREFIX}_common_flags}
$<$<BOOL:${strict}>:${${PREFIX}_strict_flags}>$<$<CONFIG:DEBUG>:
${${PREFIX}_DEBUG_flags}
$<$<BOOL:${coverage}>:${${PREFIX}_coverage_flags}>>
$<$<CONFIG:RELEASE>:
${${PREFIX}_RELEASE_flags}
$<$<BOOL:${link_time_optimization}>:${${PREFIX}_link_time_optimization_flags}>
$<$<BOOL:${nonportable_optimization}>:${${PREFIX}_nonportable_optimization_flags}>>

${CXX_appended_flags} ${leapr_appended_flags} )
target_link_libraries( contin.contin_util.start_util.test PUBLIC leapr ) 
add_test( NAME contin.contin_util.start_util COMMAND contin.contin_util.start_util.test )