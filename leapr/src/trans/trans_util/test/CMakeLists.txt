
add_executable( trans.trans_util.test
                terps.test.cpp
                sbfill.test.cpp
                s_table_generation.test.cpp
                bessel_K1.test.cpp )
target_compile_options( trans.trans_util.test PRIVATE ${${PREFIX}_common_flags}
$<$<BOOL:${strict}>:${${PREFIX}_strict_flags}>$<$<CONFIG:DEBUG>:
${${PREFIX}_DEBUG_flags}
$<$<BOOL:${coverage}>:${${PREFIX}_coverage_flags}>>
$<$<CONFIG:RELEASE>:
${${PREFIX}_RELEASE_flags}
$<$<BOOL:${link_time_optimization}>:${${PREFIX}_link_time_optimization_flags}>
$<$<BOOL:${nonportable_optimization}>:${${PREFIX}_nonportable_optimization_flags}>>

${CXX_appended_flags} ${leapr_appended_flags} )
target_link_libraries( trans.trans_util.test PUBLIC leapr ) 
add_test( NAME trans.trans_util COMMAND trans.trans_util.test )