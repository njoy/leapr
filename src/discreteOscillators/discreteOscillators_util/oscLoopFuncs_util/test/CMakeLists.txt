
add_executable( discreteOscillators.discreteOscillators_util.oscLoopFuncs_util.test bfact.test.cpp )
target_compile_options( discreteOscillators.discreteOscillators_util.oscLoopFuncs_util.test PRIVATE ${${PREFIX}_common_flags}
$<$<BOOL:${strict}>:${${PREFIX}_strict_flags}>$<$<CONFIG:DEBUG>:
${${PREFIX}_DEBUG_flags}
$<$<BOOL:${coverage}>:${${PREFIX}_coverage_flags}>>
$<$<CONFIG:RELEASE>:
${${PREFIX}_RELEASE_flags}
$<$<BOOL:${link_time_optimization}>:${${PREFIX}_link_time_optimization_flags}>
$<$<BOOL:${nonportable_optimization}>:${${PREFIX}_nonportable_optimization_flags}>>

${CXX_appended_flags} ${leapr_appended_flags} )
target_link_libraries( discreteOscillators.discreteOscillators_util.oscLoopFuncs_util.test PUBLIC leapr ) 
add_test( NAME discreteOscillators.discreteOscillators_util.oscLoopFuncs_util COMMAND discreteOscillators.discreteOscillators_util.oscLoopFuncs_util.test )
