
add_executable( simple.discreteOsc.discreteOscTools.oscillatorLoopTools.test bfact.test.cpp )
target_compile_options( simple.discreteOsc.discreteOscTools.oscillatorLoopTools.test PRIVATE ${${PREFIX}_common_flags}
$<$<BOOL:${strict}>:${${PREFIX}_strict_flags}>$<$<CONFIG:DEBUG>:
${${PREFIX}_DEBUG_flags}
$<$<BOOL:${coverage}>:${${PREFIX}_coverage_flags}>>
$<$<CONFIG:RELEASE>:
${${PREFIX}_RELEASE_flags}
$<$<BOOL:${link_time_optimization}>:${${PREFIX}_link_time_optimization_flags}>
$<$<BOOL:${nonportable_optimization}>:${${PREFIX}_nonportable_optimization_flags}>>

${CXX_appended_flags} ${leapr_appended_flags} )
target_link_libraries( simple.discreteOsc.discreteOscTools.oscillatorLoopTools.test PUBLIC leapr ) 
add_test( NAME simple.discreteOsc.discreteOscTools.oscillatorLoopTools COMMAND simple.discreteOsc.discreteOscTools.oscillatorLoopTools.test )