cmake_minimum_required( VERSION 3.2 ) 

project( catch-adapter LANGUAGES CXX )
get_directory_property( is_subproject PARENT_DIRECTORY )
include( CMakeDependentOption REQUIRED )

 if ( NOT GIT_EXECUTABLE )
     find_package( Git )
     if ( NOT GIT_FOUND )
         message( FATAL_ERROR "git installation was not found." )
     endif()
 endif()

 execute_process(
     COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
     WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
     OUTPUT_VARIABLE GIT_BRANCH
     OUTPUT_STRIP_TRAILING_WHITESPACE
)
 execute_process(
     COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
     WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
     OUTPUT_VARIABLE GIT_HASH
     OUTPUT_STRIP_TRAILING_WHITESPACE
) 

message( STATUS "" )
message( STATUS "-----------------------------------------------------------" )
message( STATUS "" )
message( STATUS "catch-adapter" )
message( STATUS "Git current branch: ${GIT_BRANCH}" )
message( STATUS "Git commit hash: ${GIT_HASH}" )
message( STATUS "" )
message( STATUS "-----------------------------------------------------------" )

add_library( catch-adapter INTERFACE )
target_include_directories( catch-adapter INTERFACE src/single_include )
target_sources( catch-adapter INTERFACE
                "${CMAKE_CURRENT_SOURCE_DIR}/src/single_include/catch.hpp" )
target_include_directories( catch-adapter INTERFACE 
                            "${CMAKE_CURRENT_SOURCE_DIR}/src/single_include" )
