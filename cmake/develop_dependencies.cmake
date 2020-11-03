cmake_minimum_required( VERSION 3.14 )
include( FetchContent )

#######################################################################
# Declare project dependencies
#######################################################################

FetchContent_Declare( ENDFtk
    GIT_REPOSITORY  https://github.com/njoy/ENDFtk
    GIT_TAG         origin/build/fetchcontent
    GIT_SHALLOW     TRUE
    )

FetchContent_Declare( catch-adapter
    GIT_REPOSITORY  https://github.com/njoy/catch-adapter
    GIT_TAG         origin/build/fetchcontent
    GIT_SHALLOW     TRUE
    )

FetchContent_Declare( range-v3-adapter
    GIT_REPOSITORY  https://github.com/njoy/range-v3-adapter
    GIT_TAG         origin/build/fetchcontent
    GIT_SHALLOW     TRUE
    )

FetchContent_Declare( lipservice 
    GIT_REPOSITORY  https://github.com/njoy/lipservice
    GIT_TAG         origin/build/fetchcontent
    GIT_SHALLOW     TRUE
    )



#######################################################################
# Load dependencies
#######################################################################

FetchContent_MakeAvailable(
    ENDFtk
    catch-adapter
    range-v3-adapter
    lipservice 
    )
