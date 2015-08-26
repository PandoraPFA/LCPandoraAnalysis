########################################################
# cmake file for building PandoraAnalysis
# @author Jan Engels, DESY
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.2 FATAL_ERROR)
########################################################


# project name
PROJECT( PandoraAnalysis )


# project version
SET( PandoraAnalysis_VERSION_MAJOR 01 )
SET( PandoraAnalysis_VERSION_MINOR 00 )
SET( PandoraAnalysis_VERSION_PATCH 01 )


### CMAKE ###################################################################

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


### SOURCE ##################################################################

# include directories
INCLUDE_DIRECTORIES( BEFORE ./include ) # use BEFORE to fix problem on mac osx (case insensitive filesystem)

AUX_SOURCE_DIRECTORY( ./src PandoraAnalysis_SRCS )

# analysis executables
SET( PERFORMANCE_PROGRAMS AnalysePerformance
    AnalysePerformanceFull )
SET( CALIBRATION_PROGRAMS ECalDigitisation_ContainedEvents
    HCalDigitisation_ContainedEvents
    HCalDigitisation_DirectionCorrectionDistribution
    PandoraPFACalibrate_MipResponse
    PandoraPFACalibrate_EMScale
    PandoraPFACalibrate_HadronicScale_ChiSquareMethod
    PandoraPFACalibrate_HadronicScale_TotalEnergyMethod
    PandoraPFACalibrate_HadronicEnergyGaussianFit
    RescaleEnergies 
    SimCaloHitEnergyDistribution )
SET( PROGRAMS ${PERFORMANCE_PROGRAMS} ${CALIBRATION_PROGRAMS} )

# require proper C++
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )

# need long long for int64 for now
ADD_DEFINITIONS( "-Wno-long-long -Wno-sign-compare -fno-strict-aliasing" )


### DEPENDENCIES ############################################################

FIND_PACKAGE( Marlin 1.0 REQUIRED )
FIND_PACKAGE( MarlinUtil REQUIRED )
FIND_PACKAGE( LCIO REQUIRED )
FIND_PACKAGE( ROOT REQUIRED )

FOREACH( pkg Marlin MarlinUtil LCIO ROOT )
    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
    ENDIF()
ENDFOREACH()


### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ELSE()
    MESSAGE( STATUS "INSTALL_DOC set to OFF" )
    SET( INSTALL_DOC OFF )
ENDIF()

INCLUDE( "${PROJECT_SOURCE_DIR}/cmakemodules/UseLATEX.cmake" )
SET( LATEX_OUTPUT_PATH build/doc )
ADD_LATEX_DOCUMENT( doc/CalibrationPandoraAnalysisExplained.tex DEFAULT_PDF )


### OUTPUT ##################################################################

# add library
ADD_SHARED_LIBRARY( PandoraAnalysis ${PandoraAnalysis_SRCS} )

FOREACH( program ${PERFORMANCE_PROGRAMS} )
    ADD_EXECUTABLE( ${program} ${PROJECT_SOURCE_DIR}/performance/${program}.cc )
ENDFOREACH( program )

FOREACH( program ${CALIBRATION_PROGRAMS} )
    ADD_EXECUTABLE( ${program} ${PROJECT_SOURCE_DIR}/calibration/${program}.cc )
ENDFOREACH( program )

FOREACH( program ${PROGRAMS} )
    TARGET_LINK_LIBRARIES( ${program} PandoraAnalysis )
ENDFOREACH( program )


### INSTALL #################################################################

# install library
INSTALL_SHARED_LIBRARY( PandoraAnalysis DESTINATION lib )

# install executables
FOREACH( program ${PROGRAMS} )
    INSTALL( TARGETS ${program} DESTINATION bin PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE )
ENDFOREACH( program )

# install calibration documentation
ADD_CUSTOM_COMMAND( TARGET pdf POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${PROJECT_SOURCE_DIR}/build/doc/CalibrationPandoraAnalysisExplained.pdf
    ${PROJECT_SOURCE_DIR}/doc/CalibrationPandoraAnalysisExplained.pdf )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()
