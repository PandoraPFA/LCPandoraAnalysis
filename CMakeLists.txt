########################################################
# cmake file for building PandoraAnalysis
# @author Jan Engels, DESY
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################


# project name
PROJECT( PandoraAnalysis )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 7 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )



### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


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

### DOCUMENTATION ############################################################
OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ELSE()
    MESSAGE( STATUS "INSTALL_DOC set to OFF" )
    SET( INSTALL_DOC OFF )
ENDIF()

### CALIBRATION DOCUMENTATION ###############################################

#include directories
SET( LATEX_OUTPUT_PATH build/doc )

INCLUDE( "${PROJECT_SOURCE_DIR}/cmakemodules/UseLATEX.cmake" )

ADD_LATEX_DOCUMENT(doc/CalibrationPandoraAnalysisExplained.tex
DEFAULT_PDF
)

ADD_CUSTOM_COMMAND(TARGET pdf POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy
                        ${PROJECT_SOURCE_DIR}/build/doc/CalibrationPandoraAnalysisExplained.pdf
                        ${PROJECT_SOURCE_DIR}/doc/CalibrationPandoraAnalysisExplained.pdf)

### LIBRARY AND TOOLS #######################################################

#include directories
INCLUDE_DIRECTORIES( BEFORE ./include ) # use BEFORE to fix problem on mac osx (case insensitive filesystem)
#INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
ADD_DEFINITIONS( "-Wno-long-long" )

AUX_SOURCE_DIRECTORY( ./src library_sources )
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION lib )


# make standalone analysis executables
SET( PROGRAMS1 AnalysePerformance AnalysePerformanceFull )

FOREACH( program ${PROGRAMS1} )
    ADD_EXECUTABLE( ${program} ${PROJECT_SOURCE_DIR}/performance/${program}.cc )
    TARGET_LINK_LIBRARIES( ${program} ${PROJECT_NAME} )

    # install executable
    INSTALL( TARGETS ${program} DESTINATION bin PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE )
ENDFOREACH( program )

# make standalone analysis executables
SET( PROGRAMS2 ECalDigitisation_ContainedEvents HCalDigitisation_ContainedEvents HCalDigitisation_DirectionCorrectionDistribution HCalDigitisation_Ring PandoraPFACalibrate_MipResponse PandoraPFACalibrate_EMScale PandoraPFACalibrate_HadronicScale_ChiSquareMethod PandoraPFACalibrate_HadronicScale_TotalEnergyMethod PandoraPFACalibrate_HadronicEnergyGaussianFit RescaleEnergies )

FOREACH( program ${PROGRAMS2} )
    ADD_EXECUTABLE( ${program} ${PROJECT_SOURCE_DIR}/calibration/${program}.cc )
    TARGET_LINK_LIBRARIES( ${program} ${PROJECT_NAME} )

    # install executable
    INSTALL( TARGETS ${program} DESTINATION bin PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE )
ENDFOREACH( program )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()
