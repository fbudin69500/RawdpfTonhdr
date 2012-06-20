cmake_minimum_required(VERSION 2.6)

project( RawdpfTonhdr )


  FIND_PACKAGE(ITK REQUIRED)
  IF(USE_ITK_FILE)
    INCLUDE(${USE_ITK_FILE})
  ENDIF(USE_ITK_FILE)

  SET(ITK_LIBRARIES
    ITKCommon
    ITKIO
    ITKNumerics
    ITKBasicFilters
    ITKAlgorithms
    )

find_package(GenerateCLP REQUIRED)
if(GenerateCLP_FOUND)
  include(${GenerateCLP_USE_FILE})
else(GenerateCLP_FOUND)
   message(FATAL_ERROR "GenerateCLP Not Found")
endif(GenerateCLP_FOUND)

set( RawdpfTonhdr_SOURCE RawdpfTonhdr.cpp)
GENERATECLP( RawdpfTonhdr_SOURCE RawdpfTonhdr.xml)


add_executable( RawdpfTonhdr ${RawdpfTonhdr_SOURCE} )
target_link_libraries( RawdpfTonhdr  ${ITK_LIBRARIES})
