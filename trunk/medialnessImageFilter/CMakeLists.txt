PROJECT (medialnessImageFilter)
cmake_minimum_required(VERSION 2.4)
IF(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "Cannot build InsightApplications without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

FIND_PACKAGE(Boost)
IF(Boost_FOUND)
  INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS})
	LINK_DIRECTORIES( ${Boost_LIBRARY_DIRS})
ELSE(Boost_FOUND)
  MESSAGE(FATAL_ERROR
          "Cannot build without boost. Please set BOOST_DIR.")
ENDIF(Boost_FOUND)

IF (WIN32)
  INCLUDE_DIRECTORIES("E:\\packages\\boost\\boost")
  LINK_DIRECTORIES("C:\\Packages\\lib64")
ENDIF (WIN32)
										
LINK_LIBRARIES(
  ITKAlgorithms
  ITKCommon
  ITKBasicFilters
  ITKIO
  ITKStatistics
)

ADD_EXECUTABLE(medialnessImageFilter medialnessImageFilter.cxx)
#INSTALL(TARGETS medialnessImageFilter DESTINATION . )

