CMAKE_MINIMUM_REQUIRED(VERSION 2.4)
PROJECT(BIGRCarotid)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
  cmake_policy(SET CMP0015 NEW)
endif(COMMAND cmake_policy)

SET(Boost_USE_STATIC_LIBS OFF)
#SET(Boost_USE_MULTITHREADED ON)
#FIND_PACKAGE(Boost REQUIRED COMPONENTS program_options filesystem system serialization regex)



IF(MSVC)
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
ENDIF(MSVC)

# MeVisLab
SET( MLAB_ROOT $ENV{MLAB_ROOT} CACHE PATH "Path to MLAB_ROOT")

#ADD_DEFINITIONS( -DBOOST_ALL_DYN_LINK )

INCLUDE_DIRECTORIES( 
  ${MLAB_ROOT}/MeVis/ThirdParty/Sources/misc/boost
) 
LINK_DIRECTORIES(
  ${MLAB_ROOT}/MeVis/ThirdParty/lib 
)

IF(UNIX) 
#  LINK_LIBRARIES( boost_program_options-gcc41-mt-1_36 ) 
  FIND_PACKAGE(Boost 1.33.0 REQUIRED) 
  IF(Boost_FOUND) 
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS}) 
    LINK_LIBRARIES(${Boost_LIBRARIES}) 
  ELSE(Boost_FOUND) 
    MESSAGE( FATAL_ERROR "Cannot build TTProgram without BOOST. Please set BOOST_DIR.")  
  ENDIF(Boost_FOUND) 
ENDIF(UNIX) 
IF(UNIX) 
  LINK_LIBRARIES( libboost_program_options.so libboost_filesystem.so libboost_date_time.so ) 
ENDIF(UNIX) 


FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "Cannot build BIGRCarotid without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)




SUBDIRS(
  utilities
  medialness
  invertImage
  multiplyImage
  #pathtrackingCMPR
  lumenIntensitySimilarity
  makeCostImage
  pathtracking
  makeMPRStack
  #medialnessNew
)

