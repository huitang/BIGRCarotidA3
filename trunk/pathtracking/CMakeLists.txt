PROJECT(pathtracking)

LINK_LIBRARIES(
  ITKCommon
  ITKIO
)


INCLUDE_DIRECTORIES( ${BIGRCarotid_SOURCE_DIR}/utilities )

LINK_DIRECTORIES( ${BIGRCarotid_BINARY_DIR}/utilities )

ADD_EXECUTABLE(pathtracking
               pathTracking.cxx )
     
TARGET_LINK_LIBRARIES(pathtracking utilities)

