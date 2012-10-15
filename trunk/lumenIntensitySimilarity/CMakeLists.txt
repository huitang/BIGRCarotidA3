PROJECT(lumenIntensitySimilarity)

INCLUDE_DIRECTORIES( ${BIGRCarotid_SOURCE_DIR}/utilities )

LINK_DIRECTORIES(${BIGRCarotid_BINARY_DIR}/utilities )


ADD_EXECUTABLE(lumenIntensitySimilarity 
  lumenIntensitySimilarity.cxx 
)

TARGET_LINK_LIBRARIES(lumenIntensitySimilarity utilities ITKStatistics)


