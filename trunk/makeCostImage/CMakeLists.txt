PROJECT(makeCostImage)

INCLUDE_DIRECTORIES( ${BIGRCarotid_SOURCE_DIR}/medialness
                     ${BIGRCarotid_SOURCE_DIR}/lumenIntensitySimilarity
                     ${BIGRCarotid_SOURCE_DIR}/invertImage
                     ${BIGRCarotid_SOURCE_DIR}/utilities
                     )

LINK_DIRECTORIES(${BIGRCarotid_BINARY_DIR}/utilities )

ADD_EXECUTABLE(makeCostImage makeCostImage.cxx)

TARGET_LINK_LIBRARIES(makeCostImage utilities ITKStatistics)
