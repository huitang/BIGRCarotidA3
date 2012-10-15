PROJECT(medialness)
										
LINK_LIBRARIES(
  ITKAlgorithms
  ITKCommon
  ITKBasicFilters
  ITKIO
  ITKStatistics
)

ADD_EXECUTABLE(medialness medialness.cxx)
