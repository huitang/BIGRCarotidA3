PROJECT (invertImage)
										
LINK_LIBRARIES(
  ITKAlgorithms
  ITKCommon
  ITKBasicFilters
  ITKIO
)

ADD_EXECUTABLE(invertImage invertImageFilter.cxx)

