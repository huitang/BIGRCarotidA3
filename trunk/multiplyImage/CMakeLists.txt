PROJECT (multiplyImage)
			
LINK_LIBRARIES(
  ITKCommon
  ITKBasicFilters
  ITKIO
)

ADD_EXECUTABLE(multiplyImage multiplyImage.cxx)
