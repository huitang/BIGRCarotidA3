PROJECT (TPGAC)
										
LINK_LIBRARIES(
  ITKAlgorithms
  ITKCommon
  ITKBasicFilters
  ITKIO
)

ADD_EXECUTABLE(TPGAC TPGAC.cxx)

