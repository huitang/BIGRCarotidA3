PROJECT(utilities)

ADD_LIBRARY( utilities STATIC 
  SeedPointFileIO.cxx
  pathTrackingLibrary.h
  pathTrackingLibrary.cxx
  path.h
  path.cxx
  splinefitter.h
  splinefitter.cpp  
)

TARGET_LINK_LIBRARIES(utilities ITKCommon ITKIO)
