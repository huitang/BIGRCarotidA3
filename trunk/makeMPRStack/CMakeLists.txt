PROJECT(makeMPRStack)

SET( EMCPACKAGE $ENV{MLAB_ROOT} CACHE PATH "Path to EMC Package")

INCLUDE_DIRECTORIES( 
  ${BIGRCarotid_SOURCE_DIR}/utilities
  ${MLAB_ROOT}/MeVisLab/Standard/Sources/ML/MLBase
  ${MLAB_ROOT}/MeVisLab/Standard/Sources/ML/MLTools
  ${MLAB_ROOT}/MeVisLab/Standard/Sources/ML/MLTools/include
  ${MLAB_ROOT}/MeVis/Foundation/Sources/SharedIncludes
  ${MLAB_ROOT}/MeVis/Foundation/Sources/ML/include/deprecated
  ${MLAB_ROOT}/MeVis/Foundation/Sources/ML/include
  ${MLAB_ROOT}/MeVis/Foundation/Sources/MLLinearAlgebra/deprecated
  ${MLAB_ROOT}/MeVis/Foundation/Sources/MLLinearAlgebra
  ${MLAB_ROOT}/MeVis/Foundation/Sources/MLProfilingManager/Include
  ${MLAB_ROOT}/MeVis/Foundation/Sources/MLUtilities
  ${MLAB_ROOT}/MeVis/Foundation/Sources/MLMemoryManager/Include
  ${MLAB_ROOT}/FMEwork/ITK/Sources/ITK/MLITK/ITKSupport/
  ${MLAB_ROOT}/MeVis/ThirdParty/Sources/misc/boost
  #${ML}/MeVis/Foundation/Sources/ColumbusCore
  #${ML}/MeVis/Foundation/Sources/MLOpenGL
  #${ML}/MeVis/Foundation/Sources/MLMacOSSupport
  #${ML}/MeVis/ThirdParty/Sources/misc/glew/include
  ${EMCPACKAGE}/
)
                     
LINK_DIRECTORIES( 
  ${BIGRCarotid_BINARY_DIR}/utilities 
  ${MLAB_ROOT}/MeVisLab/Standard/lib
  ${MLAB_ROOT}/MeVis/Foundation/lib
  ${MLAB_ROOT}/FMEwork/lib
  ${EMCPACKAGE}/General/lib
)

ADD_EXECUTABLE(makeMPRStack makeMPRStack.cxx)

SET(ML_LINK_LIBS 
  utilities 
optimized  MLBase
optimized  MLTools
optimized  ML
optimized  MLLinearAlgebra
optimized  MLMemoryManager
optimized  MLUtilities
optimized  MLDicomTree_OFFIS 
debug  MLBase_d
debug  MLTools_d
debug  ML_d
debug  MLLinearAlgebra_d
debug  MLMemoryManager_d
debug  MLUtilities_d
debug  MLDicomTree_OFFIS_d
)

TARGET_LINK_LIBRARIES( makeMPRStack ITKCommon ${ML_LINK_LIBS})

#glu32.lib
#opengl32.lib
#gdi32.lib
#user32.lib
#MLCSO.lib
#glew.lib
#ITKIO
#ML
#MLUtilities
#MLLinearAlgebra
#MLMemoryManager
#MLProfilingManager
#MLOpenGL.lib
#MLProfilingManager.lib
