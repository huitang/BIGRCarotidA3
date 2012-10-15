/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#define BOOST_EXCEPTION_DISABLE 1
#include "itkRegionBasedLevelSetFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"
#include "itkVersorRigid3DTransform.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkLevelSetDomainMapImageFilter.h"
#include "itkLevelSetContainerBase.h"
#include "itkLevelSetEquationChanAndVeseTerm.h"
#include "itkLevelSetEquationRSFTerm.h"
#include "itkHeavisideImageFilter.h"
#include "itkLevelSetEquationChanAndVeseInternalTerm.h"
#include "itkLevelSetEquationChanAndVeseExternalTerm.h"
#include "itkLevelSetEquationCurvatureTerm.h"
#include "itkLevelSetEquationGeodesicCurvatureTerm.h"
#include "itkLevelSetEquationTermContainerBase.h"
#include "itkLevelSetEquationContainerBase.h"
#include "itkAtanRegularizedHeavisideStepFunction.h"
#include "itkLevelSetEvolution.h"
#include "itkWhitakerSparseLevelSetImage.h"
#include "itkLevelSetDenseImageBase.h"
#include "itkLevelSetContainer.h"
#include "itkLevelSetEvolutionNumberOfIterationsStoppingCriterion.h"
#include "itkBinaryImageToLevelSetImageAdaptor.h"
#include "itkSeedPointFileIO.h"
//#include "itkNumericTraits.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkLevelSetEquationCurvatureTerm.h"
#include "itkLevelSetEquationAdvectionTerm.h"
//#include "itkWhitakerCommandIterationUpdate.h"
//#include "itkShiCommandIterationUpdate.h"
#include "itkGeodesicActiveContourLevelSetImageFilter.h"
#include "itkLevelSetContainer.h"
#include "itkLevelSetEquationLaplacianTerm.h"
#include "itkAtanRegularizedHeavisideStepFunction.h"
#include "itkBinaryImageToLevelSetImageAdaptor.h"
#include "itkLevelSetEvolutionNumberOfIterationsStoppingCriterion.h"

//#include "itkWhitakerCommandIterationUpdate.h"
#include "boost/program_options.hpp"
#include "boost/filesystem/operations.hpp"
#include "itkScalarChanAndVeseSparseLevelSetImageFilter.h"
#include "itkIdentityTransferFunction.h"
#include <time.h>
#include <math.h>
#include "../../common/utils.h"
#include "../../common/path.h"
#include "../../common/pathio.h"


// ------------------------------------------------------------------------
//
// Exercise: Add a curvature term for regularization.
//
// The coefficient for the new term will be provided as an argument to the
// executable.
//
// ------------------------------------------------------------------------

using namespace std;
// Boost program options
namespace po = boost::program_options;
namespace bfs = boost::filesystem;

std::vector<float> normalizeVector(const std::vector<float> inputVec)
{
	float length=sqrt(pow(inputVec[0],2)+pow(inputVec[1],2)+pow(inputVec[2],2));
std::vector<float> outputVec(3,0);
outputVec[0]=inputVec[0]/length;
outputVec[1]=inputVec[1]/length;
outputVec[2]=inputVec[2]/length;
return outputVec;
}

float length(const std::vector<float> inputVec)
{
	float length=sqrt(pow(inputVec[0],2)+pow(inputVec[1],2)+pow(inputVec[2],2));
return length;
}

float dotProduct(const std::vector<float> inputVec1,const std::vector<float> inputVec2)
{
	float length=inputVec1[0]*inputVec2[0]+inputVec1[1]*inputVec2[1]+inputVec1[2]*inputVec2[2];
	return length;
}
std::vector<float> crossProduct(const std::vector<float> inputVec1,const std::vector<float> inputVec2)
{
	 std::vector<float> outputVec(3,0);
	 outputVec[0]=(inputVec1[1]*inputVec2[2])-(inputVec1[2]*inputVec2[1]);
	 outputVec[1]=(inputVec1[2]*inputVec2[0])-(inputVec1[0]*inputVec2[2]);
	 outputVec[2]=(inputVec1[0]*inputVec2[1])-(inputVec1[1]*inputVec2[0]);
	 return outputVec;

}
template< class TImage >
std::vector<float> getShiftFromSliceCenter(typename TImage::PointType pos,std::vector<float> axisx,std::vector<float> axisy,typename TImage::Pointer CMPRImage)
{
/*vec3 voxelsize=getUpdatedInImg(0)->getVoxelSize();
 double shiftX2D=(index[0]-(getUpdatedInImg(0)->getImgExt().x/2))*voxelsize[0];
 double shiftY2D=(index[1]-(getUpdatedInImg(0)->getImgExt().y/2))*voxelsize[1];*/
InputImageType::IndexType index;
CMPRImage->TransformPhysicalPointToIndex(pos, index);
std::vector<float> centerV(index[0]/2,index[1]/2,index[2]);
std::vector<float> centerW;
getInImg(0)->TransformIndexToPhysicalPoint(centerV,centerW);
std::vector<float> shift2D;
shift2D[0]= pos[0] -centerW[0];
shift2D[1]= pos[1] -centerW[1];
shift2D[2]= pos[2] -centerW[2];
shift=shift2D[0]*axisx+shift2D[1]*axisy;
return shift;
}

int main( int argc, char* argv[] )
{
  std::vector<std::string> parameters;
  std::vector<std::string> help;
  std::vector<float>  parametersF;
  help.push_back("Convert a centerline generated in CMPR back to origin");
  help.push_back("");
  help.push_back(" Author: Hui Tang");

  po::options_description desc1("Images");
  desc1.add_options()
    ("help,h", "produce help message")
	("CMPRimage,i", po::value<std::string>(), "input CMPR image, .mhd file ")
	("CMPRCenterline,c", po::value<std::string>(), "input CMPR centerline, .mhd file ")
	("OriginalCenterline,r", po::value<std::string>(), "input original centerline (please be the same as the one used to generate the input CMPR image), .mhd file ")
	("OutputOriginalCenterline,o", po::value<std::string>(), "output original centerline");

  po::options_description desc2("Parameters");
  desc2.add_options()
	  ("parameters,P",po::value<std::vector<std::string> >(), "parameters");
  po::positional_options_description pd; 
  pd.add("parameters", -1);
  po::options_description cmddescAll("Parse");
  cmddescAll.add(desc1).add(desc2);

  po::variables_map vm;
  po::parsed_options parsed=po::command_line_parser(argc, argv).options(cmddescAll).allow_unregistered().positional(pd).run();
  //po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);

  try
  {
	  po::store(parsed, vm);
	  po::notify(vm);
  }
  catch (std::exception& e)
  {
	  cout << "Exception while parsing parameters: " << endl;
	  cout << e.what() << endl;
	  exit(-1);
  }



  // Print help message
  if (vm.count("help") || vm.size()==1) {
	  std::cout << "Convert a CMPR centerline back to original centerline" << std::endl ;
	  std::cout << " "<< std::endl;
	  std::cout << "Example:ConvertCMPRBackToOrigin.exe -i CMPR.mhd -c CMPRCenterline.mhd -r OriginalCenterline.xml -o OutputOriginalCenterline.xml"<< std::endl;
	  return EXIT_SUCCESS;
  }

  std::vector<std::string> values(vm["parameters"].as<std::vector<std::string> >()); 
  for(std::vector<std::string>::iterator iter = values.begin(); values.end() != iter; ++iter) 
  {
	  parameters.push_back(*iter);
	  std::string a = *iter;
	  std::istringstream b(a);
	  float f;
	  b >> f;
	  parametersF.push_back(f);
  }



  std::string CMPRImageN;
  std::string CMPRCenterlinePointN ;
  std::string OriginalCenterlinePointN;
  std::string OutputOriginalCMPRCenterlinePointN;
  // Get intensity file
  if (vm.count("CMPRimage")) {
	  CMPRImageN = vm["CMPRimage"].as<std::string>();
	  std::cout << "CMPRimage: " << CMPRImageN<< std::endl;
  } else {
	  std::cout << "Error: no original image provided" << std::endl;
	  return EXIT_FAILURE;
  }
  const unsigned int Dimension = 3;
  typedef float                                    InputPixelType;
  typedef itk::Image< InputPixelType, Dimension >           InputImageType;
    typedef itk::Image< InputPixelType, 2>           OutputImageType;
  typedef itk::ImageFileReader< InputImageType >            ReaderType;
 
  ReaderType::Pointer initialReader = ReaderType::New();
  initialReader->SetFileName( CMPRImageN );
  initialReader->Update();
  InputImageType::Pointer CMPRInput = initialReader->GetOutput();
  std::vector<float> imgExt(3, 0);
  imgExt[0] = CMPRInput->GetBufferedRegion().GetSize()[0];
  imgExt[1] = CMPRInput->GetBufferedRegion().GetSize()[1];
  imgExt[2] = CMPRInput->GetBufferedRegion().GetSize()[2];
  InputImageType::PointType originIn;
  originIn= CMPRInput->GetOrigin();
  std::cout<<"Origin"<< originIn[0]<< originIn[1]<< originIn[2]<<std::endl;  
  std::cout<<"Read "<< CMPRImageN<<"done"<<std::endl;
  std::cout<<"Size"<<imgExt[0]<<imgExt[1]<<imgExt[2]<<std::endl;  

  if (vm.count("CMPRCenterline")) {
	  CMPRCenterlinePointN = vm["CMPRCenterline"].as<std::string>();
	  std::cout << "CMPRCenterline: " <<CMPRCenterlinePointN << std::endl;
  } else {
	  std::cout << "Error: no CMPRCenterline name provided" << std::endl;
	  return EXIT_FAILURE;
  }

    if (vm.count("OriginalCMPRSeedPointN")) {
	  OriginalCenterlinePointN = vm["OriginalCMPRCenterlinePointN"].as<std::string>();
	  std::cout << "OriginalCMPRCenterlinePoint: " <<OriginalCenterlinePointN << std::endl;
  } else {
	  std::cout << "Error: no OriginalCMPRCenterlinePoint name provided" << std::endl;
	  return EXIT_FAILURE;
  }

	  if (vm.count("OutputOriginalCMPRCenterlinePointN")) {
	  OutputOriginalCMPRCenterlinePointN = vm["OutputOriginalCMPRCenterlinePointN"].as<std::string>();
	  std::cout << "OutputOriginalCMPRCenterlinePoint: " <<OutputOriginalCMPRCenterlinePointN << std::endl;
  } else {
	  std::cout << "Error: no output original CMPRCenterlinePoint name provided" << std::endl;
	  return EXIT_FAILURE;
  }

  itk::SeedPointFileIO::PointListType CMPRPoints;
  itk::SeedPointFileIO::Pointer seedPointReader1 = itk::SeedPointFileIO::New();
  seedPointReader1->SetFileName( CMPRCenterlinePointN );
  seedPointReader1->SetVerbose( false );
  CMPRPoints = seedPointReader1->GetPoints();
  std::cout<<"CMPR Centerline contains "<<CMPRPoints.size()<<"Points"<<std::endl;
  if ( CMPRPoints.size() < 1 ) {
	  std::cout << "At least one seed points are needed, but only " << CMPRPoints.size() << " are found." << std::endl;
	  return EXIT_FAILURE;
  }

  
  itk::SeedPointFileIO::PointListType OriginalPoints;
  itk::SeedPointFileIO::Pointer seedPointReader = itk::SeedPointFileIO::New();
  seedPointReader->SetFileName( OriginalCenterlinePointN );
  seedPointReader->SetVerbose( false );
  OriginalPoints = seedPointReader->GetPoints();
  std::cout<<"Original Centerline contains "<<OriginalPoints.size()<<"Points"<<std::endl;
  if ( OriginalPoints.size() < 1 ) {
	  std::cout << "At least one seed points are needed, but only " << OriginalPoints.size() << " are found." << std::endl;
	  return EXIT_FAILURE;
  }
  if (!(OriginalPoints.size()==imgExt[2]))
  {
	  std::cout << "Please make sure that the input original centerline is the same as the one which is used to generate the CMPR image!" << std::endl;
	  return EXIT_FAILURE;
  }

//here I first calculate all the key frame lists and save them in std::vector<InputImageType::DirectionType> keyFrams;

      std::vector<InputImageType::DirectionType> keyFrames;

	  //this is used for generate the first key frame(See Theo's email)
	  std::vector<float> tangential1(3, 0);
	  tangential1[0] = OriginalPoints[2][0]- OriginalPoints[0][0];
	  tangential1[1] = OriginalPoints[2][1]- OriginalPoints[0][1];
	  tangential1[2] = OriginalPoints[2][2]- OriginalPoints[0][2];
	  tangential1=normalizeVector(tangential1);

  for (int i=0;i<OriginalPoints.size();i++)
  {
	  InputImageType::IndexType seedIndexFore;
	  InputImageType::IndexType seedIndexBack;
	  int fore=min(i+1,(int)(OriginalPoints.size()-1));
	  int back=max(i-1,0);
	  CMPRInput->TransformPhysicalPointToIndex( OriginalPoints[fore],seedIndexFore );
	  CMPRInput->TransformPhysicalPointToIndex( OriginalPoints[back],seedIndexBack );
      std::vector<float> tangential(3, 0);
	  std::vector<float> localy;
	  std::vector<float> localx;
	  tangential[0] = OriginalPoints[fore][0]- OriginalPoints[back][0];
	  tangential[1] = OriginalPoints[fore][1]- OriginalPoints[back][1];
	  tangential[2] = OriginalPoints[fore][2]- OriginalPoints[back][2];
	  tangential=normalizeVector(tangential);
	  //tangential.normalize();
      // local_y_i = t_i x local_x_i-1   (thus current tangent cross previous local_x)
      //local_x_i = local_y_i x t_i
	  if(i==0)
	  {
		  // prevent parallel
	    localy= crossProduct(tangential,tangential1);
	    if (length(localy)>0)
	    {
	      localy=normalizeVector(localy);
	     }
	    else
	    {
	    	  std::vector<float> localyN(3,0);
	    	  localyN[0]=-localy[1];
	    	  localyN[1]=-localy[0];
	    	  localyN[2]=0;
		      localy=localyN;
	     }
	    localx=crossProduct(localy,tangential);
	    localx=normalizeVector(localx);
	  }

	  else
	  {
	  localy= crossProduct(tangential,localx);
	  localy=normalizeVector(localy);
	  localx= crossProduct(localy,tangential);
	  localx=normalizeVector(localx);
	  }
	  InputImageType::DirectionType direction;
	  direction(0,0) = localx[0];
	  direction(1,0) = localx[1];
	  direction(2,0) = localx[2];
	  direction(0,1) = localy[0];
	  direction(1,1) = localy[1];
	  direction(2,1) = localy[2];
	  direction(0,2) = tangential[0];
	  direction(1,2) = tangential[1];
	  direction(2,2) = tangential[2];
	  keyFrames[i]=direction;

  }
 pathtype outputCenterline;
 for (int i=0;i<CMPRPoints.size();i++)
  {         
	    InputImageType::IndexType index;
        InputImageType::PointType pMPR=CMPRPoints[i];
	    CMPRInput->TransformPhysicalPointToIndex( CMPRPoints[i], index );
        std::vector<float> averagedNorm;
		averagedNorm[0]=keyFrames[i](0,2);
		averagedNorm[1]=keyFrames[i](1,2);
		averagedNorm[2]=keyFrames[i](2,2);
        std::vector<float> averagedAxisY;
		averagedAxisY[0]=keyFrames[i](0,1);
		averagedAxisY[1]=keyFrames[i](1,1);
		averagedAxisY[2]=keyFrames[i](2,1);
        std::vector<float> averagedAxisX;
		averagedAxisX[0]=keyFrames[i](0,0);
		averagedAxisX[1]=keyFrames[i](1,0);
		averagedAxisX[2]=keyFrames[i](2,0);
        normalizeVector(averagedAxisX);
        normalizeVector(averagedAxisY);
        std::vector<float> averagedP;
		averagedP[0]=OriginalPoints[i][0];
		averagedP[1]=OriginalPoints[i][1];
		averagedP[2]=OriginalPoints[i][2];
        std::vector<float> shift=getShiftFromSliceCenter(pMPR,averagedAxisX,averagedAxisY,CMPRInput);
        std::vector<float> newLocation;
		newLocation[0]=averagedP[0]+shift[0];
		newLocation[1]=averagedP[1]+shift[1];
		newLocation[2]=averagedP[2]+shift[2];
		position p(newLocation[0],newLocation[1],newLocation[2]);
		outputCenterline.push_back(p);
 }
 writePathToFile(OutputOriginalCMPRCenterlinePointN, outputCenterline, false, false);
}