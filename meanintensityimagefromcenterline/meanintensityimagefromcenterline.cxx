/*=========================================================================

  Program:   Lumen intensity similarity
  Module:    $RCSfile: lumenIntensitySimilarity.cxx,v $
  Language:  C++
  Date:      $Date: 2007/06/12 17:55:24 $
  Version:   $Revision: 1.6 $

=========================================================================*/
#if defined (_MSC_VER)
#pragma warning (disable: 4786)
#endif

// Include ITK classes
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkFastMarchingImageFilter.h"
#include "itkConstShapedNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkThresholdImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "../utilities/itkSeedPointFileIO.h"
#include <math.h>

// Include boost classes
#include "boost/program_options.hpp"

// Boost program options
typedef float PixelType;
const unsigned int Dimension = 3;
typedef  itk::Image< PixelType, Dimension > ImageType;

namespace po = boost::program_options;

int main (int argc, char ** argv) {

  std::vector<std::string> help;
  help.push_back(" meanIntensityImageFromCenterline");
  help.push_back(" Calculate image with estimation of mean intensity of the foregound");
  help.push_back(" Author: Hui Tang");

  po::options_description desc("Parameters");
  desc.add_options()
    ("help,h", "produce help message")
    ("inputImage,i",       po::value<std::string>(),                                  "input image ")
	("sdm,s",       po::value<std::string>(),                                  "sdm image")
    ("centerline,c",       po::value<std::string>(),                                  "centerline file")
    ("radius,r",          po::value<float>(),                                        "radius of ROI")
	("meanRange,n",        po::value<float>(),                                        "mean+-n*sigma")
    ("outputMean,o",       po::value<std::string>(),                                  "output mean file");


  // Parse command line options
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  po::notify(vm);

  // Print help message
  if (vm.count("help") || vm.size()==0) {
    std::cout << "Author: Hui Tang" << "\n";
    std::cout << "Calculate image with estimation of mean intensity of the foregound." << std::endl ;
	std::cout<<"This application is based on Hui Tang, Theo van Walsum Lumen Segmentation for atherosclerotic carotid arteries in CTA, ISBI2012"<< std::endl;
	std::cout<<"You are welcome to use this application if you cite the above paper"<<std::endl;
    std::cout << desc << "\n";
    return EXIT_SUCCESS;
  }

  // set parameters

  // input file name
  std::string inputImageFileName;
  if (vm.count("inputImage")) {
    inputImageFileName = vm["inputImage"].as<std::string>();
    std::cout << "input image: " << inputImageFileName << std::endl;
  } else {
    std::cout << "Error: no input image provided" << std::endl;
    return EXIT_FAILURE;
  }

  std::string seedPointFileName;
  if (vm.count("centerline")) {
    seedPointFileName = vm["centerline"].as<std::string>();
    std::cout << "centerline file name: " << seedPointFileName << std::endl;
  } else {
    std::cout << "No centerlinefile provided." << std::endl;
    return EXIT_FAILURE;
  }

  std::string sdmFileName;
  if (vm.count("sdm")) {
	  sdmFileName = vm["sdm"].as<std::string>();
	  std::cout << "sdm: " << sdmFileName << std::endl;
  } else {
	  std::cout << "Error: no input image provided" << std::endl;
	  return EXIT_FAILURE;
  }

  float radius = 10;
  if (vm.count("radius")) {
    radius= vm["radius"].as<float>();
  } 
  std::cout << "radius : " << radius << std::endl;

  float meanRange = 1;
  if (vm.count("meanRange")) {
    meanRange = vm["meanRange"].as<float>();
  } 
  std::cout << "meanRange : " << meanRange << std::endl;

  std::string outputImageFileName;
  if (vm.count("outputMean")) {
	  outputImageFileName = vm["outputMean"].as<std::string>();
	  std::cout << "Output image with mean lumen intensity: " << outputImageFileName << std::endl;
  } else {
	  std::cout << "Error: no output file name provided" << std::endl;
	  //return EXIT_FAILURE;
  }
 

  // Define pixel type, dimension, image type and reader type

  typedef  itk::ImageFileReader< ImageType > ReaderType;
  
  // Load input image
  ReaderType::Pointer inputImageReader = ReaderType::New();
  std::cout << "Reading input image file..." << std::endl;
  inputImageReader->SetFileName(inputImageFileName); 
  try{
    inputImageReader->Update();
  }
  catch ( itk::ExceptionObject &err){
    std::cout << "Error while reading input file !" << std::endl;
    std::cout << err << std::endl;
    return EXIT_FAILURE;
  }
  ImageType::Pointer original=   inputImageReader->GetOutput();

  ReaderType::Pointer sdmReader = ReaderType::New();
  std::cout << "Reading input image file..." << std::endl;
  sdmReader->SetFileName(sdmFileName); 
  try{
	  sdmReader->Update();
  }
  catch ( itk::ExceptionObject &err){
	  std::cout << "Error while reading input file !" << std::endl;
	  std::cout << err << std::endl;
	  return EXIT_FAILURE;
  }
  std::cout << "Reading finished" << std::endl;
  ImageType::Pointer sdm =   sdmReader->GetOutput();

  std::vector<float> imgExt(3, 0);
  imgExt[0] = original->GetBufferedRegion().GetSize()[0];
  imgExt[1] = original->GetBufferedRegion().GetSize()[1];
  imgExt[2] = original->GetBufferedRegion().GetSize()[2];
     std::cout << "original image Extent " << imgExt[0] << ", " << imgExt[1] << ", " << imgExt[2] << std::endl;

  std::vector<float> imgExt1(3, 0);
  imgExt1[0] = sdm->GetBufferedRegion().GetSize()[0];
  imgExt1[1] = sdm->GetBufferedRegion().GetSize()[1];
  imgExt1[2] = sdm->GetBufferedRegion().GetSize()[2];
 std::cout << "sdm image Extent " << imgExt1[0] << ", " << imgExt1[1] << ", " << imgExt1[2] << std::endl;
 
  if (imgExt1[0]!=imgExt[0]||imgExt1[1]!=imgExt[1]||imgExt1[2]!=imgExt[2])
  {
	  std::cout << "image size should be the same!" << std::endl;
	  return EXIT_FAILURE;
  }

  // Get Seed Points
  itk::SeedPointFileIO::PointListType seedPoints;
  itk::SeedPointFileIO::Pointer seedPointReader = itk::SeedPointFileIO::New();
  seedPointReader->SetFileName( seedPointFileName );
  seedPointReader->SetVerbose( false );
  seedPoints = seedPointReader->GetPoints();

  /*if ( seedPoints.size() < 2 ) {
    std::cout << "At least two seed points are needed, but only " << seedPoints.size() << " are found." << std::endl;
    return EXIT_FAILURE;
  }*/

  for (int i =0; i< seedPoints.size(); i++)
  {
	  ImageType::IndexType seedIndex;
      original->TransformPhysicalPointToIndex( seedPoints[0],seedIndex );
  }

  ImageType::RegionType region;
  ImageType::IndexType index;
  ImageType::SizeType size;
  ImageType::Pointer constSpeed=ImageType::New();
  constSpeed->SetRegions( original->GetLargestPossibleRegion() );
  constSpeed->SetSpacing( original->GetSpacing());
  constSpeed->SetOrigin( original->GetOrigin());
  constSpeed->SetDirection( original->GetDirection());
  //constSpeed->CopyInformation( input );
  constSpeed->Allocate();
  constSpeed->FillBuffer( itk::NumericTraits<PixelType>::One);



  // Calculate Centerline mean
  int nElements=seedPoints.size();
  	      std::cout<<nElements<<std::endl;
		  int validN=0;
  ImageType::PixelType sum = static_cast<ImageType::PixelType>(0);
  float mean = 0.0f;
  ImageType::PixelType sSum = static_cast<ImageType::PixelType>(0);
  float stdev = 0.0f;
  float meanofSquare=0.0f;
  for (int i=0;i<nElements;i++)
  {
	  ImageType::IndexType seedIndex;
	  original->TransformPhysicalPointToIndex( seedPoints[i],seedIndex);
	  if(seedIndex[0]<imgExt[0]-1&&seedIndex[0]>0&&seedIndex[1]>0&&seedIndex[1]<imgExt[1]-1&&seedIndex[2]>0&&seedIndex[2]<imgExt[2]-1)
	  {
	  sum=sum+original->GetPixel(seedIndex);
	  sSum= sSum+original->GetPixel(seedIndex)*original->GetPixel(seedIndex);
	 // std::cout<<i<<std::endl;
	  validN++;
	  }

	
  }
    std::cout<<sum<<std::endl;
  mean=sum/validN;
  meanofSquare=sSum/validN;
  std::cout<<validN<<std::endl;
   // std::cout<<215<<std::endl;
  float var = sSum/(validN-1) - (sum/validN)*(sum/validN)*validN/(validN-1);
     // std::cout<<217<<std::endl;
  stdev = sqrt(var);


  std::cout<<"mean"<<mean<<"stdev"<<stdev<<std::endl;
  typedef itk::ImageRegionIteratorWithIndex< ImageType > InputIteratorType;
  std::vector<float> centerlineIntensity;
  for (int i=0; i<seedPoints.size();i++)
  {
	  ImageType::IndexType seedIndex;
	  original->TransformPhysicalPointToIndex( seedPoints[i],seedIndex);
	  if(seedIndex[0]<imgExt[0]-1&&seedIndex[0]>0&&seedIndex[1]>0&&seedIndex[1]<imgExt[1]-1&&seedIndex[2]>0&&seedIndex[2]<imgExt[2]-1)
	  {
		  InputIteratorType centerlineIndex( original, original->GetLargestPossibleRegion() );
		  centerlineIndex.SetIndex(seedIndex);
	  float aveIntensity=0;
	  int N=0;
	  for (int i=-2;i<=2;i++){
		  for (int j=-1;j<=2;j++){
			  for(int k=-1;k<=2;k++){
				  ImageType::IndexType seedIndexTemp;			 
				  seedIndexTemp[0]=seedIndex[0]+i;
				  seedIndexTemp[1]=seedIndex[1]+j;
				  seedIndexTemp[2]=seedIndex[2]+k;
				  if(seedIndexTemp[0]<imgExt[0]-1&&seedIndexTemp[0]>0&&seedIndexTemp[1]>0&&seedIndexTemp[1]<imgExt[1]-1&&seedIndexTemp[2]>0&&seedIndexTemp[2]<imgExt[2]-1)
				  {	
					  aveIntensity= aveIntensity+original->GetPixel(seedIndexTemp);
					  N++;
				  }
			  
			  }
		  }
	  }
		  aveIntensity=aveIntensity/N;
		  centerlineIntensity.push_back(aveIntensity);
	  }


  }

  //here to save memory, I use copy constSpeed to generate the mean image

  InputIteratorType meanIt( constSpeed, constSpeed->GetLargestPossibleRegion() );
  InputIteratorType inputIt( original, original->GetLargestPossibleRegion() );
  InputIteratorType dmIt( sdm,sdm->GetLargestPossibleRegion() );
  dmIt.GoToBegin();
  meanIt.GoToBegin();
  inputIt.GoToBegin();
  while( (!meanIt.IsAtEnd()))
  {

	  const PixelType currentInputIntensity  = inputIt.Get();
	  const PixelType currentDisToCenterline  = dmIt.Get();
	  PixelType outMean = static_cast<ImageType::PixelType>(0);
	  if (currentDisToCenterline<radius)
	  {
		  ImageType::IndexType voxelIndex= inputIt.GetIndex();
		  ImageType::PointType voxelPosition;
		  original->TransformIndexToPhysicalPoint( voxelIndex,voxelPosition);
		  float minD=1000000000000;
		  for (int i=0; i <seedPoints.size(); i++)
		  {
			  ImageType::IndexType seedIndex;
			  original->TransformPhysicalPointToIndex( seedPoints[i],seedIndex);
			  if(seedIndex[0]<imgExt[0]-1&&seedIndex[0]>0&&seedIndex[1]>0&&seedIndex[1]<imgExt[1]-1&&seedIndex[2]>0&&seedIndex[2]<imgExt[2]-1)
			  {
             // InputIteratorType centerlineIndex( original, original->GetLargestPossibleRegion() );
			//  centerlineIndex.SetIndex(seedIndex);
       //std::cout<<aveIntensity<<std::endl;
			  //float currentD =sqrt((float)(seedIndex[0]-voxelIndex[0])*(seedIndex[0]-voxelIndex[0])*original->GetSpacing()[0]*original->GetSpacing()[0]+(seedIndex[1]-voxelIndex[1])*(seedIndex[1]-voxelIndex[1])*original->GetSpacing()[1]*original->GetSpacing()[1]+(seedIndex[2]-voxelIndex[2])*(seedIndex[2]-voxelIndex[2]))*original->GetSpacing()[2]*original->GetSpacing()[2];
			  float currentD =sqrt((voxelPosition[0]-seedPoints[i][0])*(voxelPosition[0]-seedPoints[i][0])+(voxelPosition[1]-seedPoints[i][1])*(voxelPosition[1]-seedPoints[i][1])+(voxelPosition[2]-seedPoints[i][2])*(voxelPosition[2]-seedPoints[i][2]));
			  if (minD>currentD&&centerlineIntensity[i]<mean+ meanRange*stdev&& centerlineIntensity[i]>mean- meanRange*stdev)
			  {
				  minD=currentD;

				  outMean=centerlineIntensity[i];

			  }
			  }

		  }


	  }
	  else
	  {
		  outMean=mean;
	  }
	      meanIt.Set(outMean);
	      ++inputIt;
		  ++meanIt;
		  ++dmIt;

  }

  // write output image

  typedef  itk::ImageFileWriter< ImageType > WriterType;
  WriterType::Pointer outputWriter2 = WriterType::New();
  std::cout << "writing output image file..." << std::endl;
  outputWriter2->SetFileName(outputImageFileName); 
  outputWriter2->SetInput( constSpeed );
  try{
    outputWriter2->Update();
  }
  catch ( itk::ExceptionObject &err){
    std::cout << "Error while writing output file !" << std::endl;
    std::cout << err << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "writing finished" << std::endl;

  // Exit program
  return EXIT_SUCCESS;
}
