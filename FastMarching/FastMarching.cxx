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
#include "../utilities/itkSeedPointFileIO.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkConstShapedNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkFastMarchingImageFilter.h"
#include <math.h>

// Include boost classes
#include "boost/program_options.hpp"

// Boost program options
namespace po = boost::program_options;

int main (int argc, char ** argv) {

  std::vector<std::string> help;
  help.push_back(" fastMarching");
  help.push_back(" Calculate fastmarching with speed defined on input image");
  help.push_back(" Author: Hui Tang");

  po::options_description desc("Parameters");
  desc.add_options()
    ("help,h", "produce help message")
    ("inputImage,i",       po::value<std::string>(),                                  "input image ")
    ("centerline,c",       po::value<std::string>(),                                  "centerline file")
    ("constantSpeed,C",       po::value<float>(),                                  "constant speed (0/1)")
    ("stopValue,S",       po::value<float>(),                                      "stop value")
	("binaryThreshould,b",       po::value<float>(),                                      "binary threshould,0 output sdm, n, output binary image threshoulded at n")
    ("output,o",       po::value<std::string>(),                                  "output file");


  // Parse command line options
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  po::notify(vm);

  // Print help message
  if (vm.count("help") || vm.size()==0) {
    std::cout << "Calculate fastmarching with speed defined on input image." << std::endl << std::endl;
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

  bool c =false;
  if (vm.count("constantSpeed")) {
    c= vm["constantSpeed"].as<float>();
  } 
  std::cout << "constantSpeed : " << c << std::endl;

  float s=0;
  if (vm.count("stopValue")) {
	  s= vm["stopValue"].as<float>();
  } 
  std::cout << "stopValue : " << s << std::endl;

  float b=0;
  if (vm.count("binaryThreshould")) {
	  b= vm["binaryThreshould"].as<float>();
  } 
  std::cout << "binaryThreshould: " << b << std::endl;


  std::string outputImageFileName;
  if (vm.count("output")) {
	  outputImageFileName = vm["output"].as<std::string>();
	  std::cout << "Output image with mean lumen intensity: " << outputImageFileName << std::endl;
  } else {
	  std::cout << "Error: no output file name provided" << std::endl;
	  //return EXIT_FAILURE;
  }
 

  // Define pixel type, dimension, image type and reader type
  typedef float PixelType;
  const unsigned int Dimension = 3;
  typedef  itk::Image< PixelType, Dimension > ImageType;
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
  std::cout << "Reading finished" << std::endl;
  ImageType::Pointer input = inputImageReader->GetOutput();

  // Get Seed Points
  itk::SeedPointFileIO::PointListType seedPoints;
  itk::SeedPointFileIO::Pointer seedPointReader = itk::SeedPointFileIO::New();
  seedPointReader->SetFileName( seedPointFileName );
  seedPointReader->SetVerbose( false );
  seedPoints = seedPointReader->GetPoints();

  if ( seedPoints.size() < 0 ) {
    std::cout << "At least two seed points are needed, but only " << seedPoints.size() << " are found." << std::endl;
    return EXIT_FAILURE;
  }

  ImageType::RegionType region;
  ImageType::IndexType index;
  ImageType::SizeType size;
  ImageType::Pointer constSpeed=ImageType::New();
  constSpeed->SetRegions( input->GetLargestPossibleRegion() );
  constSpeed->CopyInformation( input );
  constSpeed->Allocate();
  constSpeed->FillBuffer( itk::NumericTraits<PixelType>::One);

  typedef itk::FastMarchingImageFilter< ImageType,ImageType > FastMarchingFilterType;
  typedef FastMarchingFilterType::NodeContainer NodeContainer;
  typedef FastMarchingFilterType::NodeType NodeType;
  NodeContainer::Pointer seeds = NodeContainer::New();
  seeds->Initialize();
  for (int i=0;i<seedPoints.size();i++)
  {
	    ImageType::IndexType seedIndex;
        input->TransformPhysicalPointToIndex( seedPoints[i],seedIndex);
		NodeType node;
		const double seedValue = 0.0;
		node.SetValue( seedValue );
		node.SetIndex( seedIndex );

		seeds->InsertElement( i, node );
  }
  FastMarchingFilterType::Pointer fastMarching = FastMarchingFilterType::New();
  if (c)
{
	  fastMarching->SetInput( constSpeed );
	  fastMarching->SetTrialPoints(seeds);
	  fastMarching->SetStoppingValue(s);
	  fastMarching->SetNormalizationFactor(1);
	  fastMarching->SetOutputSize( input->GetBufferedRegion().GetSize() );
	  fastMarching->Update();


	
	  if (b>0)
	  {
		  
		  typedef itk::BinaryThresholdImageFilter< ImageType,ImageType> ThreshouldType;
		  ThreshouldType::Pointer threshoulder = ThreshouldType::New();
		  threshoulder->SetInput(fastMarching->GetOutput());
		  threshoulder->SetLowerThreshold(0);
		  threshoulder->SetUpperThreshold(b);
		  threshoulder->SetOutsideValue(0);
		  threshoulder->SetInsideValue(1);
		  typedef  itk::ImageFileWriter< ImageType > WriterType;
		  WriterType::Pointer outputWriter = WriterType::New();
		  std::cout << "writing output image file..." << std::endl;
		  outputWriter->SetFileName(outputImageFileName); 
		  outputWriter->SetInput( threshoulder->GetOutput() );
		  try{
			  outputWriter->Update();
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
	  else
	  {
		  typedef  itk::ImageFileWriter< ImageType > WriterType;
		  WriterType::Pointer outputWriter = WriterType::New();
		  std::cout << "writing output image file..." << std::endl;
		  outputWriter->SetFileName(outputImageFileName); 
		  outputWriter->SetInput( fastMarching->GetOutput() );
		  try{
			  outputWriter->Update();
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
}
  else
  {
	  fastMarching->SetInput( input);
	  fastMarching->SetTrialPoints(seeds);
	  fastMarching->SetStoppingValue(s);
	  fastMarching->SetNormalizationFactor(1);
	  fastMarching->SetOutputSize( input->GetBufferedRegion().GetSize() );
	  fastMarching->Update();

	  if (b>0)
	  {

		  typedef itk::BinaryThresholdImageFilter< ImageType,ImageType> ThreshouldType;
		  ThreshouldType::Pointer threshoulder = ThreshouldType::New();
		  threshoulder->SetInput(fastMarching->GetOutput());
		  threshoulder->SetLowerThreshold(0);
		  threshoulder->SetUpperThreshold(b);
		  threshoulder->SetOutsideValue(0);
		  threshoulder->SetInsideValue(1);
		  typedef  itk::ImageFileWriter< ImageType > WriterType;
		  WriterType::Pointer outputWriter = WriterType::New();
		  std::cout << "writing output image file..." << std::endl;
		  outputWriter->SetFileName(outputImageFileName); 
		  outputWriter->SetInput( threshoulder->GetOutput() );
		  try{
			  outputWriter->Update();
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
	  else
	  {
		  typedef  itk::ImageFileWriter< ImageType > WriterType;
		  WriterType::Pointer outputWriter = WriterType::New();
		  std::cout << "writing output image file..." << std::endl;
		  outputWriter->SetFileName(outputImageFileName); 
		  outputWriter->SetInput( fastMarching->GetOutput() );
		  try{
			  outputWriter->Update();
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
      }
}
