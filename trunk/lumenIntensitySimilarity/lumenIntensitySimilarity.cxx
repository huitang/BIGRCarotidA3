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
#include "../utilities/SeedPointFileIO.h"
#include "itkLumenIntensitySimilarityFilter.h"

// Include boost classes
#include "boost/program_options.hpp"

// Boost program options
namespace po = boost::program_options;

int main (int argc, char ** argv) {

  std::vector<std::string> help;
  help.push_back(" lumenIntensitySimilarity");
  help.push_back(" Calculate the similarity to lumen seed points");
  help.push_back(" Author: Reinhard Hameeteman");

  po::options_description desc("Parameters");
  desc.add_options()
    ("help,h", "produce help message")
    ("inputImage,i",       po::value<std::string>(),                                  "input image ")
    ("seeds,s",            po::value<std::string>(),                                  "seed point file")
    ("radius0,r",          po::value<float>(),                                        "radius for seed0 statistics")
    ("radius1,R",          po::value<float>(),                                        "radius for seed1 statistics")
    ("darkLumen,d",        po::value<int>()->default_value(0),                        "lumen intensity dark (0) or bright(1)")
    ("outputLIS,o",        po::value<std::string>(),                                  "output LIS file")
    ("numberOfThreads,p",  po::value<int>(),                                          "numberOfThread");


  // Parse command line options
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  po::notify(vm);

  // Print help message
  if (vm.count("help") || vm.size()==1) {
    std::cout << "Calculate the lumen intensity similarity." << std::endl << std::endl;
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
  if (vm.count("seeds")) {
    seedPointFileName = vm["seeds"].as<std::string>();
    std::cout << "seed point file name: " << seedPointFileName << std::endl;
  } else {
    std::cout << "No seed point file provided." << std::endl;
    return EXIT_FAILURE;
  }

  int threads = 1;
  if (vm.count("numberOfThreads")) {
    threads = vm["numberOfThreads"].as<int>();
  } 
  std::cout << "numberOfThreads: " << threads << std::endl;

  float radius0 = 3.5;
  if (vm.count("radius0")) {
    radius0 = vm["radius0"].as<float>();
  } 
  std::cout << "radius0 : " << radius0 << std::endl;

  float radius1 = 2.5;
  if (vm.count("radius1")) {
    radius1 = vm["radius1"].as<float>();
  } 
  std::cout << "radius1 : " << radius1 << std::endl;

  int darkLumen = 0;
  if (vm.count("darkLumen")) {
    darkLumen = vm["darkLumen"].as<int>();
  } 
  std::cout << "darkLumen: " << darkLumen << std::endl;

  // Get outputfile
  std::string outputImageFileName;
  if (vm.count("outputLIS")) {
    outputImageFileName = vm["outputLIS"].as<std::string>();
    std::cout << "Output lumen intensity similarity: " << outputImageFileName << std::endl;
  } else {
    std::cout << "Error: no output file name provided" << std::endl;
    //return EXIT_FAILURE;
  }

  itk::MultiThreader::Pointer mt = itk::MultiThreader::New();
  mt->SetGlobalMaximumNumberOfThreads( threads );



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

  if ( seedPoints.size() < 2 ) {
    std::cout << "At least two seed points are needed, but only " << seedPoints.size() << " are found." << std::endl;
    return EXIT_FAILURE;
  }

  typedef itk::LumenIntensitySimilarityFilter< ImageType, ImageType> LisFilterType;
  LisFilterType::RadiusListType radii;
  radii.push_back( radius0 );
  radii.push_back( radius1 );
  LisFilterType::Pointer lisFilter = LisFilterType::New();
  lisFilter->SetInput( inputImageReader->GetOutput() );
  lisFilter->SetSeedList( seedPoints );
  lisFilter->SetDarkLumen( darkLumen );
  lisFilter->SetRadii( radii );
  lisFilter->SetVerbose( true );
  lisFilter->Update();


  // write output image
  typedef  itk::ImageFileWriter< ImageType > WriterType;
  WriterType::Pointer outputWriter = WriterType::New();
  std::cout << "writing output image file..." << std::endl;
  outputWriter->SetFileName(outputImageFileName); 
  outputWriter->SetInput( lisFilter->GetOutput() );
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
