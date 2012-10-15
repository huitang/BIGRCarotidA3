/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: pathTracking.cxx,v $
  Language:  C++
  Date:      $Date: 2007/06/12 17:55:24 $
  Version:   $Revision: 1.6 $

=========================================================================*/
#if defined (_MSC_VER)
#pragma warning (disable: 4786)
#endif

// Include ITK classes
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMakeCostImageFilter.h"
#include "SeedPointFileIO.h"
#include "boost/program_options.hpp"

// Boost program options
namespace po = boost::program_options;

int main (int argc, char ** argv) {

  std::vector<std::string> help;
  help.push_back(" medialnessImagefilter");
  help.push_back(" Perform slice based medialness image filter based on Robust Vessel Tree Modeling.");
  help.push_back(" Author: Reinhard Hameeteman");

  std::vector< float > radii;
  po::options_description desc("Parameters");
  desc.add_options()
    ("help,h", "produce help message")
    ("bbImage,b",          po::value<std::string>(),                                  "black blood image file name")
    ("pcImage,c",          po::value<std::string>(),                                  "phase contrast image file name")
    ("output,o",           po::value<std::string>(),                                  "output file name")
    ("seeds,s",            po::value<std::string>(),                                  "seed point file name")
    ("mask,m",             po::value<std::string>(),                                  "mask image file name")
    ("threshold,t",        po::value<float>(),                                        "threshould of mask image")
    ("darkLumen,d",        po::value<float>(),                                        "lumen intensity dark (0) or bright(1)")
    //("sigma,s",            po::value<float>(),                                        "sigma for the gradient used in the medialness calculation")
    ("numberOfAngles,a",   po::value<int>(),                                          "number of angles used in medialness calculation")
    ("minimumRadius,r",    po::value<float>(),                                        "minimum radius used in the medialness calculation")
    ("maximumRadius,R",    po::value<float>(),                                        "maximum radius used in the medialness calculation")
    ("numberOfRadius,n",   po::value<int>(),                                          "number of radii used in the medialness calculation")
    ("radii,x",            po::value< std::vector<float> >(&radii)->multitoken(),     "radius used for the statistics of the seeds used in the lis calculation")
    ("verbose,v",          po::value<bool>(),                                         "turn verbose on")
    ("numberOfThreads,p",  po::value<int>(),                                          "numberOfThread");


  // Parse command line options
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  po::notify(vm);

  // Output program header
  //programHeader("MedialnessImageFilter", "Hui Tang", 2011);
  
  // Print help message
  if (vm.count("help") || vm.size()==1) {
    std::cout << "Perform medialness image filter based on Robust Vessel Tree Modeling." << std::endl << std::endl;
    std::cout << desc << "\n";
    return EXIT_SUCCESS;
  }

  /*
   * Get command line parameters
   */

  // Set number of threads
  int threads = 1;
  if (vm.count("numberOfThreads")) {
    threads = vm["numberOfThreads"].as<int>();
  } 
  std::cout << "numberOfThreads: " << threads << std::endl;
  itk::MultiThreader::Pointer mt = itk::MultiThreader::New();
  mt->SetGlobalMaximumNumberOfThreads( threads );

  // String for in- and ouput file

  std::string bbImageName;
  if (vm.count("bbImage")) {
    bbImageName = vm["bbImage"].as<std::string>();
    std::cout << "black blood image name: " << bbImageName << std::endl;
  } else {
    std::cout << "Error: no black blood image provided" << std::endl;
    return EXIT_FAILURE;
  }

  std::string pcImageName;
  if (vm.count("pcImage")) {
    pcImageName = vm["pcImage"].as<std::string>();
    std::cout << "phase contrast image name: " << pcImageName << std::endl;
  } else {
    std::cout << "no phase contras image provided, using only black blood image" << std::endl;
  }

  std::string maskImageName;
  bool useMask = false;
  if (vm.count("mask")) {
    maskImageName = vm["mask"].as<std::string>();
    std::cout << "mask image name: " << maskImageName << std::endl;
    useMask = true;
	} else {
    std::cout << "No mask image provided using entire image." << std::endl;
  }


  // Get seed point file name
  std::string seedPointFileName;
  if (vm.count("seeds")) {
    seedPointFileName = vm["seeds"].as<std::string>();
    std::cout << "seed point file name: " << seedPointFileName << std::endl;
  } else {
    std::cout << "Error: no seed point file provided" << std::endl;
    return EXIT_FAILURE;
  }

  // Get outputfile
  std::string outputFileName;
  if (vm.count("output")) {
    outputFileName = vm["output"].as<std::string>();
    std::cout << "output image name: " << outputFileName << std::endl;
  } else {
    std::cout << "Error: no output file provided" << std::endl;
    return EXIT_FAILURE;
  }

  bool verbose = false;
  if (vm.count("verbose") ) {
    verbose = true;
    std::cout << "verbose on" << std::endl;
  }

  typedef float PixelType;
  const unsigned int Dimension = 3;
  typedef  itk::Image< PixelType, Dimension > ImageType;
  typedef  itk::ImageFileReader< ImageType > ReaderType;

  // Define pixel type, dimension, image type and reader type
  
  // Load files
  ReaderType::Pointer bbReader = ReaderType::New();
  std::cout << "Reading black blood image..." << std::endl;
  bbReader->SetFileName(bbImageName); 
  try {
    bbReader->Update();
  }
  catch ( itk::ExceptionObject &err){
    std::cout << "Error while reading black blood file !" << std::endl;
    std::cout << err << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "reading finished" << std::endl;
  const ImageType::SpacingType bbSpacing = bbReader->GetOutput()->GetSpacing();
  const ImageType::SizeType bbExtent = bbReader->GetOutput()->GetLargestPossibleRegion().GetSize();
  std::cout << "black blood image spacing: "  << bbSpacing << std::endl;
  std::cout << "black blood image extent: " << bbExtent << std::endl;

  ReaderType::Pointer pcReader = ReaderType::New();
  if ( pcImageName != "") {
    std::cout << "Reading phase contrast image ..." << std::endl;
    pcReader->SetFileName(pcImageName); 
    try{
      pcReader->Update();
    }
    catch ( itk::ExceptionObject &err){
      std::cout << "Error while reading phase contrast file !" << std::endl;
      std::cout << err << std::endl;
      return EXIT_FAILURE;
    }
    std::cout << "reading finished" << std::endl;
    const ImageType::SpacingType pcSpacing = pcReader->GetOutput()->GetSpacing();
    const ImageType::SizeType pcExtent = pcReader->GetOutput()->GetLargestPossibleRegion().GetSize();
    std::cout << "phase contrast image spacing: "  << pcSpacing << std::endl;
    std::cout << "phase contrast image extent: " << pcExtent << std::endl;
  }

  ReaderType::Pointer maskImageReader = ReaderType::New();
  if ( useMask ) {
    std::cout << "Reading mask file..." << std::endl;
    maskImageReader->SetFileName( maskImageName );
    try {
      maskImageReader->Update();
    }
    catch ( itk::ExceptionObject &err){
      std::cout << "Error while reading mask file !" << std::endl;
      std::cout << err << std::endl;
      return EXIT_FAILURE;
    }
    std::cout << "Reading  mask finished" << std::endl;

    const ImageType::SpacingType maskSpacing = maskImageReader->GetOutput()->GetSpacing();
    const ImageType::SizeType maskExtent = maskImageReader->GetOutput()->GetLargestPossibleRegion().GetSize();
    std::cout << "mask image spacing: "  << maskSpacing << std::endl;
    std::cout << "mask image extent: " << maskExtent << std::endl;

    if ( bbExtent != maskExtent) {
      std::cout << "Error:black blood and mask image should have the same size" << std::endl;
      return EXIT_FAILURE;
    }
  }

  // Get Seed Points
  itk::SeedPointFileIO::PointListType seedPoints;
  itk::SeedPointFileIO::Pointer seedPointReader = itk::SeedPointFileIO::New();
  seedPointReader->SetFileName( seedPointFileName );
  seedPointReader->SetVerbose( verbose );
  seedPoints = seedPointReader->GetPoints();
  if (seedPoints.size() < 2 && pcImageName == "" ){
    std::cout << "at least 2 seed points are required, only " << seedPoints.size() << " are given" << std::endl;
    return EXIT_FAILURE;
  }

  typedef itk::MakeCostImageFilter< ImageType, ImageType,ImageType,ImageType> CostImageFilterType;
  CostImageFilterType::Pointer costImageFilter = CostImageFilterType::New();
  costImageFilter->SetBlackBloodImage( bbReader->GetOutput() );
  if (pcImageName != "" ) {
    costImageFilter->SetPhaseContrastImage( pcReader->GetOutput() );
  } else {
    costImageFilter->SetUseOnlyBBImage( true );
  }
  costImageFilter->SetSeedList( seedPoints );
  costImageFilter->SetRadii( radii );
  costImageFilter->SetUseMask( useMask );
  if ( useMask ) {
    costImageFilter->SetMaskImage( maskImageReader->GetOutput() );
  }
  costImageFilter->SetVerbose( verbose );
  costImageFilter->UpdateParameters();

  try {
    costImageFilter->Update();
  }
  catch ( itk::ExceptionObject &err){
    std::cout << "Error while calculating cost image!" << std::endl;
    std::cout << err << std::endl;
    return EXIT_FAILURE;
  }


  typedef itk::ImageFileWriter< ImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputFileName );
  writer->SetInput( costImageFilter->GetOutput() );
  try {
    writer->Update();
  }
  catch ( itk::ExceptionObject &err){
    std::cout << "Error while writing output file !" << std::endl;
    std::cout << err << std::endl;
    return EXIT_FAILURE;
  }


  // Exit program
  return EXIT_SUCCESS;
}
