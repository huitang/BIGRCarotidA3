
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
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include <math.h>

// Include boost classes
#include "boost/program_options.hpp"

// Boost program options
namespace po = boost::program_options;

int main (int argc, char ** argv) {

  std::vector<std::string> help;
  help.push_back(" gaussianGradientMagnitude3D");
  help.push_back(" Calculate gaussian gradient magnitude");
  help.push_back(" Author: Hui Tang");

  po::options_description desc("Parameters");
  desc.add_options()
    ("help,h", "produce help message")
    ("inputImage,i",       po::value<std::string>(),                                  "input image ")
    ("sigma,s",          po::value<float>(),                                          "sigma of kernal")
    ("outputGradient,o",        po::value<std::string>(),                                  "output gradient magnitude file");

  // Parse command line options
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  po::notify(vm);

  // Print help message
  if (vm.count("help") || vm.size()==0) {
	  std::cout<<"Author: Hui Tang"<<std::endl;
    std::cout << "Calculate gaussian gradient magnitude." << std::endl << std::endl;
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

  float sigma = 1;
  if (vm.count("sigma")) {
    sigma = vm["sigma"].as<float>();
  } 
  std::cout << "sigma : " << sigma << std::endl;

  std::string outputImageFileName;
  if (vm.count("outputGradient")) {
	  outputImageFileName = vm["outputGradient"].as<std::string>();
	  std::cout << "Output gradient: " << outputImageFileName << std::endl;
  } else {
	  std::cout << "Error: no output file name provided" << std::endl;
	  //return EXIT_FAILURE;
  }

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

  typedef itk::GradientMagnitudeRecursiveGaussianImageFilter <ImageType> GradientMagnitudeType;
  GradientMagnitudeType::Pointer gradientMagnitude =GradientMagnitudeType::New();
  gradientMagnitude->SetInput(input);
  gradientMagnitude->SetSigma(sigma);
  gradientMagnitude->Update();

  
  // write output image
  typedef  itk::ImageFileWriter< ImageType > WriterType;
  WriterType::Pointer outputWriter = WriterType::New();
  std::cout << "writing output image file..." << std::endl;
  outputWriter->SetFileName(outputImageFileName); 
  outputWriter->SetInput(  gradientMagnitude->GetOutput() );
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
