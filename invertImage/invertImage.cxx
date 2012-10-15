
#include "itkInvertImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include <iostream>
#include <string>
#include "boost/program_options.hpp"
#include "itkImageConstIterator.h"
using namespace std;
// Boost program options
namespace po = boost::program_options;
/** Usage:
 * testSubtracter <inputImageFileName> <outputImageFileName>
 */
int main( int argc, char** argv )
{

  /** Read command line argument */
  std::vector<std::string> help;
  help.push_back(" invertImage");
  help.push_back(" Perform image inverting and constructing cost image from medialness or vesselness");
  help.push_back(" Author: Hui Tang");

  po::options_description desc("Parameters");
  desc.add_options()
    ("help,h", "produce help message")
    ("inputImage,i", po::value<std::string>(), "input image, .mhd file ")
    ("outputImage,o", po::value<std::string>()->default_value("cost.mhd"), "output image")
    ("exponential,e", po::value<float>(), "exponential")
    ("constant,c", po::value<double>(), "constant");


  // Parse command line options
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  po::notify(vm);

  // Output program header
  //programHeader("MedialnessImageFilter", "Hui Tang", 2011);

  // Print help message
  if (vm.count("help") || vm.size()==1) {
    std::cout << "invert image intensity" << std::endl << std::endl;
    std::cout << desc << "\n";
    return EXIT_SUCCESS;
  }

  // String for in- and ouput file
  std::string inputImageN;
  std::string outputImageN;

  // Get intensity file
  if (vm.count("inputImage")) {
    inputImageN = vm["inputImage"].as<std::string>();
    std::cout << "inputImage: " << inputImageN << std::endl;
  } else {
    std::cout << "Error: no image provided" << std::endl;
    return EXIT_FAILURE;
  }

  // Get outputfile
  if (vm.count("outputImage")) {
    outputImageN = vm["outputImage"].as<std::string>();
    std::cout << "Output Image: " << outputImageN << std::endl;
  } else {
    std::cout << "Error: no output file provided" << std::endl;
    return EXIT_FAILURE;
  }

  float e = 1.0f;
  if (vm.count("v")) {
    e = vm["exponential"].as< float >();
    std::cout << "exponential: " << e << std::endl;
  } else {
    std::cout << "Error: no exponential provided" << std::endl;
    return EXIT_FAILURE;
  }

  double c = 0.0001f;
  if (vm.count("v")) {
    c = vm["constant"].as< double >();
    std::cout << "constant: " << c << std::endl;
  } else {
    std::cout << "Error: no constant provided" << std::endl;
    return EXIT_FAILURE;
  }

  /** Define image type */
  const unsigned int Dimension = 3;
  typedef float PixelType;
  typedef itk::Image< PixelType, Dimension >   ImageType;

  /** Define image reader, writer and image filter */
  typedef itk::ImageFileReader< ImageType >    ReaderType;
  typedef itk::ImageFileWriter< ImageType >    WriterType;
  typedef itk::InvertImageFilter<
    ImageType, ImageType >                     InverterType;

  /** Create image reader, writer and filter */
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();
  InverterType::Pointer inverter = InverterType::New();

  /** Set up reader and filter */
  reader->SetFileName( inputImageN );
  inverter->SetInput( reader->GetOutput() );
  inverter->SetConstant( c );
  inverter->SetExponential( static_cast<int>(e) );
  writer->SetFileName( outputImageN);
  writer->SetInput( inverter->GetOutput() );

  try
  {
    writer->Update();
  }
  catch ( itk::ExceptionObject & err )
  {
    std::cerr << "ERROR: " << err << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;

} // end main
