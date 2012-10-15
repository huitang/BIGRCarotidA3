/*=========================================================================

  Program:   Path tracking
  Module:    $RCSfile: pathTracking.cxx,v $
  Language:  C++
  Date:      $Date: 2007/06/12 17:55:24 $
  Version:   $Revision: 1.6 $

=========================================================================*/
#if defined (_MSC_VER)
#pragma warning (disable: 4786)
#endif



//#undef ML_DISABLE_CPP
//#define ML_EXPORTS

// ML
#ifdef WIN32
  #define ML_LIBRARY_EXPORT_ATTRIBUTE __declspec(dllexport)
  #define ML_LIBRARY_IMPORT_ATTRIBUTE __declspec(dllimport)
//   #define _WINDOWS
//   #define UNICODE
//   #define WIN32
//   #define QT_LARGEFILE_SUPPORT
//   #define MEVIS_64BIT
//   #define MEVISLAB_VERSION 202
//   #define ML_VERSION_2
//   #define BOOST_ALL_DYN_LINK
//   #define DEBUG
//   #define MeVisLab
//   #define MEVISLAB
//   #define __LITTLE_ENDIAN__
//   #define NOMINMAX
//   #define _CRT_SECURE_NO_DEPRECATE
#else
  #define ML_LIBRARY_EXPORT_ATTRIBUTE __attribute__((__visibility__("default")))
  #define ML_LIBRARY_IMPORT_ATTRIBUTE 
#endif

#define MEVIS_TARGET BIGRCarotid

#include "itkMakeMPRStackImageFilter.h"

#include "boost/program_options.hpp"

// Boost program options
// Include boost classes
namespace po = boost::program_options;

// Include ITK classes
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"


int main (int argc, char ** argv) {

  // Directories
  std::string inputdir, outputdir;
  size_t id;

  // Declare the supported options.
  po::options_description desc("Parameters");
  desc.add_options()
    ("help,h", "produce help message")
    ("inputImage,i", po::value<std::string>(), "input image filename")
    ("outputfile,o", po::value<std::string>(), "output image file")
    ("path,s",       po::value<std::string>(), "path point file");

  // Parse command line options
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).
          options(desc).run(), vm);
  po::notify(vm);

  // Print help message
  if (vm.count("help") || vm.size()==1) {
    std::cout << "Perform pathtracking using cost image in CMPR, no seed points needed as input, the seed points are selected in the slice center of front and end slice ." << std::endl << std::endl;
   std::cout << desc << "\n";
    return EXIT_SUCCESS;
  }

  std::string costfile;
  // Get intensity file
  if (vm.count("inputImage")) {
    costfile = vm["inputImage"].as<std::string>();
    std::cout << "Input image file: " << costfile << std::endl;
  } else {
    std::cout << "Error: no input image file provided" << std::endl;
    return EXIT_FAILURE;
  }

  // Seeds filename
  std::string pathFileName;
  if (vm.count("path")) {
    pathFileName = vm["path"].as<std::string>();
    std::cout << "path point file name: " << pathFileName << std::endl;
  } else {
    std::cout << "No path point file provided" << std::endl;
    return EXIT_FAILURE;
  }

  itk::SeedPointFileIO::PointListType pathPoints;
  itk::SeedPointFileIO::Pointer pathPointReader = itk::SeedPointFileIO::New();
  pathPointReader->SetFileName( pathFileName );
  pathPointReader->SetVerbose( true );
  pathPoints = pathPointReader->GetPoints();

  // Get outputfile
  std::string outputfile;
  if (vm.count("outputfile")) {
    outputfile = vm["outputfile"].as<std::string>();
    std::cout << "Output file: " << outputfile << std::endl;
  } else {
    std::cout << "Error: no output file provided" << std::endl;
    return EXIT_FAILURE;
  }

  // Define pixel type, dimension, imagetype and reader type
  typedef float PixelType;
  const unsigned int Dimension = 3;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;

  // Load files
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( costfile.c_str() );
  reader->Update();
  std::cout << "Reading file..." << std::endl;

 
  typedef itk::MakeMPRStackImageFilter< ImageType, ImageType> MPRFilterType;
  MPRFilterType::Pointer mprFilter = MPRFilterType::New();
  mprFilter->SetInput( reader->GetOutput() );
  mprFilter->SetPath( pathPoints );
  mprFilter->Update();


  // Exit program
  return 0;
}

