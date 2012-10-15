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

// Include ITK classes
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"

// Include path classes
#include "pathTrackingLibrary.h"
#include "SeedPointFileIO.h"

// Include boost classes
#include "boost/program_options.hpp"

// Boost program options
namespace po = boost::program_options;

int main (int argc, char ** argv) {
  // Directories
  std::string inputdir, outputdir;

  // Declare the supported options.
  po::options_description desc("Parameters");
  desc.add_options()
    ("help,h", "produce help message")
    ("costimage,i",  po::value<std::string>(), "cost image filename")
    ("outputfile,o", po::value<std::string>(), "output file")
    ("seeds,s",      po::value<std::string>(), "seed point file");

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

  // String for in- and ouput file
  std::string costfile, outputfile;

  // Get intensity file
  if (vm.count("costimage")) {
    costfile = vm["costimage"].as<std::string>();
    std::cout << "Cost image file: " << costfile << std::endl;
  } else {
    std::cout << "Error: no cost image file provided" << std::endl;
    return EXIT_FAILURE;
  }

  // Seeds filename
  std::string seedPointFileName;
  if (vm.count("seeds")) {
    seedPointFileName = vm["seeds"].as<std::string>();
    std::cout << "seed point file name: " << seedPointFileName << std::endl;
  } else {
    std::cout << "No seed point file provided, using slice centers." << std::endl;
  }

  // Get outputfile
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
  
  // ROI
  ImageType::RegionType desiredRegion;
  // Start of ROI
  ImageType::IndexType start;
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  desiredRegion.SetIndex( start );

  // Create vector of image extents
  std::vector<int> imgExt(3, 0);

  // Load files
  ReaderType::Pointer reader = ReaderType::New();
  float* data = NULL;
  std::cout << "Reading file..." << std::endl;
  
  // Intensity
  typedef itk::ExtractImageFilter< ImageType, ImageType > FilterType;
  FilterType::Pointer filter0 = FilterType::New();
  reader->SetFileName( costfile.c_str() );
  reader->Update();
  ImageType::SizeType size = reader->GetOutput()->GetBufferedRegion().GetSize();
  desiredRegion.SetSize( size );
  // Get ROI from image
  filter0->SetExtractionRegion( desiredRegion );
  filter0->SetInput( reader->GetOutput() );
  filter0->Update();
  filter0->GetOutput()->GetPixelContainer()->SetContainerManageMemory(false);
  // Set pointer to data
  data = filter0->GetOutput()->GetPixelContainer()->GetImportPointer();
  // Get image extent
  imgExt[0] = reader->GetOutput()->GetBufferedRegion().GetSize()[0];
  imgExt[1] = reader->GetOutput()->GetBufferedRegion().GetSize()[1];
  imgExt[2] = reader->GetOutput()->GetBufferedRegion().GetSize()[2];

  // Output image extents
  std::cout << "Image extent: " << imgExt[0] << ", " << imgExt[1] << ", " << imgExt[2] << std::endl;

  // Get voxelsizes
  const ImageType::SpacingType& sp = reader->GetOutput()->GetSpacing();
  std::cout << "Spacing: "  << sp[0] << ", " << sp[1] << ", " << sp[2] << std::endl;
  std::vector<float> voxelSizes (3, 1.0f);
  voxelSizes[0] = sp[0];
  voxelSizes[1] = sp[1];
  voxelSizes[2] = sp[2];

  // Get Seed Points
  float sx, sy, sz, ex, ey, ez;

  if ( seedPointFileName == "" ) {
    sx=floor((float)imgExt[0]/2.0f);
    sy=floor((float)imgExt[1]/2.0f);
    sz=floor((float)0.0);
    ex=floor((float)imgExt[0]/2.0f);
    ey=floor((float)imgExt[1]/2.0f);
    ez=floor((float)imgExt[2]-1.0f);
  } else{
    itk::SeedPointFileIO::PointListType seedPoints;
    itk::SeedPointFileIO::Pointer seedPointReader = itk::SeedPointFileIO::New();
    seedPointReader->SetFileName( seedPointFileName );
    seedPointReader->SetVerbose( false );
    seedPoints = seedPointReader->GetPoints();

    ImageType::IndexType startSeed;
    ImageType::IndexType endSeed;
    reader->GetOutput()->TransformPhysicalPointToIndex( seedPoints[0], startSeed );
    reader->GetOutput()->TransformPhysicalPointToIndex( seedPoints[1], endSeed );
    std::cout << "start point (in wc): " << startSeed << std::endl;
    std::cout << "end point (in wc):   " << endSeed << std::endl;
    sx = startSeed[0];
    sy = startSeed[1];
    sz = startSeed[2];
    ex = endSeed[0];
    ey = endSeed[1];
    ez = endSeed[2];
  }


  std::cout << "start point (in voxels): " << sx << ", " << sy << ", " << sz << std::endl;
  std::cout << "end point (in voxels): " << ex << ", " << ey << ", " << ez << std::endl;


  std::cout << "Perform path tracking..." << std::endl;
  PathTracking pt(imgExt, data, voxelSizes);
  // Do path tracking for these parameter settings
  PathTracking::PositionType startPoint (sx, sy, sz);
  PathTracking::PositionType endPoint (ex, ey, ez);
  pt.start(startPoint, endPoint);
  std::cout << "Number of path voxels: " << pt.getPath()->size() << std::endl;

  PathTracking::writePathToFile ( outputfile.c_str(), 
                                  pt.getPathPointer(), 
                                  reader->GetOutput() );

  // Free allocated memory
  delete[] data;

  // Exit program
  return EXIT_SUCCESS;
}

