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
#include "itkExtractImageFilter.h"

// Include path classes
#include "pathTrackingLibrary.cpp"
#include "../common/pathio.h"

// Include boost classes
#include "boost/program_options.hpp"

// Boost program options
namespace po = boost::program_options;

// Forward declaration
void doPathTracking (float sx, float sy, float sz, float ex, float ey, float ez,
                     std::string outputfile, const std::vector<int> imgExt,
                     float* &data, const std::vector<float> &voxelSizes);

int main (int argc, char ** argv) {
  // Directories
  std::string inputdir, outputdir;

  // Declare the supported options.
  po::options_description desc("Parameters");
  desc.add_options()
    ("help,h", "produce help message")
    ("costimage,i", po::value<std::string>(), "cost image filename")
    ("outputfile,o", po::value<std::string>()->default_value("path.txt"), "output file")
    ("startx,x", po::value<float>(), "start point x-coordinate\n(in voxel coordinates)")
    ("starty,y", po::value<float>(), "start point y-coordinate\n(in voxel coordinates)")
    ("startz,z", po::value<float>(), "start point z-coordinate\n(in voxel coordinates)")
    ("endx,X", po::value<float>(), "end point x-coordinate\n(in voxel coordinates)")
    ("endy,Y", po::value<float>(), "end point y-coordinate\n(in voxel coordinates)")
    ("endz,Z", po::value<float>(), "end point z-coordinate\n(in voxel coordinates)");

  // Parse command line options
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).
          options(desc).run(), vm);
  po::notify(vm);

  // Output program header
  programHeader("PathTracking", "Coert Metz", 2008);
  
  // Print help message
  if (vm.count("help") || vm.size()==1) {
    std::cout << "Perform pathtracking using cost image." << std::endl << std::endl;
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

  // Floats for start- and endpoint
  float sx, sy, sz, ex, ey, ez;
  
  // Get startpoint
  if (vm.count("startx") * vm.count("starty") * vm.count("startz")) {
    sx = floor(vm["startx"].as<float>()); // /voxelSizes[0]);
    sy = floor(vm["starty"].as<float>()); // /voxelSizes[1]);
    sz = floor(vm["startz"].as<float>()); // /voxelSizes[2]);
    std::cout << "Start point (in voxels): " << sx << ", " << sy << ", " << sz << std::endl;
  } else {
    std::cout << "Start coordinate missing" << std::endl;
    return EXIT_FAILURE;
  }

  // Get endpoint
  if (vm.count("endx") * vm.count("endy") * vm.count("endz")) {
    ex = floor(vm["endx"].as<float>()); // /voxelSizes[0]);
    ey = floor(vm["endy"].as<float>()); // /voxelSizes[1]);
    ez = floor(vm["endz"].as<float>()); // /voxelSizes[2]);
    std::cout << "End point (in voxels): " << ex << ", " << ey << ", " << ez << std::endl;
  } else {
    std::cout << "End coordinate missing" << std::endl;
    return EXIT_FAILURE;
  }

  doPathTracking(sx, sy, sz, ex, ey, ez, outputfile, imgExt, data, voxelSizes);

  // Free allocated memory
  delete[] data;

  // Exit program
  return EXIT_SUCCESS;
}

void doPathTracking (float sx, float sy, float sz, float ex, float ey, float ez, 
                     std::string outputfile, const std::vector<int> imgExt, 
                     float* &data, const std::vector<float> &voxelSizes) {
  PathTracking pt (imgExt, data, voxelSizes);
  std::cout << "Perform path tracking..." << std::endl;
  // Do path tracking for these parameter settings
  position start (sx, sy, sz);
  position end (ex, ey, ez);
  pt.start(start, end);
  std::cout << "Number of path voxels: " << pt.getPath()->size() << std::endl;
  writePathToFile (outputfile.c_str(), pt.getPathPointer(), false, true, voxelSizes[0], voxelSizes[1], voxelSizes[2]);
}
