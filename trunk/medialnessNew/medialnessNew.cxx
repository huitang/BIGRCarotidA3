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
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include "vtksys/CommandLineArguments.hxx"
#include "vtkMetaImageReader.h"
#include "vtkImageData.h"
#include "vtkMetaImageWriter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkNthElementImageAdaptor.h"
#include "itkIndex.h"
#include "itkImageAdaptor.h"
#include "itkImageRegionIterator.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkMedialnessImageFilterNew.h"

//#include "../common/pathio.h"
// Include boost classes
#include "boost/program_options.hpp"
#include "itkImageConstIterator.h"
using namespace std;
// Boost program options
namespace po = boost::program_options;

int main (int argc, char ** argv) {

  std::vector<std::string> help;
  help.push_back(" medialnessImagefilter");
  help.push_back(" Perform slice based medialness image filter based on Robust Vessel Tree Modeling.");
  help.push_back(" Author: Hui Tang");

  po::options_description desc("Parameters");
  desc.add_options()
    ("help,h", "produce help message")
    ("inputImage,i",       po::value<std::string>(),                                  "input image, .mhd file ")
    ("mask,m",             po::value<std::string>(),                                  "mask, .mhd file")
    ("outputMedialness,o", po::value<std::string>()->default_value("Medialness.dcm"), "output medialness")
    ("darkLumen,d",        po::value<float>()->default_value(0),                      "lumen intensity dark or bright(0/1),2 if neither")
    ("sigma,s",            po::value<float>()->default_value(1),                      "sigma for gradient calculation, default 1")
    ("threshold,t",        po::value<float>()->default_value(0),                      "threshould of mask image")
    ("numberOfAngles,a",   po::value<int>(),                                          "numberOfAngles")
    ("minimumRadius,r",    po::value<float>(),                                        "minimumRadius")
    ("maximumRadius,R",    po::value<float>(),                                        "maximumRadius")
    ("numberOfRadius,n",   po::value<int>(),                                          "numberOfRadius")
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

  // Get intensity file
  int threads = 1;
  if (vm.count("numberOfThreads")) {
    threads = vm["numberOfThreads"].as<int>();
  } 
  std::cout << "numberOfThreads: " << threads << std::endl;
  itk::MultiThreader::Pointer mt = itk::MultiThreader::New();
  mt->SetGlobalMaximumNumberOfThreads( threads );

  // String for in- and ouput file
  std::string inputImageN;
  std::string maskN;
  std::string outputMedialnessN;

  // Get intensity file
  if (vm.count("inputImage")) {
    inputImageN = vm["inputImage"].as<std::string>();
    std::cout << "inputImage: " << inputImageN << std::endl;
  } else {
    std::cout << "Error: no image provided" << std::endl;
    return EXIT_FAILURE;
  }



  bool useMask = false;
  if (vm.count("mask")) {
    maskN = vm["mask"].as<std::string>();
    std::cout << "Mask used: " << maskN << std::endl;
    useMask = true;
	} else {
    std::cout << "No mask image provided" << std::endl;
  }


  // Get outputfile
  if (vm.count("outputMedialness")) {
    outputMedialnessN = vm["outputMedialness"].as<std::string>();
    std::cout << "Output medialness: " << outputMedialnessN << std::endl;
  } else {
    std::cout << "Error: no output file provided" << std::endl;
    return EXIT_FAILURE;
  }

  // Define pixel type, dimension, image type and reader type
  typedef float PixelType;
  const unsigned int Dimension = 3;
  typedef  itk::Image< PixelType, Dimension > ImageType;
  typedef  itk::ImageFileReader< ImageType > ReaderType;
  
  // ROI
  ImageType::RegionType desiredRegion;
  

  // Load files
  ReaderType::Pointer inputImageReader = ReaderType::New();

  std::cout << "Reading input image file..." << std::endl;
  inputImageReader->SetFileName(inputImageN); 
  inputImageReader->Update();
  std::cout << "Reading finished" << std::endl;

  // Get voxelsizes
  const  ImageType::SpacingType & sp = inputImageReader->GetOutput()->GetSpacing();
  std::cout << "Spacing: "  << sp[0] << ", " << sp[1] << ", " << sp[2] << std::endl;

  std::vector<float> inputExtent(3, 0);
  inputExtent[0] = inputImageReader->GetOutput()->GetBufferedRegion().GetSize()[0];
  inputExtent[1] = inputImageReader->GetOutput()->GetBufferedRegion().GetSize()[1];
  inputExtent[2] = inputImageReader->GetOutput()->GetBufferedRegion().GetSize()[2];
  std::cout << "Orginal image extent: " << inputExtent[0] << ", " << inputExtent[1] << ", " << inputExtent[2] << std::endl;

  ReaderType::Pointer maskImageReader = ReaderType::New();
  if ( useMask ) {
    std::cout << "Reading mask file..." << std::endl;
    maskImageReader->SetFileName( maskN);
    maskImageReader->Update();
    std::cout << "Reading finished" << std::endl;

    std::vector<float> maskExtent(3, 0);
    maskExtent[0] = maskImageReader->GetOutput()->GetBufferedRegion().GetSize()[0];
    maskExtent[1] = maskImageReader->GetOutput()->GetBufferedRegion().GetSize()[1];
    maskExtent[2] = maskImageReader->GetOutput()->GetBufferedRegion().GetSize()[2];
    std::cout << "Mask image extent: " << maskExtent[0] << ", " << maskExtent[1] << ", " << maskExtent[2] << std::endl;

    bool equalSize = false;
    if ( (inputExtent[0]==maskExtent[0]) &&
      (inputExtent[1]==maskExtent[1]) &&
      (inputExtent[2]==maskExtent[2]) )
    {
      equalSize = true;
    }
    else
    {
      std::cout<<"Error:original and mask image should have the same size"<<endl;
      return EXIT_FAILURE;
    }
  }

  if (  vm.count("threshold")* 
        vm.count("sigma") * 
        vm.count("darkLumen") * 
        vm.count("numberOfAngles") * 
        vm.count("minimumRadius") *
        vm.count("maximumRadius") *
        vm.count("numberOfRadius")) {

      int   a = vm["numberOfAngles"].as<int>();
      std::cout << "numberOfAngles: " << a << std::endl;
      const float r = floor(vm["minimumRadius"].as<float>());
      std::cout << "minimumRadius: " << r << std::endl;
      const float R = floor(vm["maximumRadius"].as<float>());
      std::cout << "maximumRadius: " << R << std::endl;
      const int   n = vm["numberOfRadius"].as<int>();
      std::cout << "numberOfRadius: " << n << std::endl;
      const float d = floor(vm["darkLumen"].as<float>());
      std::cout << "darkLumen: " << d << std::endl;
      const float s = floor(vm["sigma"].as<float>());
      std::cout << "sigma: " << s << std::endl;
      const float t = floor(vm["threshold"].as<float>());
      std::cout << "threshold: " << t << std::endl;

      typedef  itk::MedialnessImageFilter<ImageType,ImageType,ImageType> MedialnessImageFilterType;
      MedialnessImageFilterType::Pointer medialnessImageFilter = MedialnessImageFilterType::New();
      medialnessImageFilter->SetInput(inputImageReader->GetOutput());
      if ( useMask ) {
        medialnessImageFilter->SetUseMask( true );
        medialnessImageFilter->SetMaskImage(maskImageReader->GetOutput());
      } else {
        medialnessImageFilter->SetUseMask( false );
      }
      medialnessImageFilter->SetNumberOfAngles(a);
      medialnessImageFilter->SetDarkLumen(d);
      medialnessImageFilter->SetMinimalRadius(r);
      medialnessImageFilter->SetMaximalRadius(R);
      medialnessImageFilter->SetThreshold(t);
      medialnessImageFilter->SetNumberOfRadii(n);
      medialnessImageFilter->SetSigma(s);

      std::cout << "medialness calculating...." << std::endl;
      try 
      {
        medialnessImageFilter->Update(); // If you are reusing a pipeline you should call the method UpdateLargestPossibleRegion() in the last filter of the pipeline, instead of just using: Update() Otherwise, the pipeline expects that the new input image has the same size of the previous input image. 
      }
      catch ( itk::ExceptionObject &err)
      {
        std::cout << "error in medialness calculation !" << std::endl;
        std::cout << err << std::endl;
        std::cin.get();
        return EXIT_FAILURE;
      }

      typedef  itk::ImageFileWriter< ImageType > WriterType;
      WriterType::Pointer writer = WriterType::New();
      writer->SetFileName( outputMedialnessN.c_str());
      writer->SetInput(medialnessImageFilter->GetOutput());
      try
      {
          writer->Update();
          cout<<"The medialness image of "<<vm["inputImage"].as<std::string>()<<" is saved as a name of" << vm["outputMedialness"].as<std::string>()<<endl;
      }
      catch ( itk::ExceptionObject &err)
      {
          std::cout << "can not update output image, exception object caught !" << std::endl;
          std::cout << err << std::endl;
          return -1;
      }
      
  } else {
      std::cout << "missing some parameter(s) or input image is not valid" << std::endl;
      return EXIT_FAILURE;
      
  }

 // Exit program
  return EXIT_SUCCESS;
}
