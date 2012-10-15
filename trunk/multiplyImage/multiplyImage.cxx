/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Language:  C++
  Author:    Hui Tang

=========================================================================*/
#if defined (_MSC_VER)
#pragma warning (disable: 4786)
#endif

// Include ITK classes
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMultiplyImageFilter.h"

//#include "../common/pathio.h"
// Include boost classes
#include "boost/program_options.hpp"
#include "itkImageConstIterator.h"
using namespace std;
// Boost program options
namespace po = boost::program_options;

int main (int argc, char ** argv) {

  
      std::vector<std::string> help;
      help.push_back(" MultiplyImageFilter");
      help.push_back(" Author: Hui Tang");
  
    po::options_description desc("Parameters");
    desc.add_options()
      ("help,h", "produce help message")
      ("inputImage1,m", po::value<std::string>(), "input image1, .mhd file ")
      ("inputImage2,n", po::value<std::string>(), "input image2, .mhd file")
      ("outputImage,o", po::value<std::string>(), "output multiplied image");
  
  
    // Parse command line options
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    po::notify(vm);
  
    // Print help message
    if (vm.count("help") || vm.size()==1) {
      std::cout << "Multiply two images." << std::endl << std::endl;
      std::cout << desc << "\n";
      return EXIT_SUCCESS;
    }
  
    // String for in- and ouput file
    std::string inputImage1;
    std::string inputImage2;
    std::string outputImage;
  
    // Get intensity file
    if (vm.count("inputImage1")) {
      inputImage1 = vm["inputImage1"].as<std::string>();
      std::cout << "inputImage1: " << inputImage1 << std::endl;
    } else {
      std::cout << "Error: no image1 provided" << std::endl;
      return EXIT_FAILURE;
    }
  
  
  
    if (vm.count("inputImage2")) {
  	    inputImage2 = vm["inputImage2"].as<std::string>();
          std::cout << "inputImage2: " << inputImage2 << std::endl;
  	} else {
  	    std::cout << "Error: no image2 provided" << std::endl;
  	    return EXIT_FAILURE;
    }
  
  
    // Get outputfile
    if (vm.count("outputImage")) {
      outputImage = vm["outputImage"].as<std::string>();
      std::cout << "Output Image: " << outputImage << std::endl;
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
    ReaderType::Pointer image1Reader = ReaderType::New();
    std::cout << "Reading input image file..." << std::endl;
    image1Reader->SetFileName(inputImage1); 
    try{
      image1Reader->Update();
    }
    catch ( itk::ExceptionObject &err){
      std::cout << "Error while reading input file 1 !" << std::endl;
      std::cout << err << std::endl;
      return EXIT_FAILURE;
    }
  
    ReaderType::Pointer image2Reader = ReaderType::New();
    image2Reader->SetFileName( inputImage2);
    try {
      image2Reader->Update();
    }
    catch ( itk::ExceptionObject &err){
      std::cout << "Error while reading input file 2 !" << std::endl;
      std::cout << err << std::endl;
      return EXIT_FAILURE;
    }
    std::cout << "Reading finished" << std::endl;
  
    // Check sizes
    const  ImageType::SpacingType & sp = image1Reader->GetOutput()->GetSpacing();
    std::cout << "Spacing: "  << sp << std::endl;
  
    const ImageType::SizeType size1 = image1Reader->GetOutput()->GetLargestPossibleRegion().GetSize();
    std::cout << "image1 extent: " << size1 << std::endl;
  
    const ImageType::SizeType size2 = image2Reader->GetOutput()->GetLargestPossibleRegion().GetSize();
    std::cout << "image2 extent: " << size2 << std::endl;
  
    if ( size1 != size2 )
    {
      std::cout << "Error: two images should have the same size" << std::endl;
      return EXIT_FAILURE;
    }
  
    // Multiply images
    typedef  itk::MultiplyImageFilter<ImageType,ImageType,ImageType> ImageFilterType;
    ImageFilterType::Pointer imageFilter = ImageFilterType::New();
    imageFilter->SetInput1(image1Reader->GetOutput());
    imageFilter->SetInput2(image2Reader->GetOutput());
    imageFilter->Update(); 
  
    typedef  itk::ImageFileWriter< ImageType > WriterType;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( outputImage.c_str());
    writer->SetInput(imageFilter->GetOutput());
    try
    {
      writer->Update();
      std::cout << "The multiplication image of " 
                << inputImage1 <<" and " << inputImage2 
                << " is saved as " << outputImage << std::endl;
    }
    catch ( itk::ExceptionObject &err)
    {
      std::cout << "can not update output image, exception object caught !" << std::endl;
      std::cout << err << std::endl;
      return EXIT_FAILURE;
    }
        
  
   // Exit program
    return EXIT_SUCCESS;
}
