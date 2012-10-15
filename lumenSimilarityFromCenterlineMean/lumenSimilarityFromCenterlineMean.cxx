/*=========================================================================

Program:   Supress calcium intensity
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
#include "itkNeighborhoodAlgorithm.h"
#include "itkFastMarchingImageFilter.h"
#include "itkConstShapedNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"
#include "itkGaussianDistributionImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "../utilities/itkSeedPointFileIO.h"
#include <math.h>

// Include boost classes
#include "boost/program_options.hpp"

// Boost program options
namespace po = boost::program_options;

int main (int argc, char ** argv) {

	std::vector<std::string> help;
	help.push_back(" meanIntensityImageFromCenterline");
	help.push_back(" calculate lumen intensity similarity from mean image, mean image is from meanIntensityImageFromCenterline");
	help.push_back(" Author: Hui Tang");

	po::options_description desc("Parameters");
	desc.add_options()
		("help,h", "produce help message")
		("inputImage,i",       po::value<std::string>(),                                  "input image ")
		("meanImage,m",       po::value<std::string>(),                                  "centerline file")	
		("standarddeviation,s", po::value<float>(),                                     "standard deviation")
		("outputImage,o",       po::value<std::string>(),                                  "output image file");


	// Parse command line options
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
	po::notify(vm);

	// Print help message
	if (vm.count("help") || vm.size()==0) {
		std::cout << "calculate lumen intensity similarity from mean image, ean image is from meanIntensityImageFromCenterline." << std::endl << std::endl;
		std::cout << "Author: Hui Tang" << "\n";
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
	std::string meanImageFileName;
	if (vm.count("meanImage")) {
		meanImageFileName = vm["meanImage"].as<std::string>();
		std::cout << "mean image: " << meanImageFileName << std::endl;
	} else {
		std::cout << "Error: no mean image provided" << std::endl;
		return EXIT_FAILURE;
	}


	float standarddeviation = 8100;
	if (vm.count("standarddeviation")) {
		standarddeviation = vm["standarddeviation"].as<float>();
	} 
	std::cout << "standard deviation : " << standarddeviation << std::endl;

	std::string outputImageFileName;
	if (vm.count("outputImage")) {
		outputImageFileName = vm["outputImage"].as<std::string>();
		std::cout << "Lumen similarity is : " << outputImageFileName << std::endl;
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
	ImageType::Pointer original=   inputImageReader->GetOutput();


	std::vector<float> imgExt(3, 0);
	imgExt[0] = original->GetBufferedRegion().GetSize()[0];
	imgExt[1] = original->GetBufferedRegion().GetSize()[1];
	imgExt[2] = original->GetBufferedRegion().GetSize()[2];
	std::cout << "original image Extent " << imgExt[0] << ", " << imgExt[1] << ", " << imgExt[2] << std::endl;



	// Load input image
	ReaderType::Pointer meanImageReader = ReaderType::New();
	std::cout << "Reading input image file..." << std::endl;
	meanImageReader->SetFileName(meanImageFileName); 
	try{
		meanImageReader->Update();
	}
	catch ( itk::ExceptionObject &err){
		std::cout << "Error while reading input file !" << std::endl;
		std::cout << err << std::endl;
		return EXIT_FAILURE;
	}
	ImageType::Pointer mean=   meanImageReader->GetOutput();


	std::vector<float> imgExt1(3, 0);
	imgExt1[0] = mean->GetBufferedRegion().GetSize()[0];
	imgExt1[1] = mean->GetBufferedRegion().GetSize()[1];
	imgExt1[2] = mean->GetBufferedRegion().GetSize()[2];
	std::cout << "mean image Extent " << imgExt1[0] << ", " << imgExt1[1] << ", " << imgExt1[2] << std::endl;
	if (imgExt1[0]!=imgExt[0]||imgExt1[1]!=imgExt[1]||imgExt1[2]!=imgExt[2])
	{
		std::cout << "image size should be the same!" << std::endl;
		return EXIT_FAILURE;
	}



	//here to save memory, I use copy constSpeed to generate the mean image
	typedef itk::ImageRegionIteratorWithIndex< ImageType > InputIteratorType;
	InputIteratorType inputIt( original,original->GetLargestPossibleRegion() );
    InputIteratorType meanIt( mean,mean->GetLargestPossibleRegion() );
	inputIt.GoToBegin();
	while( (!inputIt.IsAtEnd()))
	{

		const PixelType currentInputIntensity  = inputIt.Get();
		const PixelType currentMeanIntensity  = meanIt.Get();
		inputIt.Set(exp(-(currentInputIntensity-currentMeanIntensity)*(currentInputIntensity-currentMeanIntensity)/standarddeviation));


	    ++inputIt;
		++meanIt;


	}

// write output image

typedef  itk::ImageFileWriter< ImageType > WriterType;
WriterType::Pointer outputWriter2 = WriterType::New();
std::cout << "writing output image file..." << std::endl;
outputWriter2->SetFileName(outputImageFileName); 
outputWriter2->SetInput( original );
try{
	outputWriter2->Update();
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
