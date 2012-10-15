/*=========================================================================
*
*  Copyright Insight Software Consortium
*
*  Licensed under the Apache License, Version 2.0 (the "License");
*  you may not use this file except in compliance with the License.
*  You may obtain a copy of the License at
*
*         http://www.apache.org/licenses/LICENSE-2.0.txt
*
*  Unless required by applicable law or agreed to in writing, software
*  distributed under the License is distributed on an "AS IS" BASIS,
*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*  See the License for the specific language governing permissions and
*  limitations under the License.
*
*=========================================================================*/
#define BOOST_EXCEPTION_DISABLE 1
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
//#include "itkLevelSetDomainMapImageFilter.h"
//#include "itkLevelSetContainerBase.h"
#include "itkLevelSetEquationCurvatureTerm.h"
#include "itkSinRegularizedHeavisideStepFunction.h"
//#include "itkLevelSetSparseEvolutionBase.h"
//#include "itkBinaryImageToWhitakerSparseLevelSetAdaptor.h"
//#include "itkLevelSetEvolutionNumberOfIterationsStoppingCriterion.h"
#include "itkGeodesicActiveContourAndRegionLevelSetImageFilter.h"
#include "itkGeodesicActiveContourLevelSetImageFilter.h"
#include "itkNumericTraits.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkInvertImageFilter.h"

//#include "itkWhitakerCommandIterationUpdate.h"
#include "boost/program_options.hpp"
#include "boost/filesystem/operations.hpp"
#include <time.h>
#include <math.h>


// ------------------------------------------------------------------------
//
// Exercise: Add a curvature term for regularization.
//
// The coefficient for the new term will be provided as an argument to the
// executable.
//
// ------------------------------------------------------------------------

using namespace std;
// Boost program options
namespace po = boost::program_options;
namespace bfs = boost::filesystem;
int main( int argc, char* argv[] )
{
	std::vector<std::string> parameters;
	std::vector<std::string> help;
	std::vector<float>  parametersF;
	help.push_back(" Geodesic Active Contour with Regional Information");
	help.push_back(" Perform Carotid Segmentation using Geodesic Active Contour with Regional Information");
	help.push_back(" Author: Hui Tang");

	po::options_description desc1("Images");
	desc1.add_options()
		("help,h", "produce help message")
		("initialImage,i", po::value<std::string>(), "inital image, .mhd file ")
		("similarityImage,s", po::value<std::string>(), "similarity image, .mhd file")
		("originalImage,o", po::value<std::string>(), "original image")
		("potentialImage,p", po::value<std::string>(), "potential image")
		//("propergationImage,f", po::value<std::string>(), "propergation image")
		("outputImage,O",po::value<std::string>(),"outputimage");

	po::options_description desc2("Parameters");
	desc2.add_options()
		("parameters,P",po::value<std::vector<std::string> >(), "parameters");

	po::positional_options_description pd; 
	pd.add("parameters", -1);

	po::options_description cmddescAll("Parse");
	cmddescAll.add(desc1).add(desc2);



	po::variables_map vm;
	po::parsed_options parsed=po::command_line_parser(argc, argv).options(cmddescAll).allow_unregistered().positional(pd).run();
	//po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);

	try
	{
		po::store(parsed, vm);
		po::notify(vm);
	}
	catch (std::exception& e)
	{
		cout << "Exception while parsing parameters: " << endl;
		cout << e.what() << endl;
		exit(-1);
	}

	const unsigned int Dimension = 3;

	typedef float                                    InputPixelType;
	typedef itk::Image< InputPixelType, Dimension >           InputImageType;



	// Print help message
	if (vm.count("help") || vm.size()==0) {
		std::cout<<"Author: Hui Tang"<<std::endl;
		std::cout<<"This application is based on Hui Tang, Theo van Walsum Lumen Segmentation for atherosclerotic carotid arteries in CTA, ISBI2012"<< std::endl;
		std::cout<<"You are welcome to use this application if you cite the above paper"<<std::endl;
		std::cout << desc1 << "\n";
		std::cout << "Example:GACRSegmentation.exe -i initial.mhd -p potential.mhd -s similarity.mhd -o original.mhd -O output.mhd ExternalRegionWeight(1) InternalRegionWeight(0.05) AdvectionWeight(20) CurvatureWeight(1-10) iterationTimes(1000) RSM(0.001) minCurvature(1) isovalue(1-2mm)"<< std::endl;

		return EXIT_SUCCESS;
	}
	std::vector<std::string> values(vm["parameters"].as<std::vector<std::string> >()); 
	for(std::vector<std::string>::iterator iter = values.begin(); values.end() != iter; ++iter) 
	{
		parameters.push_back(*iter);
		std::string a = *iter;
		std::istringstream b(a);
		float f;
		b >> f;
		parametersF.push_back(f);
	}
	if (parameters.size()<8)
	{
		std::cout << "Error: not enough parameters provided, check help file" << std::endl;
		return EXIT_FAILURE;
	}

	// String for in- and ouput file
	std::string initialImageN;
	std::string similarityImageN;
	std::string originalImageN;
	std::string potentialImageN;
	std::string outputImageN;
	//std::string propergationImageN;
	// Get intensity file
	if (vm.count("initialImage")) {
		initialImageN = vm["initialImage"].as<std::string>();
		std::cout << "inputImage: " << initialImageN<< std::endl;
	} else {
		std::cout << "Error: no initial image provided" << std::endl;
		return EXIT_FAILURE;
	}
	typedef itk::ImageFileReader< InputImageType >            ReaderType;
	ReaderType::Pointer initialReader = ReaderType::New();
	initialReader->SetFileName( initialImageN );
	initialReader->Update();
	InputImageType::Pointer initial = initialReader->GetOutput();
	std::vector<float> imgExt(3, 0);
	imgExt[0] = initial->GetBufferedRegion().GetSize()[0];
	imgExt[1] = initial->GetBufferedRegion().GetSize()[1];
	imgExt[2] = initial->GetBufferedRegion().GetSize()[2];
	std::cout << "initial Image Extent " << imgExt[0] << ", " << imgExt[1] << ", " << imgExt[2] << std::endl;
	// Get intensity file
	if (vm.count("similarityImage")) {
		similarityImageN = vm["similarityImage"].as<std::string>();
		std::cout << "similarityImage: " << similarityImageN<< std::endl;
	} else {
		std::cout << "Error: no similarity provided" << std::endl;
		return EXIT_FAILURE;
	}
	ReaderType::Pointer similarityReader = ReaderType::New();
	similarityReader->SetFileName( similarityImageN );
	similarityReader->Update();
	InputImageType::Pointer similarityI= similarityReader->GetOutput();
	std::vector<float> imgExt1(3, 0);
	imgExt1[0] = similarityI->GetBufferedRegion().GetSize()[0];
	imgExt1[1] = similarityI->GetBufferedRegion().GetSize()[1];
	imgExt1[2] = similarityI->GetBufferedRegion().GetSize()[2];
	std::cout << "similarity image Extent " << imgExt1[0] << ", " << imgExt1[1] << ", " << imgExt1[2] << std::endl;
	if (imgExt1[0]!=imgExt[0]||imgExt1[1]!=imgExt[1]||imgExt1[2]!=imgExt[2])
	{
		std::cout << "image size should be the same!" << std::endl;
		return EXIT_FAILURE;
	}

	// Get intensity file
	if (vm.count("originalImage")) {
		originalImageN = vm["originalImage"].as<std::string>();
		std::cout << "originalImage: " << originalImageN<< std::endl;
	} else {
		std::cout << "Error: no original image provided" << std::endl;
		return EXIT_FAILURE;
	}
	ReaderType::Pointer originalReader = ReaderType::New();
	originalReader->SetFileName( originalImageN);
	originalReader->Update();

	InputImageType::Pointer original = originalReader->GetOutput();
	std::vector<float> imgExt2(3, 0);
	imgExt2[0] = original->GetBufferedRegion().GetSize()[0];
	imgExt2[1] = original->GetBufferedRegion().GetSize()[1];
	imgExt2[2] = original->GetBufferedRegion().GetSize()[2];
	std::cout << "original image Extent " << imgExt2[0] << ", " << imgExt2[1] << ", " << imgExt2[2] << std::endl;
	if (imgExt[0]!=imgExt2[0]||imgExt[1]!=imgExt2[1]||imgExt[2]!=imgExt2[2])
	{
		std::cout << "image size should be the same!" << std::endl;
		return EXIT_FAILURE;
	}

	InputImageType::Pointer potential= InputImageType::New();
	// Get intensity file
	if (vm.count("potentialImage")) {
		potentialImageN = vm["potentialImage"].as<std::string>();
		std::cout << "potentialImage: " << potentialImageN<< std::endl;
		ReaderType::Pointer potentialReader = ReaderType::New();
		potentialReader->SetFileName( potentialImageN );
		potentialReader->Update();
		potential = potentialReader->GetOutput();
		std::vector<float> imgExt3(3, 0);
		imgExt3[0] = potential->GetBufferedRegion().GetSize()[0];
		imgExt3[1] = potential->GetBufferedRegion().GetSize()[1];
		imgExt3[2] = potential->GetBufferedRegion().GetSize()[2];
		std::cout << "potential image Extent " << imgExt3[0] << ", " << imgExt3[1] << ", " << imgExt3[2] << std::endl;
		if (imgExt2[0]!=imgExt3[0]||imgExt2[1]!=imgExt3[1]||imgExt2[2]!=imgExt3[2])
		{
			std::cout << "image size should be the same!" << std::endl;
			return EXIT_FAILURE;
		}
	} else {
		std::cout << "calculate potential image from original image using 1/(1+|grad(I)|)" << std::endl;
		typedef itk::GradientMagnitudeRecursiveGaussianImageFilter< InputImageType, InputImageType > DerivativeFilterType;
		DerivativeFilterType::Pointer derivative = DerivativeFilterType::New();
		derivative->SetInput( original );
		derivative->SetSigma(0.5);
		derivative->Update();
		typedef itk::InvertImageFilter< InputImageType, InputImageType > InvertFilterType;
		InvertFilterType::Pointer inverter= InvertFilterType::New();
		inverter->SetInput(derivative->GetOutput());
		inverter->SetExponential(1);
		inverter->SetConstant(0.01);
		inverter->Update();
		potential =  inverter->GetOutput();

	}
	/* if (vm.count("propergationImage")) {
	propergationImageN = vm["propergationImage"].as<std::string>();
	std::cout << "propergationImage: " << propergationImageN<< std::endl;
	} else {
	std::cout << "Error: no propergationImage provided" << std::endl;
	return EXIT_FAILURE;
	}
	ReaderType::Pointer propergationReader = ReaderType::New();
	propergationReader->SetFileName( propergationImageN );
	propergationReader->Update();
	InputImageType::Pointer propergationI= propergationReader->GetOutput();*/

	if (vm.count("outputImage")) {
		outputImageN = vm["outputImage"].as<std::string>();
		std::cout << "outputImage: " <<outputImageN<< std::endl;
	} else {
		std::cout << "Error: no outputImage image name provided" << std::endl;
		return EXIT_FAILURE;
	}
	std::cout<<"original"<<original->GetOrigin()[0]<<original->GetOrigin()[1]<<original->GetOrigin()[2]<<std::endl;


	std::cout<<"potential before"<<potential->GetOrigin()[0]<<potential->GetOrigin()[1]<<potential->GetOrigin()[2]<<std::endl;
	potential->SetOrigin(original->GetOrigin() );
	std::cout<<"potential after"<<potential->GetOrigin()[0]<<potential->GetOrigin()[1]<<potential->GetOrigin()[2]<<std::endl;
	std::cout<<"similarityI before"<<similarityI->GetOrigin()[0]<<similarityI->GetOrigin()[1]<<similarityI->GetOrigin()[2]<<std::endl;
	similarityI->SetOrigin( original->GetOrigin()  );
	std::cout<<"similarityI after"<<similarityI->GetOrigin()[0]<<similarityI->GetOrigin()[1]<<similarityI->GetOrigin()[2]<<std::endl;
	std::cout<<"initial before"<<initial->GetOrigin()[0]<<initial->GetOrigin()[1]<<initial->GetOrigin()[2]<<std::endl;
	initial->SetOrigin( original->GetOrigin() );
	std::cout<<"initial after"<<initial->GetOrigin()[0]<<initial->GetOrigin()[1]<<initial->GetOrigin()[2]<<std::endl;




	// Image Dimension


	// Read input image (to be processed).








	typedef itk::GeodesicActiveContourAndRegionLevelSetImageFilter<InputImageType,InputImageType,InputImageType> GACRLevelSetType;
	GACRLevelSetType::Pointer gacrLevelSet=GACRLevelSetType::New();
	gacrLevelSet->SetInput(initial);
	//gacrLevelSet->SetsimilarityImage(similarityI);
	gacrLevelSet->SetFeatureImage(potential);
	gacrLevelSet->SetIntensityImage(original);
	gacrLevelSet->SetSimilarityImage(similarityI);
	gacrLevelSet->SetDerivativeSigma(1);
	gacrLevelSet->SetExternalRegionWeight(parametersF[1]);
	std::cout<<"ExternalRegionWeight:"<<parametersF[1]<<std::endl;
	gacrLevelSet->SetInternalRegionWeight(parametersF[2]);
	std::cout<<"InternalRegionWeight:"<<parametersF[2]<<std::endl;
	gacrLevelSet->SetHeavisideEpsilon(0.5);
	gacrLevelSet->SetAdvectionScaling(parametersF[3]);
	std::cout<<"AdvectionScaling:"<<parametersF[3]<<std::endl;
	gacrLevelSet->SetCurvatureScaling(parametersF[4]);
	std:: cout<<"CurvatureScaling:"<<parametersF[4]<<std::endl;
	gacrLevelSet->SetInterpolateSurfaceLocation(true);
	gacrLevelSet->SetNumberOfIterations(parametersF[5]);
	std::cout<<"NumberOfIterations:"<<parametersF[5]<<std::endl;
	gacrLevelSet->SetMaximumRMSError(parametersF[6]);
	std:: cout<<"MaximumRMSError:"<<parametersF[6]<<std::endl;
	gacrLevelSet->SetUseMinimalCurvature((bool)parametersF[7]);
	std:: cout<<"UseMinimalCurvature"<<parametersF[7]<<std::endl;
	gacrLevelSet->SetIsoSurfaceValue(parametersF[8]);
	std::cout<<"IsoSurfaceValue:"<<parametersF[8]<<std::endl; 
	gacrLevelSet->SetMaximumCurvatureTimeStep(0.0166667);
	gacrLevelSet->SetMaximumPropagationTimeStep(0.0166667);
	gacrLevelSet->Update();


	typedef itk::Image< float, Dimension > OutputImageType;
	OutputImageType::Pointer outputImage = OutputImageType::New();
	outputImage->SetRegions( initial->GetLargestPossibleRegion() );
	outputImage->SetOrigin( original->GetOrigin() );
	outputImage->Allocate();
	outputImage->FillBuffer( 0 );

	typedef itk::ImageRegionIteratorWithIndex< OutputImageType > OutputIteratorType;
	OutputIteratorType oIt( outputImage, outputImage->GetLargestPossibleRegion() );
	oIt.GoToBegin();

	OutputImageType::IndexType idx;

	while( !oIt.IsAtEnd() )
	{
		idx = oIt.GetIndex();
		oIt.Set( gacrLevelSet->GetOutput()->GetPixel(idx) );
		++oIt;
	}

	std::cout << "writing output image file..." << std::endl;
	typedef itk::ImageFileWriter< InputImageType >     OutputWriterType;
	OutputWriterType::Pointer writer = OutputWriterType::New();
	writer->SetFileName(outputImageN);
	writer->SetInput(gacrLevelSet->GetOutput());
	writer->SetUseCompression(true);
	try
	{
		writer->Update();
		std::cout << "Image was written as" << outputImageN<<std::endl;

	}
	catch ( itk::ExceptionObject& err )
	{
		std::cout << err << std::endl;
	}

	return EXIT_SUCCESS;
}
