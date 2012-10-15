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
#include <windows.h>
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

//#include "../common/pathio.h"
// Include boost classes
#include "boost/program_options.hpp"
#include "itkImageConstIterator.h"
using namespace std;
// Boost program options
namespace po = boost::program_options;

void medialnessCal(float* dataGradientx,float* dataGradienty,float* dataMask,float* output,float a,float d,float r,float R,float n,const std::vector<float> imgExt,const std::vector<float> &voxelSizes);
int main (int argc, char ** argv) {

    std::vector<std::string> help;
    help.push_back(" medialnessImagefilter");
    help.push_back(" Perform medialness image filter based on Robust Vessel Tree Modeling.");

  po::options_description desc("Parameters");
  desc.add_options()
    ("help,h", "produce help message")
    ("inputImage,i", po::value<std::string>(), "input image, .mhd file ")
    ("mask,m", po::value<std::string>(), "mask, .mhd file")
    ("outputMedialness,o", po::value<std::string>()->default_value("Medialness.dcm"), "output medialness")
    ("darkLumen,d", po::value<float>()->default_value(0), "lumen intensity dark or bright(0/1)")
    ("sigma,s", po::value<float>()->default_value(1),"sigma for gradient calculation, default 1")
    ("numberOfAngles,a", po::value<float>(), "numberOfAngles")
    ("minimumRadius,r", po::value<float>(), "minimumRadius")
    ("maximumRadius,R", po::value<float>(), "maximumRadius")
    ("numberOfRadius,n", po::value<float>(), "numberOfRadius");


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



  if (vm.count("mask")) {
	    maskN = vm["mask"].as<std::string>();
        std::cout << "mask used: " << maskN << std::endl;
	} else {
	    std::cout << "Error: no mask image provided" << std::endl;
	    return EXIT_FAILURE;
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
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  
  // ROI
  ImageType::RegionType desiredRegion;
  std::vector<float> imgExt(3, 0);
  

  // Load files
  ReaderType::Pointer inputImageReader = ReaderType::New();

  std::cout << "Reading input image file..." << std::endl;
  inputImageReader->SetFileName(inputImageN); 
  inputImageReader->Update();
  std::cout << "Reading finished" << std::endl;
  ReaderType::Pointer maskReader = ReaderType::New();

  std::cout << "Reading mask file..." << std::endl;

  maskReader->SetFileName( maskN);
  maskReader->Update();
  std::cout << "Reading finished" << std::endl;

  //calculate gradient vector
  typedef float ComponentType;
  typedef itk::CovariantVector< ComponentType, Dimension > OutGradientVectorPixelType;
  typedef itk::Image <OutGradientVectorPixelType,3> OutGradientVectorImageType;
  typedef itk::GradientRecursiveGaussianImageFilter<ImageType,OutGradientVectorImageType> GradientFilterType;
  GradientFilterType::Pointer gradientFilter= GradientFilterType::New();
  gradientFilter->SetSigma(floor(vm["sigma"].as<float>()));
  gradientFilter->SetInput(inputImageReader->GetOutput());
  gradientFilter->Update();
  cout<<"Gradient Calculation Done"<<endl;
  OutGradientVectorImageType::Pointer vectorImage = gradientFilter->GetOutput();
  typedef itk::VectorIndexSelectionCastImageFilter<OutGradientVectorImageType, ImageType> IndexSelectionFilterType;
  IndexSelectionFilterType::Pointer indexSelectionFilterX = IndexSelectionFilterType::New();
  indexSelectionFilterX->SetInput(vectorImage);
  indexSelectionFilterX->SetIndex(0);
  indexSelectionFilterX->Update();
  ImageType::Pointer gradientx=indexSelectionFilterX->GetOutput();
  IndexSelectionFilterType::Pointer indexSelectionFilterY = IndexSelectionFilterType::New();
  indexSelectionFilterY->SetInput(vectorImage);
  indexSelectionFilterY->SetIndex(1);
  indexSelectionFilterY->Update();
  ImageType::Pointer gradienty=indexSelectionFilterY->GetOutput();

  //std::cout << gradientFilter->GetOutput()->GetPixel(index) << std::endl;
  //std::cout << adaptor->GetPixel(index) << std::endl;
  //GradientVectorImageType::Pointer gradient = GradientVectorImageType::New();
  /*typedef itk::ExtractImageFilter< ImageType, ImageType > FilterType;
  FilterType::Pointer filter0 = FilterType::New();
  FilterType::Pointer filter1 = FilterType::New();
  FilterType::Pointer filter2 = FilterType::New();
  typedef itk::ImageIOBase                        IOBaseType;
  ImageType::IndexType start;
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  desiredRegion.SetIndex( start );
  ImageType::SizeType size = inputImageReader->GetOutput()->GetBufferedRegion().GetSize();
  desiredRegion.SetSize( size );


  filter0->SetExtractionRegion( desiredRegion );
  filter0->SetInput( indexSelectionFilterX->GetOutput() );
  filter0->Update();
  filter0->GetOutput()->GetPixelContainer()->SetContainerManageMemory(false);

  filter1->SetExtractionRegion( desiredRegion );
  filter1->SetInput( indexSelectionFilterY->GetOutput() );
  filter1->Update();
  filter1->GetOutput()->GetPixelContainer()->SetContainerManageMemory(false);

  filter2->SetExtractionRegion( desiredRegion );
  filter2->SetInput( maskReader->GetOutput() );
  filter2->Update();
  filter2->GetOutput()->GetPixelContainer()->SetContainerManageMemory(false);*/

  //Set pointer to data

  imgExt[0] = inputImageReader->GetOutput()->GetBufferedRegion().GetSize()[0];
  imgExt[1] = inputImageReader->GetOutput()->GetBufferedRegion().GetSize()[1];
  imgExt[2] = inputImageReader->GetOutput()->GetBufferedRegion().GetSize()[2];
  float* dataGradientX = new float[imgExt[0],imgExt[1],imgExt[2]];
  float* dataGradientY = new float[imgExt[0],imgExt[1],imgExt[2]];
  float* dataMask =new float[imgExt[0],imgExt[1],imgExt[2]];
  float* output = new float[imgExt[0],imgExt[1],imgExt[2]];


  /*dataGradientX= filter0->GetOutput()->GetPixelContainer()->GetImportPointer();
  dataGradientY= filter1->GetOutput()->GetPixelContainer()->GetImportPointer();
  dataMask = filter2->GetOutput()->GetPixelContainer()->GetImportPointer();*/


  typedef itk::ImageRegionIteratorWithIndex<ImageType> ImageIterator;
  ImageIterator it0( gradientx, gradientx->GetRequestedRegion() );
  for (it0.GoToBegin(); !it0.IsAtEnd(); ++it0)
  { 
      ImageType::IndexType idx=it0.GetIndex();
       dataGradientX[idx[0],idx[1],idx[2]]=it0.Get();
       cout<<idx[0]<<idx[1]<<idx[2]<<"  "<<it0.Get()<<"  ";
       idx[0]=0;
       idx[1]=0;
       idx[2]=0;
      it0.SetIndex(idx);
       cout<<idx[0]<<idx[1]<<idx[2]<<"  "<<it0.Get()<<"  ";
       //cout<<idx[0]<<idx[1]<<idx[2]<<dataGradientX[idx[0],idx[1],idx[2]]<<"  "<<it0.Get()<<"  ";

  }
  ImageIterator it1( gradienty, gradienty->GetRequestedRegion() );
  for (it1.GoToBegin(); !it1.IsAtEnd(); ++it1)
  { 
      ImageType::IndexType idx = it1.GetIndex();
      dataGradientY[idx[0],idx[1],idx[2]]=it1.Get();

  }
  ImageIterator it2( maskReader->GetOutput(),maskReader->GetOutput()->GetRequestedRegion() );
  for (it2.GoToBegin(); !it2.IsAtEnd(); ++it2)
  { 
      ImageType::IndexType idx = it2.GetIndex();
      dataMask[idx[0],idx[1],idx[2]]=it2.Get();

  }

  /*dataGradientX= gradientx->GetBufferPointer();
  dataGradientY= gradienty->GetBufferPointer();
  dataMask = maskReader->GetOutput()->GetBufferPointer();*/
  for (int i=0;i<512;i++)
  {
      for (int j=0;j<512;j++)
      {
           cout<<dataGradientX[i,j,0]<<endl;
      }
  }
  cout<<dataGradientX[515,515,0]<<endl;
  cout<<dataGradientX[1,1,0]<<endl;
  cout<<dataGradientX[2,3,0]<<endl;
  cout<<dataGradientX[4,5,0]<<endl;
  cout<<dataGradientY[511,511,0]<<endl;
  cout<<dataMask[511,511,0]<<endl;





  std::cout << "Image extent: " << imgExt[0] << ", " << imgExt[1] << ", " << imgExt[2] << std::endl;

  // Get voxelsizes
  const ImageType::SpacingType& sp = inputImageReader->GetOutput()->GetSpacing();
  std::cout << "Spacing: "  << sp[0] << ", " << sp[1] << ", " << sp[2] << std::endl;
  std::vector<float> voxelSizes (3, 1.0f);
  voxelSizes[0] = sp[0];
  voxelSizes[1] = sp[1];
  voxelSizes[2] = sp[2];

  // Floats for start- and endpoint

  if (vm.count("darkLumen") * vm.count("numberOfAngles") * vm.count("minimumRadius")*vm.count("maximumRadius")*vm.count("numberOfRadius")&&!dataGradientX==NULL&&!dataGradientY==NULL&&!dataMask==NULL) {
      float a;
      float r;
      float R;
      float n;
      float d;
      d = floor(vm["darkLumen"].as<float>());
      a = floor(vm["numberOfAngles"].as<float>());
      r = floor(vm["minimumRadius"].as<float>());
      R = floor(vm["maximumRadius"].as<float>());
      n = floor(vm["numberOfRadius"].as<float>());
      medialnessCal(dataGradientX,dataGradientY,dataMask,output,a,d,r,R,n,imgExt,voxelSizes);
      typedef itk::ImageFileWriter< ImageType > WriterType;
      WriterType::Pointer writer = WriterType::New();
      ImageType:: Pointer outputimage = ImageType::New();
      outputimage->SetRegions( inputImageReader->GetOutput()->GetLargestPossibleRegion() );
      outputimage->Allocate();
      outputimage->SetSpacing( inputImageReader->GetOutput()->GetSpacing() );
      outputimage->SetOrigin( inputImageReader->GetOutput()->GetOrigin() );


      typedef itk::ImageRegionIteratorWithIndex<ImageType> ImageIterator;
      ImageIterator it( outputimage, outputimage->GetRequestedRegion() );

      for (it.GoToBegin(); !it.IsAtEnd(); ++it)
      { 
          ImageType::IndexType idx = it.GetIndex();
          if ( !output==NULL)
          {
               it.Set( output[idx[0],idx[1],idx[2] ]);
          } 
          else
          {
               it.Set(0);
          }
         
      }

      writer->SetFileName( outputMedialnessN.c_str());
      writer->SetInput(outputimage);
      try
      {
          writer->Update();
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
      typedef itk::ImageFileWriter< ImageType > WriterType;
      WriterType::Pointer writer = WriterType::New();
      ImageType:: Pointer outputimage = ImageType::New();
      outputimage->SetRegions( inputImageReader->GetOutput()->GetLargestPossibleRegion() );
      outputimage->Allocate();
      outputimage->SetSpacing( inputImageReader->GetOutput()->GetSpacing() );
      outputimage->SetOrigin( inputImageReader->GetOutput()->GetOrigin() );


      typedef itk::ImageRegionIteratorWithIndex<ImageType> ImageIterator;
      ImageIterator it( outputimage, outputimage->GetRequestedRegion() );

      for (it.GoToBegin(); !it.IsAtEnd(); ++it)
      { 
          ImageType::IndexType idx = it.GetIndex();
          it.Set(0);

      }

      writer->SetFileName( outputMedialnessN.c_str());
      writer->SetInput(outputimage);
      try
      {
          writer->Update();
      }
      catch ( itk::ExceptionObject &err)
      {
          std::cout << "can not update output image, exception object caught !" << std::endl;
          std::cout << err << std::endl;
          return -1;
      }
  }


  
  
  // Free allocated memory
  delete[] dataGradientX;
  delete[] dataGradientY;
  delete[] dataMask;
  delete[] output;


  // Exit program
  return EXIT_SUCCESS;
}
 


void medialnessCal (float* dataGradientx,float* dataGradienty,float* dataMask,float* output,float a,float d,float r,float R,float n,const std::vector<float> imgExt,const std::vector<float> &voxelSizes)
{
    for (int zi = 0; zi<imgExt[2];  zi++){
        for (int yi = 0;  yi<imgExt[1];  yi++){
            for (int xi = 0; xi<imgExt[0];  ++xi){   
                if (dataMask[xi,yi,zi])
                {
                    std::vector<std::vector<float>> allBoundaryMeasure;
                    std::vector<std::vector<float>> allIntensityProfile;
                    std::vector<float> maxBoundaryMeasure(a);
                    std::vector<int> maxBoundaryMeasureIndex(a);
                    float stepOfAngle =2*3.1415926/a;
                    for (int i=0; i<a;i++)
                    {
                        std::vector<float> boundaryMeasureAtiAngle;
                        float deltaRadius= (R-r)/(n-1);
                        for (int ri=0;ri<n;ri++)                                      
                        {    
                            //vec3 currentSV(x,y,z);
                            //vec3 currentSW;
                            //getInImg(0)->transformToWorldCoord(currentSV,currentSW);
                            //float xW = currentSW[0] + (r+r*deltaRadius)*cos(i*stepOfAngle);
                            //float yW = currentSW[1] + (r+r*deltaRadius)*sin(i*stepOfAngle);
                            //float zW = currentSW[2];  

                            //vec3 currentRW(xW,yW,zW);
                            //vec3 currentRV;
                            //getInImg(0)->transformToVoxelCoord(currentRW,currentRV);

                            int x=min(max(xi+((ri+ri*deltaRadius)*cos(i*stepOfAngle)/voxelSizes[0]),(float)0),(float)(imgExt[0]-1));
                            int y=min(max(yi+((ri+ri*deltaRadius)*sin(i*stepOfAngle)/voxelSizes[1]),(float)0),(float)(imgExt[1]-1));
                            int z=zi;
                            x=(int)x;
                            y=(int)y;
                            z=(int)z;

                            if (d)
                            {
                                boundaryMeasureAtiAngle.push_back((float)max(dataGradientx[x,y,z]*cos(i*stepOfAngle)+dataGradienty[x,y,z]*sin(i*stepOfAngle),(float)0));

                            }
                            else
                            { 
                                boundaryMeasureAtiAngle.push_back((float)min(dataGradientx[x,y,z]*cos(i*stepOfAngle)+dataGradienty[x,y,z]*sin(i*stepOfAngle),(float)0));
                            }

                            if (abs((float)(dataGradientx[x,y,z]*cos(i*stepOfAngle)+dataGradienty[x,y,z]*sin(i*stepOfAngle)))>maxBoundaryMeasure[i])
                            {
                                maxBoundaryMeasure[i]=abs((float)(dataGradientx[x,y,z]*cos(i*stepOfAngle)+dataGradienty[x,y,z]*sin(i*stepOfAngle)));
                                maxBoundaryMeasureIndex[i]=ri;
                            }


                        }
                        allBoundaryMeasure.push_back(boundaryMeasureAtiAngle);

                    }

                    std::vector<std::vector<float> > allBoundaryMeasureN;     
                    for (int i=0; i<a;i++)
                    {
                        std::vector<float> boundaryMeasureAtiAngleN;
                        for (int r=0;r<n;r++)                                      
                        {  
                            std::vector<float> currentBoundaryMeasureFrom0ToR;
                            float currentMax=maxBoundaryMeasure[i];
                            if (abs(currentMax)>0.005)
                            {
                                boundaryMeasureAtiAngleN.push_back((float)abs((float)allBoundaryMeasure[i][r])/abs(currentMax));
                                if ((float)abs((float)allBoundaryMeasure[i][r])/abs(currentMax)>1)
                                {
                                    std::cout<<abs(currentMax)<<std::endl;
                                }
                            }
                            else
                            {
                                boundaryMeasureAtiAngleN.push_back((float)0);
                            }

                        }

                        allBoundaryMeasureN.push_back(boundaryMeasureAtiAngleN);
                    }
                    float maxMedialness=0;
                    for (int ri=0;ri<n;r++)                                      
                    {   
                        float currentMedialness=0;
                        for (int i=0; i<a;i++)
                        {  
                            currentMedialness=currentMedialness+allBoundaryMeasureN[i][r];
                        }
                        if (currentMedialness>maxMedialness)
                        {
                            maxMedialness=currentMedialness;
                        }
                    }	
                    output[xi,yi,zi]=(float)maxMedialness/a;
                }
            }
        }
    }
} 
