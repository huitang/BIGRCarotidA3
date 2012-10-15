//
//

#ifndef __MEDIALNESSIMAGEFILTER__TXX
#define __MEDIALNESSIMAGEFILTER__TXX

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkNthElementImageAdaptor.h"
#include "itkIndex.h"
#include "itkImageAdaptor.h"
#include "itkImageRegionIterator.h"
#include "itkVectorIndexSelectionCastImageFilter.h"

namespace itk {

#define mymin(a,b) ((a) < (b) ? (a) : (b))
#define mymax(a,b) ((a) > (b) ? (a) : (b))

template< class TInputImage, class TMaskImage, class TOutputImage > 
MedialnessImageFilter< TInputImage, TMaskImage, TOutputImage >
::MedialnessImageFilter(){
  m_UseMask = true;
}

template< class TInputImage, class TMaskImage, class TOutputImage > 
void MedialnessImageFilter< TInputImage, TMaskImage, TOutputImage >
::BeforeThreadedGenerateData(){
  typedef float ComponentType;
  typedef typename itk::CovariantVector <ComponentType,InputImageType::ImageDimension> OutGradientVectorPixelType;
  typedef typename itk::Image <OutGradientVectorPixelType,InputImageType::ImageDimension> OutGradientVectorImageType;
  typedef typename itk::GradientRecursiveGaussianImageFilter <InputImageType,OutGradientVectorImageType> GradientFilterType;
  typename GradientFilterType::Pointer gradientFilter =  GradientFilterType::New();
  gradientFilter->SetSigma(this->m_Sigma);
  gradientFilter->SetInput(this->GetInput());
  gradientFilter->UpdateLargestPossibleRegion();
  typename OutGradientVectorImageType::Pointer vectorImage = gradientFilter->GetOutput();

  typedef typename itk::VectorIndexSelectionCastImageFilter <OutGradientVectorImageType, InputImageType> IndexSelectionFilterType;
  typename IndexSelectionFilterType::Pointer indexSelectionFilterX =  IndexSelectionFilterType::New();
  indexSelectionFilterX->SetInput(vectorImage);
  indexSelectionFilterX->SetIndex(0);
  indexSelectionFilterX->UpdateLargestPossibleRegion();
  m_XGradientImage=indexSelectionFilterX->GetOutput();

  typename IndexSelectionFilterType::Pointer indexSelectionFilterY =  IndexSelectionFilterType::New();
  indexSelectionFilterY->SetInput(vectorImage);
  indexSelectionFilterY->SetIndex(1);
  indexSelectionFilterY->UpdateLargestPossibleRegion();
  m_YGradientImage=indexSelectionFilterY->GetOutput();
}

template< class TInputImage, class TMaskImage, class TOutputImage > 
void MedialnessImageFilter< TInputImage, TMaskImage, TOutputImage >
::ThreadedGenerateData(const typename ImageToImageFilter<TInputImage,TOutputImage>::OutputImageRegionType& outputRegionForThread, int threadId)
{

  typename OutputImageType::Pointer output = this->GetOutput();
  typedef typename InputImageType::SizeType InputSizeType;
  const unsigned int Dimension = 3;

  //calculate gradient of input image

  //get image spacing and size
  std::vector<int> imgExt(3, 0);
  imgExt[0] = outputRegionForThread.GetSize()[0];
  imgExt[1] = outputRegionForThread.GetSize()[1];
  imgExt[2] = outputRegionForThread.GetSize()[2];
  itk::MultiThreader::Pointer mt = itk::MultiThreader::New();
  int threads = mt->GetGlobalDefaultNumberOfThreads();

  const float totalSize = static_cast<float>( imgExt[0]*imgExt[1]*imgExt[2] )/threads;

  const typename InputImageType::SpacingType & sp = this->GetOutput()->GetSpacing();
  std::vector<float> voxelSizes (3, 1.0f);
  voxelSizes[0] = sp[0];
  voxelSizes[1] = sp[1];
  voxelSizes[2] = sp[2];

  typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< InputImageType > FaceCalculatorType;
  FaceCalculatorType faceCalculatorX;
  FaceCalculatorType::FaceListType faceListX;
  FaceCalculatorType::FaceListType::iterator fitX;
  faceListX = faceCalculator( m_XGradientImage,
                              outputRegionForThread,
                              m_MaximalRadius );
  FaceCalculatorType faceCalculatorY;
  FaceCalculatorType::FaceListType faceListY;
  FaceCalculatorType::FaceListType::iterator fitY;
  faceListY = faceCalculator( m_YGradientImage,
                              outputRegionForThread,
                              m_MaximalRadius );

  //Region Iterator
  typedef itk::ConstShapedNeighborhoodIterator<InputImageType> ShapedNeighborhoodIteratorType;
  ShapedNeighborhoodIteratorType::RadiusType radius;
  radius.Fill( m_MaximalRadius );


  typedef typename itk::ImageRegionIterator<InputImageType> ImageIterator;


  //calculatemedialness
  const float deltaRadius= (this->m_MaximalRadius-this->m_MinimalRadius)/(this->m_NumberOfRadii-1);
  const float stepOfAngle =2*3.1415926/this->m_NumberOfAngles;
  //std::vector<std::vector<float> > shiftx;
  //std::vector<std::vector<float> > shifty;

  for ( fit=faceList.begin(); fit != faceList.end(); ++fit ){
    ShapedNeighborhoodIteratorType gradientImageXIt(radius, m_XGradientImage, *fit);
    ShapedNeighborhoodIteratorType gradientImageYIt(radius, m_YGradientImage, *fit);
    ImageIterator outputIt(this->GetOutput(),*fit);
    ImageIterator maskIt;
    if ( m_UseMask && this->GetMaskImage() ){
      maskIt = ImageIterator(this->GetMaskImage(),*fit);
      maskIt.GoToBegin();
    }

    // LUT for radial positions
    for (int i=0; i<this->m_NumberOfAngles;i++)
    {
      std::vector<float> shiftxR;
      std::vector<float> shiftyR;
      for (int ri=0;ri<this->m_NumberOfRadii;ri++)                                      
      { 
        float x = (this->m_MinimalRadius+ri*deltaRadius)*cos(i*stepOfAngle)/voxelSizes[0];
        float y = (this->m_MinimalRadius+ri*deltaRadius)*sin(i*stepOfAngle)/voxelSizes[1];
        //shiftxR.push_back((this->m_MinimalRadius+ri*deltaRadius)*cos(i*stepOfAngle)/voxelSizes[0]);
        //shiftyR.push_back((this->m_MinimalRadius+ri*deltaRadius)*sin(i*stepOfAngle)/voxelSizes[1]);

        ShapedNeighborhoodIteratorType::OffsetType off;
        
        off[0] = static_cast<int>(x);
        off[1] = static_cast<int>(y);
        it.ActivateOffset(off);
      }
      //shiftx.push_back(shiftxR);
      //shifty.push_back(shiftyR);
    }

    gradientImageXIt.GoToBegin(); 
    gradientImageYIt.GoToBegin();  
    outputIt.GoToBegin();
    for (; !outputIt.IsAtEnd(); ++gradientImageXIt, gradientImageYIt, ++outputIt ){
      bool useCurrentVoxel = false;
      if (!m_UseMask ) { 
        useCurrentVoxel = true;
      } else {
        ++maskIt;
        if ( maskIt.Get()>m_Threshold ) {
          useCurrentVoxel = true;
        }
      }
      if ( useCurrentVoxel ){
        ShapedNeighborhoodIteratorType::ConstIterator ciX = gradientImageXIt.GoToBegin();
        ShapedNeighborhoodIteratorType::ConstIterator ciY = gradientImageYIt.GoToBegin();
        for (; ciX != gradientImageXIt.End(); ciX++, ciY++ ){
          std::vector<std::vector<float> > allBoundaryMeasure;
          std::vector<std::vector<float> > allIntensityProfile;
          std::vector<float> maxBoundaryMeasure(this->m_NumberOfAngles);
          std::vector<int> maxBoundaryMeasureIndex(this->m_NumberOfAngles);

          float graX = ciX.Get();
          float graY = ciY.Get();
          float gra = graX*cos(i*stepOfAngle)+graY*sin(i*stepOfAngle);// here it is actually <G,n>

          // negtive gradient will not be taken into account, when the vessel voxel is black. so all gradient vectors are outward, so <G,n>is negtive
          if (this->m_DarkLumen==1)
          {
            boundaryMeasureAtiAngle.push_back(mymax(gra,0.0f));

          }
          if (this->m_DarkLumen==0)
          { 
            boundaryMeasureAtiAngle.push_back(mymin(gra,0.0f));
          }
          //for those cases whose lumen intensity is not black or bright but in between
          if (this->m_DarkLumen>1)
          { 
            boundaryMeasureAtiAngle.push_back( fabsf(gra) );
          }


          if ( fabsf(gra)>maxBoundaryMeasure[i] )
          {
            maxBoundaryMeasure[i]=fabsf(gra);
            maxBoundaryMeasureIndex[i]=ri;
          }

          allBoundaryMeasure.push_back(boundaryMeasureAtiAngle);
        } //ciX

        std::vector<std::vector<float> > allBoundaryMeasureN;     
        for (int i=0; i<this->m_NumberOfAngles;i++)
        {
          std::vector<float> boundaryMeasureAtiAngleN;
          for (int r=0;r<this->m_NumberOfRadii;r++)                                      
          {  
            std::vector<float> currentBoundaryMeasureFrom0ToR;
            float currentMax=maxBoundaryMeasure[i];
            if ( fabsf(currentMax)>0.005 )
            {
              boundaryMeasureAtiAngleN.push_back( fabsf(allBoundaryMeasure[i][r]/currentMax) );
              if ( fabsf(allBoundaryMeasure[i][r]/currentMax) > 1.0f )
              {
                std::cout << fabsf(currentMax) << std::endl;
              }
            }
            else
            {
              boundaryMeasureAtiAngleN.push_back( 0.0f );
            }

          } // radii

          allBoundaryMeasureN.push_back(boundaryMeasureAtiAngleN);
        } // angles

        float maxMedialness=0;
        for (int ri=0;ri<this->m_NumberOfRadii;ri++)                                      
        {   
          float currentMedialness=0;
          for (int i=0; i<this->m_NumberOfAngles;i++)
          {  
            currentMedialness=currentMedialness+allBoundaryMeasureN[i][ri];
          }// angles 
          if (currentMedialness>maxMedialness)
          {
            maxMedialness=currentMedialness;
          }
        }	// radii
        
        outputIt.Set((OutputPixelType)maxMedialness/m_NumberOfAngles);
      } // use current voxel
      else
      {
        //outputIt.SetIndex(idx);
        outputIt.Set((OutputPixelType)0);
      } 
  }// faceList
}

} // namespace itk
#endif
