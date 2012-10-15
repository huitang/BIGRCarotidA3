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
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

#define myMin(a,b) ((a) < (b) ? (a) : (b))
#define myMax(a,b) ((a) > (b) ? (a) : (b))

namespace itk {

template< class TInputImage, class TMaskImage, class TOutputImage > 
MedialnessImageFilter< TInputImage, TMaskImage, TOutputImage >
::MedialnessImageFilter(){
  m_UseMask = true;
  m_Verbose = true;
}

template< class TInputImage, class TMaskImage, class TOutputImage > 
void MedialnessImageFilter< TInputImage, TMaskImage, TOutputImage >
::BeforeThreadedGenerateData(){
  if ( m_Verbose ) {
    std::cout << "Start calculating medialness." << std::endl;
    std::cout << "Start calculating gradients." << std::endl;
  }
  typename InputImageType::IndexType originVoxel;
  typename InputImageType::IndexType xAxisVoxel;
  typename InputImageType::IndexType yAxisVoxel;
  originVoxel.Fill(0);
  xAxisVoxel[0]=1;
  xAxisVoxel[1]=0;
  xAxisVoxel[2]=0;
  yAxisVoxel[0]=0;
  yAxisVoxel[1]=1;
  yAxisVoxel[2]=0;
  itk::ImageBase<3>::PointType originWorld;
  itk::ImageBase<3>::PointType xAxisWorld;
  itk::ImageBase<3>::PointType yAxisWorld;
  this->GetInput()->TransformIndexToPhysicalPoint(originVoxel,originWorld);
  this->GetInput()->TransformIndexToPhysicalPoint(xAxisVoxel,xAxisWorld);
  this->GetInput()->TransformIndexToPhysicalPoint(yAxisVoxel,yAxisWorld);
  m_XAxis = xAxisWorld - originWorld;
  m_YAxis = yAxisWorld - originWorld;
  m_XAxis.Normalize();
  m_YAxis.Normalize();
  if ( m_Verbose ) {
    std::cout << "origin: " << originWorld << std::endl;
    std::cout << "xAxis: " << m_XAxis << std::endl;
    std::cout << "yAxis: " << m_YAxis << std::endl;
  }

  typedef typename itk::GradientRecursiveGaussianImageFilter <InputImageType,GradientImageType> GradientFilterType;
  typename GradientFilterType::Pointer gradientFilter =  GradientFilterType::New();
  gradientFilter->SetSigma( this->m_Sigma );
  gradientFilter->SetInput( this->GetInput() );
  gradientFilter->Update();
  m_GradientImage = gradientFilter->GetOutput();

  if ( m_Verbose ) {
    std::cout << "End calculating gradients." << std::endl;
  }
}

template< class TInputImage, class TMaskImage, class TOutputImage > 
void MedialnessImageFilter< TInputImage, TMaskImage, TOutputImage >
::AfterThreadedGenerateData(){
  if ( m_Verbose ) {
    std::cout << "End calculating medialness." << std::endl;
  }
}

template< class TInputImage, class TMaskImage, class TOutputImage > 
void MedialnessImageFilter< TInputImage, TMaskImage, TOutputImage >
::ThreadedGenerateData(const typename ImageToImageFilter<TInputImage,TOutputImage>::OutputImageRegionType& outputRegionForThread, int threadId)
{
  typename OutputImageType::Pointer output = this->GetOutput();


    //get image spacing and size
  typename InputImageType::IndexType idx;
  typename InputImageType::IndexType neibourIdx;
  std::vector<float> regionExt(3, 0);
  std::vector<float> imgExt(3, 0);
  regionExt[0] = outputRegionForThread.GetSize()[0];
  regionExt[1] = outputRegionForThread.GetSize()[1];
  regionExt[2] = outputRegionForThread.GetSize()[2];
  imgExt[0] = this->GetInput()->GetBufferedRegion().GetSize()[0];
  imgExt[1] = this->GetInput()->GetBufferedRegion().GetSize()[1];
  imgExt[2] = this->GetInput()->GetBufferedRegion().GetSize()[2];
  const typename InputImageType::SpacingType & sp = this->GetOutput()->GetSpacing();
  std::vector<float> voxelSizes (3, 1.0f);
  voxelSizes[0] = sp[0];
  voxelSizes[1] = sp[1];
  voxelSizes[2] = sp[2];

  //Region Iterator
  typedef typename itk::ImageRegionIteratorWithIndex <GradientImageType> GradientIteratorType;
  GradientIteratorType gradIt( m_GradientImage, output->GetRequestedRegion() );
  gradIt.GoToBegin(); 

  typedef typename itk::ImageRegionIteratorWithIndex <InputImageType> ImageIteratorType;
  ImageIteratorType maskIt;
  if ( m_UseMask && this->GetMaskImage() ) {
    maskIt = ImageIteratorType(this->GetMaskImage(),output->GetRequestedRegion());
    maskIt.GoToBegin();
  }

  ImageIteratorType outputIt(this->GetOutput(),outputRegionForThread);
  outputIt.GoToBegin();

  //calculatemedialness
  float deltaRadius= (this->m_MaximalRadius-this->m_MinimalRadius)/(this->m_NumberOfRadii-1);
  float stepOfAngle =2.0f*3.14159265358979323846 / static_cast<float>(this->m_NumberOfAngles);
  std::vector<std::vector<float> > shiftx(m_NumberOfAngles,std::vector<float>(m_NumberOfRadii,0.0f));
  std::vector<std::vector<float> > shifty(m_NumberOfAngles,std::vector<float>(m_NumberOfRadii,0.0f));

  // Detrmine indices
  for (int i=0; i<this->m_NumberOfAngles;i++)
  {
    const float currentAngle = i*stepOfAngle;
    for (int ri=0;ri<this->m_NumberOfRadii;ri++)                                      
    { 
      const float currentRadius = this->m_MinimalRadius+ri*deltaRadius;;
      shiftx[i][ri] = currentRadius*cos(currentAngle)/voxelSizes[0];
      shifty[i][ri] = currentRadius*sin(currentAngle)/voxelSizes[1];
    }
  }

  // Loop over region
  typename ImageToImageFilter<TInputImage,TOutputImage>::OutputImageRegionType::IndexType upperIndex;
  upperIndex = outputRegionForThread.GetIndex()+outputRegionForThread.GetSize();
  for (int zi = outputRegionForThread.GetIndex()[2]; zi<upperIndex[2];  ++zi){
    for (float yi = outputRegionForThread.GetIndex()[1];  yi<upperIndex[1];  ++yi){
      for (float xi = outputRegionForThread.GetIndex()[0]; xi<upperIndex[0];  ++xi){  
        idx[0]=static_cast<int>(xi);
        idx[1]=static_cast<int>(yi);
        idx[2]=zi;
        gradIt.SetIndex(idx);

        bool useCurrentVoxel = true;
        if ( m_UseMask ) {
          maskIt.SetIndex(idx);
          if ( maskIt.Get() <= m_Threshold ){
            useCurrentVoxel = false;
          }
        }

        if ( useCurrentVoxel )
        {
          std::vector<std::vector<float> > allBoundaryMeasure( m_NumberOfAngles,std::vector<float>(m_NumberOfRadii,0) );
          std::vector<std::vector<float> > allIntensityProfile;
          std::vector<float> maxBoundaryMeasure(this->m_NumberOfAngles,0.0f);

          for (int i=0; i<this->m_NumberOfAngles;++i)
          {
            const float currentAngle = static_cast<float>(i)*stepOfAngle;
            for (int ri=0;ri<this->m_NumberOfRadii;++ri)                                      
            {    
              // Get gradient
              float rX = shiftx[i][ri];
              float rY = shifty[i][ri];
              neibourIdx[0]=static_cast<int>( myMin( myMax(xi+rX,0.0f), imgExt[0]-1.0f) );
              neibourIdx[1]=static_cast<int>( myMin( myMax(yi+rY,0.0f), imgExt[1]-1.0f) );
              neibourIdx[2]=zi;
              gradIt.SetIndex(neibourIdx);
              VectorPixelType grad = gradIt.Get();
              const float graX = (grad*m_XAxis);
              const float graY = (grad*m_YAxis);
              //const float rS = sqrt(rX*rX+rY*rY);
              //rX /= rS;
              //rY /= rS;
              const float gra = graX*cos(currentAngle)+graY*sin(currentAngle);// here it is actually <G,n>
              //const float gra = graX*rX+graY*rY;// here it is actually <G,n>
              
              // negative gradient will not be taken into account, when the vessel voxel is black. so all gradient vectors are outward, so <G,n>is negtive
              if (this->m_DarkLumen==1.0f)
              {
                allBoundaryMeasure[i][ri] = myMax(gra,0.0f);
              }

              if (this->m_DarkLumen==0.0f)
              { 
                allBoundaryMeasure[i][ri] = myMin(gra,0.0f);
              }
              
              //for those cases whose lumen intensity is not black or bright but in between
              if (this->m_DarkLumen > 1.0f)
              { 
                allBoundaryMeasure[i][ri] = gra;
              }

              if (fabsf(gra) > maxBoundaryMeasure[i] )
              {
                maxBoundaryMeasure[i]=fabsf(gra);
              }
            }
          }

          std::vector<std::vector<float> > allBoundaryMeasureN( m_NumberOfAngles,std::vector<float>(m_NumberOfRadii,0) );     
          for (int i=0; i<this->m_NumberOfAngles;i++)
          {
            for (int r=0;r<this->m_NumberOfRadii;r++)                                      
            {  
              allBoundaryMeasureN[i][r] = (fabs(maxBoundaryMeasure[i])>0.005 ? fabsf( allBoundaryMeasure[i][r] / maxBoundaryMeasure[i]) : 0.0f);
            }
          }

          float maxMedialness=0.0f;
          for (int ri=0;ri<this->m_NumberOfRadii;ri++)                                      
          {   
            float currentMedialness=0.0f;
            for (int i=0; i<this->m_NumberOfAngles;i++)
            {  
              currentMedialness+=allBoundaryMeasureN[i][ri];
            }
            maxMedialness = myMax(maxMedialness, currentMedialness);
          }	
          outputIt.SetIndex(idx);
          outputIt.Set((OutputPixelType)maxMedialness/m_NumberOfAngles);
        }
        else
        {
          outputIt.SetIndex(idx);
          outputIt.Set(static_cast<OutputPixelType>(0));
        }
      }
    }
  }
}
} // namespace itk
#endif
