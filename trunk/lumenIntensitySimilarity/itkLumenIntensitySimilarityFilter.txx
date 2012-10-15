//
//

#ifndef __LumenIntensitySimilarityFilter__TXX
#define __LumenIntensitySimilarityFilter__TXX

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkConstShapedNeighborhoodIterator.h"
#include "itkThresholdImageFilter.h"
#include "itkGaussianDistributionImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
namespace itk {


template< class TInputImage, class TOutputImage > 
LumenIntensitySimilarityFilter< TInputImage, TOutputImage >
::LumenIntensitySimilarityFilter(){
}


template< class TInputImage, class TOutputImage > 
void
LumenIntensitySimilarityFilter< TInputImage, TOutputImage >
::SetSeedList( typename LumenIntensitySimilarityFilter< TInputImage, TOutputImage >::SeedListType seedList ){
  m_SeedList = seedList;
}

template< class TInputImage, class TOutputImage > 
typename LumenIntensitySimilarityFilter< TInputImage, TOutputImage >::SeedListType
LumenIntensitySimilarityFilter< TInputImage, TOutputImage >
::GetSeedList() {
  return m_SeedList;
}

template< class TInputImage, class TOutputImage > 
void
LumenIntensitySimilarityFilter< TInputImage, TOutputImage >
::SetRadii( typename LumenIntensitySimilarityFilter< TInputImage, TOutputImage >::RadiusListType radii ){
  m_RadiusList = radii;
}

template< class TInputImage, class TOutputImage > 
typename LumenIntensitySimilarityFilter< TInputImage, TOutputImage >::RadiusListType
LumenIntensitySimilarityFilter< TInputImage, TOutputImage >
::GetRadii() {
  return m_RadiusList;
}


template< class TInputImage, class TOutputImage > 
void LumenIntensitySimilarityFilter< TInputImage, TOutputImage >
::GenerateData()
{

  if ( m_Verbose ) {
    std::cout << "Start calculating lumen intenisty similarity" << std::endl;
  }

  typename InputImageType::Pointer  input  = const_cast<InputImageType*>( this->GetInput() );
  typename OutputImageType::Pointer output = this->GetOutput();

  if ( m_SeedList.size() < 2 ) {
    std::cout << "At least two seed points are needed, but only " << m_SeedList.size() << " are found." << std::endl;
    return;
  }

  typedef std::vector< typename InputImageType::IndexType > IndexListType;
  IndexListType indexList( m_SeedList.size() );
  typename IndexListType::iterator idxIt = indexList.begin();
  typename SeedListType::iterator seedIt = m_SeedList.begin();
  for (; seedIt < m_SeedList.end(); ++idxIt, ++seedIt){
    input->TransformPhysicalPointToIndex( *seedIt,*idxIt );
  }


  // Get image values at seed points
  InternalPrecision sum = static_cast<InternalPrecision>(0);
  InternalPrecision sSum = static_cast<InternalPrecision>(0);
  size_t nElements=0;
  InternalPrecision mean = 0.0f;
  InternalPrecision var = 0.0f;
  InternalPrecision stdev = 0.0f;

  typedef itk::ConstShapedNeighborhoodIterator< InputImageType > ShapedNeighborhoodIteratorType;
  RadiusListType::iterator radiusIt = m_RadiusList.begin();
  for ( idxIt=indexList.begin(); idxIt < indexList.end(); ++ idxIt, ++radiusIt ){
    unsigned int voxelRadius = static_cast<int>( (*radiusIt) / input->GetSpacing()[0]);
    typename ShapedNeighborhoodIteratorType::RadiusType elementRadius;
    elementRadius.Fill(voxelRadius);

    if (m_Verbose ){
      std::cout << "radius in physical units: " << *radiusIt << std::endl;
      std::cout << "element radius: " << elementRadius << std::endl;
    }

    ShapedNeighborhoodIteratorType it( elementRadius, input, input->GetLargestPossibleRegion() );
    float rad = static_cast<float>(voxelRadius);
    for (float z = -rad; z <= rad; z++){
      const float wz = input->GetSpacing()[2]*z;
      for (float y = -rad; y <= rad; y++){
        const float wy = input->GetSpacing()[1]*y;
        for (float x = -rad; x <= rad; x++){
          typename ShapedNeighborhoodIteratorType::OffsetType off;

          const float wx = input->GetSpacing()[0]*x;
          float dis = vcl_sqrt( wx*wx + wy*wy + wz*wz);
          if (dis <= *radiusIt ){
            off[0] = static_cast<int>(x);
            off[1] = static_cast<int>(y);
            off[2] = static_cast<int>(z);
            it.ActivateOffset(off);
          }
        }
      }
    }

    it.SetLocation( *idxIt );
    typename ShapedNeighborhoodIteratorType::ConstIterator ci;
    for (ci = it.Begin(); !ci.IsAtEnd(); ci++){
      sum += ci.Get();
      sSum += ci.Get()*ci.Get();
      nElements++;
    }
  }

  mean = sum/nElements;
  var = 1./(nElements-1.) * (sSum - nElements * mean* mean);
  stdev = sqrt(var);
  if (m_Verbose ){
    std::cout << "Seed point statistics (mean, stdev): " << mean << " , " << stdev << std::endl;
  }

  typedef double InternalPrecisionType;
  typedef typename itk::Image< InternalPrecisionType, 3> PrecisionImageType;

  typedef typename itk::ThresholdImageFilter< InputImageType > ThresholdFilterType;
  typename ThresholdFilterType::Pointer thFilter = ThresholdFilterType::New();
  thFilter->SetInput( input );
  thFilter->SetOutsideValue( mean );
  if ( m_DarkLumen == 0 ) {
    thFilter->ThresholdAbove( mean );
  } else {
    thFilter->ThresholdBelow( mean );
  }
  thFilter->Update();

  typedef itk::GaussianDistributionImageFilter< InputImageType,PrecisionImageType> GaussianDistributionFilterType;
  typename GaussianDistributionFilterType::Pointer gdFilter = GaussianDistributionFilterType::New();

  gdFilter->SetMean( mean );
  gdFilter->SetVariance( var );
  gdFilter->SetInput( thFilter->GetOutput() );
  gdFilter->Update();

  typedef itk::RescaleIntensityImageFilter< PrecisionImageType, OutputImageType > RescaleFilterType;
  typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput(gdFilter->GetOutput());
  rescaleFilter->SetOutputMinimum(0.0);
  rescaleFilter->SetOutputMaximum(1.0);

  rescaleFilter->GraftOutput( this->GetOutput() );
  rescaleFilter->Update();
  this->GraftOutput(rescaleFilter->GetOutput());

  if ( m_Verbose ) {
    std::cout << "End calculating lumen intenisty similarity" << std::endl;
  }

}

} // namespace itk
#endif
