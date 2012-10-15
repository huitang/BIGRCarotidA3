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
#ifndef __itkGeodesicActiveContourAndRegionLevelSetFunction_hxx
#define __itkGeodesicActiveContourAndRegionLevelSetFunction_hxx

#include "itkGeodesicActiveContourAndRegionLevelSetFunction.h"
#include "itkImageRegionIterator.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkGradientImageFilter.h"
#include "itkVectorCastImageFilter.h"
#define	PI 3.1415926

namespace itk
{

template< class TImageType, class TFeatureImageType, class TSimilarityImageType >
void GeodesicActiveContourAndRegionLevelSetFunction< TImageType, TFeatureImageType,TSimilarityImageType >
::CalculateSpeedImage()
{
  /* copy the feature image into the speed image */
  ImageRegionConstIterator< FeatureImageType >
  fit( this->GetFeatureImage(), this->GetFeatureImage()->GetRequestedRegion() );
  ImageRegionIterator< ImageType >
  sit( this->GetSpeedImage(), this->GetFeatureImage()->GetRequestedRegion() );
  for ( fit.GoToBegin(), sit.GoToBegin(); !fit.IsAtEnd(); ++sit, ++fit )
  {
	  sit.Set( static_cast< ScalarValueType >( fit.Get() ) );
  }
}
template< class TImageType, class TFeatureImageType, class TSimilarityImageType >
typename GeodesicActiveContourAndRegionLevelSetFunction< TImageType, TFeatureImageType,TSimilarityImageType>::ScalarValueType  
GeodesicActiveContourAndRegionLevelSetFunction< TImageType, TFeatureImageType,TSimilarityImageType>
::PropagationSpeed(const NeighborhoodType & neighborhood,
				   const FloatOffsetType & offset, GlobalDataStruct *) const
{
	IndexType idx = neighborhood.GetIndex();
	//std::cout<<"PropergationSpeed()"<<idx[0]<<idx[1]<<idx[2]<<std::endl;
	ContinuousIndexType cdx;
	for ( unsigned i = 0; i < ImageDimension; ++i )
	{
		cdx[i] = static_cast< double >( idx[i] ) - offset[i];
	}
	//std::cout<<"PropergationSpeed() continous"<<cdx[0]<<cdx[1]<<cdx[2]<<std::endl;
    //ScalarValueType fIn = this->ComputeInternalRegionalSpeed(cdx);
	//ScalarValueType fEx = this->ComputeExternalRegionalSpeed(cdx);
	//std::cout<<fIn<< "  "<<fEx<<std::endl;



	if (this->m_SimilarityImageInterpolator->IsInsideBuffer(cdx) )
	{
		//return ( static_cast< ScalarValueType >( m_GACPropagationWeight*this->m_Interpolator->EvaluateAtContinuousIndex(cdx) + m_RegionWeight *f))  ;
		//return ( static_cast< ScalarValueType >( m_ExternalRegionWeight *fEx - m_InternalRegionWeight*fIn))  ;
		return ( static_cast< ScalarValueType >( m_ExternalRegionWeight *this->m_SimilarityImageInterpolator->EvaluateAtContinuousIndex(cdx) - m_InternalRegionWeight * (1-this->m_SimilarityImageInterpolator->EvaluateAtContinuousIndex(cdx) ))) ;

	}
	else { 
		//return ( static_cast< ScalarValueType >( m_GACPropagationWeight*this->m_SpeedImage->GetPixel(idx)  + m_RegionWeight * f  ) ); 
		//return ( static_cast< ScalarValueType >( m_ExternalRegionWeight *fEx  - m_InternalRegionWeight*fIn  ) ); 
		return ( static_cast< ScalarValueType >( m_ExternalRegionWeight *this->m_SimilarityImage->GetPixel(idx) - m_InternalRegionWeight * (1-this->m_SimilarityImage->GetPixel(idx))))  ;

	}
}
template< class TImageType, class TFeatureImageType, class TSimilarityImageType >
typename GeodesicActiveContourAndRegionLevelSetFunction< TImageType, TFeatureImageType,TSimilarityImageType>::ScalarValueType  
GeodesicActiveContourAndRegionLevelSetFunction< TImageType, TFeatureImageType,TSimilarityImageType>
::GACPropagationSpeed(const NeighborhoodType & neighborhood,
				   const FloatOffsetType & offset, GlobalDataStruct *) const
{
	IndexType idx = neighborhood.GetIndex();

	ContinuousIndexType cdx;

	for ( unsigned i = 0; i < ImageDimension; ++i )
	{
		cdx[i] = static_cast< double >( idx[i] ) - offset[i];
	}

	if ( this->m_Interpolator->IsInsideBuffer(cdx) )
	{
		return ( static_cast< ScalarValueType >(
			this->m_Interpolator->EvaluateAtContinuousIndex(cdx)))  ;
			//this->m_PropergationInterpolator->EvaluateAtContinuousIndex(cdx)));
	}
	else { 
		return ( static_cast< ScalarValueType >( this->m_SpeedImage->GetPixel(idx) ) ); 
		//return ( static_cast< ScalarValueType >( this->m_PropergationImage->GetPixel(idx) ) ); 
		std::cout<<"GACPropagationSpeed used"<<std::endl;
	}
}
template< class TImageType, class TFeatureImageType,class TSimilarityImageType >
void GeodesicActiveContourAndRegionLevelSetFunction< TImageType, TFeatureImageType,TSimilarityImageType  >
::CalculateAdvectionImage()
{
  /* compute the gradient of the feature image. */

  typename VectorImageType::Pointer gradientImage;

  if ( m_DerivativeSigma != NumericTraits< float >::Zero )
    {
    typedef GradientRecursiveGaussianImageFilter< FeatureImageType, VectorImageType >
    DerivativeFilterType;

    typename DerivativeFilterType::Pointer derivative = DerivativeFilterType::New();
    derivative->SetInput( this->GetFeatureImage() );
    derivative->SetSigma(m_DerivativeSigma);
    derivative->Update();

    gradientImage = derivative->GetOutput();
    }
  else
    {
    typedef GradientImageFilter< FeatureImageType > DerivativeFilterType;

    typename DerivativeFilterType::Pointer derivative = DerivativeFilterType::New();
    derivative->SetInput( this->GetFeatureImage() );
    derivative->SetUseImageSpacingOn();
    derivative->Update();

    typedef typename DerivativeFilterType::OutputImageType                      DerivativeOutputImageType;
    typedef VectorCastImageFilter< DerivativeOutputImageType, VectorImageType > GradientCasterType;

    typename GradientCasterType::Pointer caster = GradientCasterType::New();
    caster->SetInput( derivative->GetOutput() );
    caster->Update();

    gradientImage = caster->GetOutput();
    }

  /* copy negative gradient into the advection image. */
  ImageRegionIterator< VectorImageType >
  dit( gradientImage, this->GetFeatureImage()->GetRequestedRegion() );
  ImageRegionIterator< VectorImageType >
  ait( this->GetAdvectionImage(), this->GetFeatureImage()->GetRequestedRegion() );

  for ( dit.GoToBegin(), ait.GoToBegin(); !dit.IsAtEnd(); ++dit, ++ait )
    {
    typename VectorImageType::PixelType v = dit.Get();
    for ( unsigned int j = 0; j < ImageDimension; j++ )
      {
      v[j] *= -1.0L;
      }
    ait.Set(v);
    }
}
} // end namespace itk

#endif
