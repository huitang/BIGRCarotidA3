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
#ifndef __itkGeodesicActiveContourAndRegionLevelSetImageFilter_hxx
#define __itkGeodesicActiveContourAndRegionLevelSetImageFilter_hxx

#include "itkGeodesicActiveContourAndRegionLevelSetImageFilter.h"

namespace itk
{
template< class TInputImage, class TFeatureImage,class TSimilarityImage, class TOutputType >
GeodesicActiveContourAndRegionLevelSetImageFilter< TInputImage, TFeatureImage, TSimilarityImage, TOutputType >
::GeodesicActiveContourAndRegionLevelSetImageFilter()
{
  /* Instantiate a geodesic active contour function and set it as the
    segmentation function. */
  m_GeodesicActiveContourAndRegionFunction = GeodesicActiveContourAndRegionFunctionType::New();

  this->SetSegmentationFunction(m_GeodesicActiveContourAndRegionFunction);

  /* Turn off interpolation. */
  this->InterpolateSurfaceLocationOff();
}

template< class TInputImage, class TFeatureImage,class TSimilarityImage, class TOutputType >
void
GeodesicActiveContourAndRegionLevelSetImageFilter< TInputImage, TFeatureImage,  TSimilarityImage, TOutputType >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << "GeodesicActiveContourAndRegionFunction: " <<  m_GeodesicActiveContourAndRegionFunction.GetPointer();
}

template< class TInputImage, class TFeatureImage, class TSimilarityImage, class TOutputType >
void
GeodesicActiveContourAndRegionLevelSetImageFilter< TInputImage, TFeatureImage,TSimilarityImage, TOutputType >
::GenerateData()
{
	/*if (this->m_SimilarityImage==NULL)
	{
		std::cout<<"No Similarity image defined"<<std::endl;
		 return EXIT_FAILURE;
	}
	if (this->m_IntensityImage==NULL)
	{
		std::cout<<"No original image defined"<<std::endl;
		 return EXIT_FAILURE;
	}*/
  // Make sure the SpeedImage is setup for the case when PropagationScaling
  // is zero
  if ( this->GetSegmentationFunction()
       && this->GetSegmentationFunction()->GetPropagationWeight() == 0 )
    {
    this->GetSegmentationFunction()->AllocateSpeedImage();
    this->GetSegmentationFunction()->CalculateSpeedImage();
    }

  // Continue with Superclass implementation
  Superclass::GenerateData();
}
} // end namespace itk

#endif
