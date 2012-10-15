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
#include "itkLinearInterpolateImageFunction.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#ifndef __itkGeodesicActiveContourAndRegionLevelSetFunction_h
#define __itkGeodesicActiveContourAndRegionLevelSetFunction_h
#define PI 3.1415926
#include "itkSegmentationLevelSetFunction.h"
namespace itk
{
	/** \class GeodesicActiveContourAndRegionLevelSetFunction
	*
	* \brief This function is used in GeodesicActiveContourLevelSetImageFilter to
	* segment structures in an image based on a user supplied edge potential map.
	*
	* \par IMPORTANT
	* The LevelSetFunction class contain additional information necessary
	* to gain full understanding of how to use this function.
	*
	* GeodesicActiveContourAndRegionLevelSetFunction is a subclass of the generic LevelSetFunction.
	* It is used to segment structures in an image based on a user supplied
	* edge potential map \f$ g(I) \f$, which
	* has values close to zero in regions near edges (or high image gradient) and values
	* close to one in regions with relatively constant intensity. Typically, the edge
	* potential map is a function of the gradient, for example:
	*
	* \f[ g(I) = 1 / ( 1 + | (\nabla * G)(I)| ) \f]
	* \f[ g(I) = \exp^{-|(\nabla * G)(I)|} \f]
	*
	* where \f$ I \f$ is image intensity and
	* \f$ (\nabla * G) \f$ is the derivative of Gaussian operator.
	*
	* The edge potential image is set via the SetFeatureImage() method.
	*
	* In this function both the propagation term \f$ P(\mathbf{x}) \f$
	* and the curvature spatial modifier term \f$ Z(\mathbf{x}) \f$ are taken directly
	* from the edge potential image such that:
	*
	* \f[ P(\mathbf{x}) = g(\mathbf{x}) \f]
	* \f[ Z(\mathbf{x}) = g(\mathbf{x}) \f]
	*
	* An advection term \f$ \mathbf{A}(\mathbf{x}) \f$ is constructed
	* from the negative gradient of the edge potential image.
	*
	* \f[ \mathbf{A}(\mathbf{x}) = -\nabla g(\mathbf{x}) \f]
	*
	* This term behaves like a doublet attracting the contour to the edges.
	*
	* This implementation is based on:
	* "Geodesic Active Contours",
	* V. Caselles, R. Kimmel and G. Sapiro.
	* International Journal on Computer Vision,
	* Vol 22, No. 1, pp 61-97, 1997
	*
	* \sa LevelSetFunction
	* \sa SegmentationLevelSetImageFunction
	* \sa GeodesicActiveContourLevelSetImageFilter
	*
	* \ingroup FiniteDifferenceFunctions
	* \ingroup ITKLevelSets
	*/
	template< class TImageType, class TFeatureImageType = TImageType,class TSimilarityImageType = TImageType >
	class ITK_EXPORT GeodesicActiveContourAndRegionLevelSetFunction:
		public SegmentationLevelSetFunction< TImageType, TFeatureImageType >
	{
	public:
		/** Standard class typedefs. */
		typedef GeodesicActiveContourAndRegionLevelSetFunction Self;
		typedef SegmentationLevelSetFunction< TImageType, TFeatureImageType >
			Superclass;
		typedef SmartPointer< Self >       Pointer;
		typedef SmartPointer< const Self > ConstPointer;
		typedef TFeatureImageType          FeatureImageType;

		/** Method for creation through the object factory. */
		itkNewMacro(Self);
		/** Run-time type information (and related methods) */
		itkTypeMacro(GeodesicActiveContourAndRegionLevelSetFunction, SegmentationLevelSetFunction);

		/** Extract some parameters from the superclass. */
		typedef typename Superclass::ImageType         ImageType;
		typedef typename TImageType::PixelType         InputPixelType;
		typedef typename TImageType::IndexType         IndexType;
		typedef typename Superclass::ContinuousIndexType  ContinuousIndexType;
		typedef typename Superclass::NeighborhoodType  NeighborhoodType;
		typedef typename Superclass::ScalarValueType   ScalarValueType;
		typedef typename Superclass::FeatureScalarType FeatureScalarType;
		typedef typename Superclass::RadiusType        RadiusType;
		typedef typename Superclass::FloatOffsetType   FloatOffsetType;
		typedef typename Superclass::VectorImageType   VectorImageType;
		typedef typename Superclass::GlobalDataStruct  GlobalDataStruct;
		typedef typename Superclass::TimeStepType TimeStepType;
		typedef typename Superclass::InterpolatorType InterpolatorType;
		struct ACGlobalDataStruct : public GlobalDataStruct {
			ScalarValueType m_MaxGlobalChange;
		};

		/** Extract some parameters from the superclass. */
		itkStaticConstMacro(ImageDimension, unsigned int,
			Superclass::ImageDimension);

		/** Compute speed image from feature image */
		virtual void CalculateSpeedImage();

		/** Compute the advection field from feature image. */
		virtual void CalculateAdvectionImage();


		/** The curvature speed is same as the propagation speed. */
		virtual ScalarValueType CurvatureSpeed(const NeighborhoodType & neighborhood,
			const FloatOffsetType & offset, GlobalDataStruct *gd) const
		{
			return this->GACPropagationSpeed(neighborhood, offset, gd);
		}
		virtual ScalarValueType PropagationSpeed(const NeighborhoodType & neighborhood, const FloatOffsetType & offset, GlobalDataStruct *gd) const;
		virtual ScalarValueType GACPropagationSpeed(const NeighborhoodType & neighborhood,
			const FloatOffsetType & offset, GlobalDataStruct *gd) const;

		/** Set/Get the sigma for the Gaussian kernel used to compute the gradient
		* of the feature image needed for the advection term of the equation. */


		void SetDerivativeSigma(const double v)
		{ m_DerivativeSigma = v; }
		double GetDerivativeSigma()
		{ return m_DerivativeSigma; }

		virtual void Initialize(const RadiusType & r)
		{
			Superclass::Initialize(r);

			this->SetAdvectionWeight(NumericTraits< ScalarValueType >::One);
			this->SetPropagationWeight(NumericTraits< ScalarValueType >::One);
			this->SetCurvatureWeight(NumericTraits< ScalarValueType >::One);
		}

		void SetHeavisideEpsilon(const double v)
		{ m_HeavisideEpsilon = v; }
		double GetHeavisideEpsilon()const
		{ return m_HeavisideEpsilon; }


		void SetInternalRegionWeight(const double v) 
		{ m_InternalRegionWeight = v; }
		double GetInternalRegionWeight()const
		{ return m_InternalRegionWeight; }

		void SetExternalRegionWeight(const double v)
		{ m_ExternalRegionWeight = v; }
		double GetExternalRegionWeight()const
		{ return m_ExternalRegionWeight; }

		virtual const ImageType * GetSimilarityImage() const
		{ return m_SimilarityImage.GetPointer(); }
		virtual void SetSimilarityImage(const ImageType *f)
		{     
			m_SimilarityImage=f;
			m_SimilarityImageInterpolator->SetInputImage(m_SimilarityImage);
			/* ImageRegionConstIterator< ImageType >
			fit( f, this->GetFeatureImage()->GetRequestedRegion() );
			ImageRegionIterator< ImageType >
			sit( m_SimilarityImage, this->GetFeatureImage()->GetRequestedRegion() );


			for ( fit.GoToBegin(), sit.GoToBegin(); !fit.IsAtEnd(); ++sit, ++fit )
			{
			sit.Set( static_cast< ScalarValueType >( fit.Get() ) );
			}  */

		}
		virtual const ImageType * GetIntensityImage() const
		{ return m_IntensityImage.GetPointer(); }
		virtual void SetIntensityImage(const ImageType *f)
		{  
			m_IntensityImage=f;
			m_IntensityImageInterpolator->SetInputImage(m_IntensityImage);
		}

		/*virtual const ImageType * GetPropergationImage() const
		{ return m_PropergationImage.GetPointer(); }*/
		/*virtual void SetPropergationImage(const ImageType *f)
		{  
		m_PropergationImage=f;
		m_PropergationInterpolator->SetInputImage(m_PropergationImage);
		}*/

	protected:
		GeodesicActiveContourAndRegionLevelSetFunction()
		{
			//m_SpeedImage = ImageType::New();
			// m_PropergationInterpolator = InterpolatorType::New();
			//m_PropergationInterpolator = InterpolatorType::New();
			m_IntensityImageInterpolator = InterpolatorType::New();
			m_SimilarityImageInterpolator = InterpolatorType::New();
			this->SetAdvectionWeight(NumericTraits< ScalarValueType >::One);
			this->SetPropagationWeight(NumericTraits< ScalarValueType >::One);
			this->SetCurvatureWeight(NumericTraits< ScalarValueType >::One);
			//m_Interpolator = InterpolatorType::New();
			m_DerivativeSigma = 1.0;
			m_SimilarityImage = ImageType::New();
			m_IntensityImage = ImageType::New();

		}
		/*virtual ImageType * GetSpeedImage()
		{ return m_SpeedImage.GetPointer(); }
		void SetSpeedImage(ImageType *s);*/

		virtual ~GeodesicActiveContourAndRegionLevelSetFunction() {}

		GeodesicActiveContourAndRegionLevelSetFunction(const Self &); //purposely not
		// implemented
		void operator=(const Self &);                        //purposely not
		// implemented

		void PrintSelf(std::ostream & os, Indent indent) const
		{
			Superclass::PrintSelf(os, indent);
			os << indent << "DerivativeSigma: " << m_DerivativeSigma << std::endl;
		}


		ScalarValueType calculateH(const ScalarValueType & z) const
		{
			return 0.5*(1.0 + (2.0/	PI) * atan(-z / this->GetHeavisideEpsilon()));
		}
		ScalarValueType calculatedH(const ScalarValueType & z) const
		{
			return (1.0/PI)*(1.0/(1.0 + (z*z) / (this->GetHeavisideEpsilon()*this->GetHeavisideEpsilon()) ));
		}

		/*ScalarValueType ComputeInternalRegionalSpeed(const ContinuousIndexType &cdx) const
		{ 

			//std::cout<<"ComputeInternalRegionalSpeed()"<<cdx[0]<<cdx[1]<<cdx[2]<<std::endl;
			ScalarValueType intensity = m_IntensityImageInterpolator->EvaluateAtContinuousIndex(cdx);
			ScalarValueType mean = m_MeanImageInterpolator->EvaluateAtContinuousIndex(cdx);
			ScalarValueType intensityDifference=(intensity-mean)*(intensity-mean);	
			ScalarValueType normalization=1-exp(0-intensityDifference/8100); //the regional term thus is normalized
			return ( static_cast< ScalarValueType > (normalization*normalization) );
		}*/
		/*ScalarValueType ComputeExternalRegionalSpeed(const ContinuousIndexType &cdx) const
		{

			// std::cout<<"ComputeExternalRegionalSpeed()"<<cdx[0]<<cdx[1]<<cdx[2]<<std::endl;
			ScalarValueType intensity = m_IntensityImageInterpolator->EvaluateAtContinuousIndex(cdx);
			ScalarValueType mean = m_MeanImageInterpolator->EvaluateAtContinuousIndex(cdx);
			ScalarValueType intensityDifference=(intensity-mean)*(intensity-mean);
			ScalarValueType normalization=exp(0-intensityDifference/8100); //the regional term thus is normalized
			return ( static_cast< ScalarValueType > (normalization*normalization));
		}*/


	private:
		 //typename InterpolatorType::Pointer m_Interpolator;
		 //typename ImageType::Pointer m_SpeedImage;
		//  typename InterpolatorType::Pointer m_PropergationInterpolator;
		typename InterpolatorType::Pointer m_IntensityImageInterpolator;
		typename InterpolatorType::Pointer m_SimilarityImageInterpolator;
		InputPixelType m_DerivativeSigma;
		InputPixelType m_HeavisideEpsilon;
		InputPixelType m_InternalRegionWeight;
		InputPixelType m_ExternalRegionWeight;
		typename ImageType::ConstPointer m_SimilarityImage;
		typename ImageType::ConstPointer m_IntensityImage;

	};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGeodesicActiveContourAndRegionLevelSetFunction.hxx"
#endif

#endif
