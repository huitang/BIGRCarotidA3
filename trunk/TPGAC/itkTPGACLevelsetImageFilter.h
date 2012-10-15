
//--------------------------------------------------------------------------------

#ifndef __itkTPGACLevelSetImageFilter_h_
#define __itkTPGACLevelSetImageFilter_h_

#include "itkGeodesicActiveContourLevelSetImageFilter.h"

namespace itk {

template <class TInputImage,
          class TFeatureImage,
          class TOutputPixelType = float >
class ITK_EXPORT TPGACLevelSetImageFilter
  : public GeodesicActiveContourLevelSetImageFilter< TInputImage,
                                                     TFeatureImage, 
                                                     TOutputPixelType>
{
public:
  /** Standard class typedefs */
  typedef TPGACLevelSetImageFilter Self;
  typedef GeodesicActiveContourLevelSetImageFilter<
    TInputImage, TFeatureImage, TOutputPixelType> Superclass;
  typedef SmartPointer<Self>  Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Inherited typedef from the superclass. */
  typedef typename Superclass::ValueType        ValueType;
  typedef typename Superclass::OutputImageType  OutputImageType;
  typedef typename Superclass::FeatureImageType FeatureImageType;
  typedef typename Superclass::IndexType        IndexType;
  typedef typename Superclass::TimeStepType     TimeStepType;
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(TPGACLevelSetImageFilter,
               GeodesicActiveContourLevelSetImageFilter);

  /** Method for creation through the object factory */
  itkNewMacro(Self);
     
    
protected:
  ~TPGACLevelSetImageFilter() {}
  TPGACLevelSetImageFilter();

  virtual void PrintSelf(std::ostream &os, Indent indent) const; 

  TPGACLevelSetImageFilter(const Self &); // purposely not implemented
  void operator=(const Self&); //purposely not implemented


  /** Overridden from the parent class to indroduce a constraint on
   *  surface flow under certain conditions. */
  inline virtual ValueType CalculateUpdateValue(const IndexType &idx,
                                                const TimeStepType &dt,
                                                const ValueType &value,
                                                const ValueType &change);
};

} // end namespace itk



#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTPGACLevelSetImageFilter.txx"
#endif

#endif
