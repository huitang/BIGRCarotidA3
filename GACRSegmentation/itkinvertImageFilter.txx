#ifndef __itkInvertImageFilter_txx_
#define __itkInvertImageFilter_txx_
/** The above makes sure the file is never read twice */

#include "itkInvertImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageAdaptor.h"
#include "itkImageRegionIterator.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkCastImageFilter.h"

namespace itk
{

  /** Constructor */
  template < class TInputImage, class TOutputImage >
  InvertImageFilter<TInputImage, TOutputImage>
    ::InvertImageFilter()
  {
    // Initialize constant to 0
    this->m_Constant = 0;
	 this->m_Exponential = 1;
  }

  /** Destructor */
  template < class TInputImage, class TOutputImage >
  InvertImageFilter<TInputImage, TOutputImage>
    ::~InvertImageFilter()
  {
    //nothing
  }

  /** GenerateData */
  template < class TInputImage, class TOutputImage >
  void InvertImageFilter<TInputImage, TOutputImage>
    ::GenerateData( void )
  {
    /** Allocate memory */
    typename OutputImageType::Pointer output = this->GetOutput();
    output->SetBufferedRegion( output->GetLargestPossibleRegion() );
    output->Allocate();                  
    typedef ImageRegionConstIterator< InputImageType >  InputIteratorType;
    typedef ImageRegionIterator< OutputImageType >  OutputIteratorType;
        
    /** Create iterators */
    InputIteratorType inIt( this->GetInput(), this->GetOutput()->GetLargestPossibleRegion() );
    OutputIteratorType outIt( this->GetOutput(), this->GetOutput()->GetLargestPossibleRegion() );

    /** out = in - constant; */
    inIt.GoToBegin();
    outIt.GoToBegin();
    while ( !inIt.IsAtEnd() )
    {
      const OutputPixelType outputValue = static_cast<OutputPixelType>(1.0/(pow(inIt.Get(),this->m_Exponential) + this->m_Constant));
      outIt.Set( outputValue );
      ++inIt;
      ++outIt;
      }
    }   

} // end namespace itk

#endif // #ifndef __itkInvertImageFilter_txx_
 