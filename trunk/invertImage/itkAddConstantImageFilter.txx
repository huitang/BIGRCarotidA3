#ifndef __itkInvertImageFilter_txx_
#define __itkInvertImageFilter_txx_
/** The above makes sure the file is never read twice */

#include "itkInvertImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

namespace itk
{

  /** Constructor */
  template < class TInputImage, class TOutputImage >
  SubtractConstantImageFilter<TInputImage, TOutputImage>
    ::SubtractConstantImageFilter()
  {
    // Initialize constant to 0
    this->m_Constant = 0;
  }

  /** Destructor */
  template < class TInputImage, class TOutputImage >
  SubtractConstantImageFilter<TInputImage, TOutputImage>
    ::~SubtractConstantImageFilter()
  {
    //nothing
  }

  /** GenerateData */
  template < class TInputImage, class TOutputImage >
  void SubtractConstantImageFilter<TInputImage, TOutputImage>
    ::GenerateData( void )
  {
    /** Allocate memory */
    typename OutputImageType::Pointer output = this->GetOutput();
    output->SetBufferedRegion( output->GetRequestedRegion() );
    output->Allocate();

    /** Typedef iterators */
    typedef ImageRegionConstIterator< InputImageType >  InputIteratorType;
    typedef ImageRegionIterator< OutputImageType >  OutputIteratorType;
        
    /** Create iterators */
    InputIteratorType inIt( this->GetInput(), this->GetOutput()->GetRequestedRegion() );
    OutputIteratorType outIt( this->GetOutput(), this->GetOutput()->GetRequestedRegion() );

    /** out = in - constant; */
    inIt.GoToBegin();
    outIt.GoToBegin();
    while ( !inIt.IsAtEnd() )
    {
      const OutputPixelType outputValue = 
        static_cast<OutputPixelType>( inIt.Get() + this->m_Constant )^this->m_Exponential;
      outIt.Set( outputValue );
      ++inIt;
      ++outIt;
    }   
  }

} // end namespace itk

#endif // #ifndef __itkSubtractConstantImageFilter_txx_
 