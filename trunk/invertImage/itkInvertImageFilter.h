#ifndef __itkInvertImageFilter_H_
#define __itkInvertImageFilter_H_
/** The above makes sure the file is never read twice */

#include "itkImageToImageFilter.h"

namespace itk
{

  /**
   * \class InvertImageFilter
   *
   * \ingroup CourseExamples 
   */

template < class TInputImage, class TOutputImage >
class InvertImageFilter : 
  public ImageToImageFilter< TInputImage, TOutputImage >
{
public:

  /** Standard ITK typedefs.*/
  typedef InvertImageFilter<  TInputImage, 
                              TOutputImage > Self;
  typedef ImageToImageFilter< TInputImage,  
                              TOutputImage > Superclass;
  typedef SmartPointer< Self >               Pointer;
  typedef SmartPointer< const Self >         ConstPointer;

  itkNewMacro( Self );
  itkTypeMacro( InvertImageFilter, ImageToImageFilter );

  /** Alias for template parameters */
  typedef TInputImage   InputImageType;
  typedef TOutputImage  OutputImageType;

  /** Useful additional typedefs */
  typedef typename InputImageType::PixelType    InputPixelType;
  typedef typename OutputImageType::PixelType   OutputPixelType;
      
  itkGetMacro( Constant, float );
  itkSetMacro( Constant, float );
  itkGetMacro( Exponential, int );
  itkSetMacro( Exponential, int );

protected:

  InvertImageFilter();
  virtual ~InvertImageFilter();

  //virtual void BeforeThreadedGenerateData();
  virtual void ThreadedGenerateData(const typename Superclass::OutputImageRegionType& outputRegionForThread, int threadId);

private:

  InvertImageFilter( const Self& );      // purposely not implemented
  void operator=( const Self& );  // purposely not implemented

  /** A private member variable that stores the constant to subtract. */
  int m_Exponential;
  InputPixelType m_Constant;

  }; // end class SubtractConstantImageFilter

} // end namespace itk

/** Include the txx file. This is necessary for templated classes */
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkinvertImageFilter.txx"
#endif

#endif //__itkInvertImageFilter_H_
