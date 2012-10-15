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
    typedef InvertImageFilter<
      TInputImage, TOutputImage>             Self;
    typedef ImageToImageFilter<
      TInputImage, TOutputImage>             Superclass;
    typedef SmartPointer<Self>               Pointer;
    typedef SmartPointer<const Self>         ConstPointer;

    /** Method for creation through the object factory.
     * This macro can be found in itkMacro.h.
     * It declares the ::New() function. */
    itkNewMacro( Self );

    /** Run-time type information. This macro is defined
     * in itkMacro.h. It declares a method that returns
     * the name of the class. It always has as arguments the
     * name of the class itself and the name of the superclass. */
    itkTypeMacro( InvertImageFilter, ImageToImageFilter );

    /** Alias for template parameters */
    typedef TInputImage            InputImageType;
    typedef TOutputImage           OutputImageType;

    /** Useful additional typedefs */
    typedef typename InputImageType::PixelType    InputPixelType;
    typedef typename OutputImageType::PixelType   OutputPixelType;
        
void SetConstant(float constant)
{
m_Constant=constant;
}

void GetConstant()
{
return m_Constant;
}

void SetExponential(int exponential)
{
m_Exponential=exponential;
}

void GetExponential()
{
return m_Exponential;
}

  protected:

    /** Constructor 
     * Note that it is declared protected. This makes sure that 
     * users can only instantiate this class with the ::New() method. */
    InvertImageFilter();

    /** Destructor */
    virtual ~InvertImageFilter();

    /** GenerateData
     * This method does the work. It is automatically invoked
     * when the user calls Update(); */
    virtual void GenerateData(void);

  private:

    /** The following constructor and assignment operator are 
     * not implemented, because otherwise the SmartPointer
     * mechanism wouldn't work anymore */
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

#endif // end #ifndef __itkSubtractConstantImageFilter_H_