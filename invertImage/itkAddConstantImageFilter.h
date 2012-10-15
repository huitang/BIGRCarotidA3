#ifndef __itkInvertImageFilter_H_
#define __itkInvertImageFilter_H_
/** The above makes sure the file is never read twice */

#include "itkImageToImageFilter.h"

namespace itk
{

  /**
   * \class SubtractConstantImageFilter
   * \brief A class that subtracts a constant from an image.
   *
   * This class is a very basic example of an ITK filter.
   * It inherits from itk::ImageToImageFilter and
   * implements the GenerateData method.
   *
   * \ingroup CourseExamples 
   */

  template < class TInputImage, class TOutputImage >
  class SubtractConstantImageFilter : 
    public ImageToImageFilter< TInputImage, TOutputImage >
  {
  public:

    /** Standard ITK typedefs.*/
    typedef SubtractConstantImageFilter<
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
    itkTypeMacro( SubtractConstantImageFilter, ImageToImageFilter );

    /** Alias for template parameters */
    typedef TInputImage            InputImageType;
    typedef TOutputImage           OutputImageType;

    /** Useful additional typedefs */
    typedef typename InputImageType::PixelType    InputPixelType;
    typedef typename OutputImageType::PixelType   OutputPixelType;
        
    /** Set/Get the constant to subtract. */
    itkSetMacro( Constant, InputPixelType );
    itkGetConstMacro( Constant, InputPixelType );
	
	void GetExponential(int exponential)
	{
	m_exponential=exponential;
	}
	
	int GetExponential()
	{
	return m_exponential;
	}

	void GetConstant(InputPixelType constant)
	{
	m_constant=constant;
	}
	
	int GetConstant()
	{
	return m_constant;
	}
  protected:

    /** Constructor 
     * Note that it is declared protected. This makes sure that 
     * users can only instantiate this class with the ::New() method. */
    SubtractConstantImageFilter();

    /** Destructor */
    virtual ~SubtractConstantImageFilter();

    /** GenerateData
     * This method does the work. It is automatically invoked
     * when the user calls Update(); */
    virtual void GenerateData(void);

  private:

    /** The following constructor and assignment operator are 
     * not implemented, because otherwise the SmartPointer
     * mechanism wouldn't work anymore */
    SubtractConstantImageFilter( const Self& );      // purposely not implemented
    void operator=( const Self& );  // purposely not implemented

    /** A private member variable that stores the constant to subtract. */
	int m_exponential;
    InputPixelType m_Constant;

  }; // end class SubtractConstantImageFilter

} // end namespace itk

/** Include the txx file. This is necessary for templated classes */
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSubtractConstantImageFilter.txx"
#endif

#endif // end #ifndef __itkSubtractConstantImageFilter_H_