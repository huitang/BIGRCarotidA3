//
//

#ifndef __MEDIALNESSIMAGEFILTER__H
#define __MEDIALNESSIMAGEFILTER__H
#include "itkImageSource.h"
#include "itkConceptChecking.h"
#include "itkImageToImageFilterDetail.h"
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"


namespace itk {

template< class TInputImage, class TMaskImage, class TOutputImage > 
class MedialnessImageFilter: public ImageToImageFilter<TInputImage,TOutputImage>
{
public:

    typedef MedialnessImageFilter
        Self;
    typedef ImageToImageFilter<TInputImage,TOutputImage>
        Superclass;
    typedef SmartPointer<Self>        Pointer;
    typedef SmartPointer<const Self>  ConstPointer;

    itkTypeMacro(MedialnessImageFilter, ImageToImageFilter);
    itkNewMacro(Self);
    
    typedef TInputImage  InputImageType;
    typedef typename InputImageType::ValueType        ValueType;    
    typedef typename InputImageType::IndexType IndexType;
    //typedef typename InputImageType::TimeStepType TimeStepType; 
    typedef typename InputImageType::SizeType    InputSizeType;
    typedef typename InputImageType::SpacingType InputSpacingType;
    typedef typename InputImageType::PointType   InputPointType;
    typedef typename InputImageType::PixelType   InputPixelType;
    //typedef typename InputImageType::Dimension   Dimension;


    typedef TMaskImage                          MaskImageType;
    typedef typename MaskImageType::Pointer     MaskImagePointer;
    typedef typename MaskImageType::RegionType  MaskRegionType;
    typedef typename MaskImageType::SizeType    MaskSizeType;
    typedef typename MaskImageType::SpacingType MaskSpacingType;
    typedef typename MaskImageType::PointType   MaskPointType;
    typedef typename MaskImageType::PixelType   MaskPixelType;
   // typedef typename MaskImageType::Dimension   Dimension;

    typedef typename itk::CovariantVector< float,InputImageType::ImageDimension>  VectorPixelType;
    typedef typename itk::Image< VectorPixelType, InputImageType::ImageDimension> GradientImageType;

    typedef TOutputImage                             OutputImageType;
   // typedef typename OutputImageType::Pointer        OutputImagePointer;
    typedef typename OutputImageType::PixelType      OutputPixelType;
    typedef typename OutputImageType::RegionType     OutputRegionType;
    typedef typename OutputImageType::SizeType       OutputSizeType;
    typedef typename OutputImageType::SizeValueType  OutputSizeValueType;
    typedef typename OutputImageType::IndexType      OutputIndexType;
    typedef typename OutputImageType::IndexValueType OutputIndexValueType;
    //typedef typename OutputImageType::Dimension   Dimension;

    itkGetMacro( UseMask, bool );
    itkSetMacro( UseMask, bool );
    itkBooleanMacro( UseMask );
    itkGetMacro( NumberOfAngles, int );
    itkSetMacro( NumberOfAngles, int );
    itkGetMacro( Sigma, float );
    itkSetMacro( Sigma, float );
    itkGetMacro( Threshold, float );
    itkSetMacro( Threshold, float );
    itkGetMacro( MinimalRadius, float );
    itkSetMacro( MinimalRadius, float );
    itkGetMacro( MaximalRadius, float );
    itkSetMacro( MaximalRadius, float );
    itkGetMacro( NumberOfRadii, int );
    itkSetMacro( NumberOfRadii, int );
    itkGetMacro( DarkLumen, float );
    itkSetMacro( DarkLumen, float );
    itkGetMacro( Verbose, bool );
    itkSetMacro( Verbose, bool );

    void SetMaskImage( MaskImageType * maskImage)
    {

    this->ProcessObject::SetNthInput( 1, const_cast< MaskImageType * >( maskImage ) );
    }


    MaskImageType * GetMaskImage()
    {
        return ( static_cast< MaskImageType* >( this->ProcessObject::GetInput(1) ) );
    }

private:
  typename GradientImageType::Pointer m_GradientImage;

  bool  m_UseMask;
  int   m_NumberOfAngles;
  float m_DarkLumen;
  float m_Sigma;
  float m_MinimalRadius;
  float m_MaximalRadius;
  int   m_NumberOfRadii;
  float m_Threshold;
  bool  m_Verbose;

  typename itk::Vector< float, 3> m_XAxis;
  typename itk::Vector< float, 3> m_YAxis;

protected:
  void operator=( const Self& ); 
   MedialnessImageFilter(const Self &);

  /** Constructor 
   * Note that it is declared protected. This makes sure that 
   * users can only instantiate this class with the ::New() method. */
  MedialnessImageFilter();

  /** Destructor */
  virtual ~MedialnessImageFilter() {}

  /** GenerateData
   * This method does the work. It is automatically invoked
   * when the user calls Update(); */
  virtual void BeforeThreadedGenerateData();
  virtual void AfterThreadedGenerateData();
  virtual void ThreadedGenerateData(const typename Superclass::OutputImageRegionType& outputRegionForThread, int threadId);
  //virtual void GenerateData(const typename Superclass::OutputImageRegionType& outputRegionForThread, int threadId);
};

} // ikt namespace

#include "itkMedialnessImageFilter.txx"

#endif //__MEDIALNESSIMAGEFILTER__H
