//
//
#ifndef __GaussianDistributionImageFilter__H
#define __GaussianDistributionImageFilter__H
#include "itkImageSource.h"
#include "itkConceptChecking.h"
#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkGaussianDistribution.h"


namespace itk {

template < typename TInputImage, typename TOutputImage > 
class GaussianDistributionImageFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:

  typedef GaussianDistributionImageFilter Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  itkTypeMacro(GaussianDistributionImageFilter, ImageToImageFilter);
  itkNewMacro(Self);

  typedef TInputImage  InputImageType;
  typedef typename InputImageType::ValueType   ValueType;    
  typedef typename InputImageType::IndexType   IndexType;

  typedef typename InputImageType::SizeType    InputSizeType;
  typedef typename InputImageType::SpacingType InputSpacingType;
  typedef typename InputImageType::PointType   InputPointType;
  typedef typename InputImageType::PixelType   InputPixelType;

  typedef TOutputImage                             OutputImageType;
  typedef typename OutputImageType::PixelType      OutputPixelType;
  typedef typename OutputImageType::RegionType     OutputRegionType;
  typedef typename OutputImageType::SizeType       OutputSizeType;
  typedef typename OutputImageType::SizeValueType  OutputSizeValueType;
  typedef typename OutputImageType::IndexType      OutputIndexType;
  typedef typename OutputImageType::IndexValueType OutputIndexValueType;

  typedef itk::Statistics::GaussianDistribution GaussianType;

  virtual void   SetMean(double mean);
  virtual double GetMean();

  virtual void   SetVariance( double var);
  virtual double GetVariance();

  virtual void   SetSigma( double sig);
  virtual double GetSigma();

  itkGetMacro(Verbose, bool);
  itkSetMacro(Verbose, bool);

private:
  double m_Mean;
  double m_Variance;
  GaussianType::Pointer m_Gaussian;
  bool   m_Verbose;

protected:
  void operator=( const Self& ); 

  GaussianDistributionImageFilter(const Self &);

  /** Constructor 
  * Note that it is declared protected. This makes sure that 
  * users can only instantiate this class with the ::New() method. */
  GaussianDistributionImageFilter();

  /** Destructor */
  virtual ~GaussianDistributionImageFilter() {}

  virtual void BeforeThreadedGenerateData();
  virtual void AfterThreadedGenerateData();
  virtual void ThreadedGenerateData(const typename Superclass::OutputImageRegionType& outputRegionForThread, int threadId);
};

} // namespace itk

#include "itkGaussianDistributionImageFilter.txx"

#endif //__GaussianDistributionImageFilter__H
