//
//
#ifndef __LumenIntensitySimilarityFilter__H
#define __LumenIntensitySimilarityFilter__H
#include "itkImageSource.h"
#include "itkConceptChecking.h"
#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkGaussianDistribution.h"


namespace itk {

template < typename TInputImage, typename TOutputImage > 
class LumenIntensitySimilarityFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:

  typedef LumenIntensitySimilarityFilter Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  itkTypeMacro(LumenIntensitySimilarityFilter, ImageToImageFilter);
  itkNewMacro(Self);

  typedef TInputImage  InputImageType;
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

  typedef float InternalPrecision;
  typedef Point<InternalPrecision,InputImageType::ImageDimension>      PrecisionPointType;
  typedef std::vector<PrecisionPointType>   SeedListType;
  typedef std::vector< InternalPrecision >  RadiusListType;

  itkGetMacro(DarkLumen, int);
  itkSetMacro(DarkLumen, int);
  itkSetMacro(Verbose, bool);
  itkGetMacro(Verbose, bool);

  virtual void SetSeedList( SeedListType seedList );
  virtual SeedListType GetSeedList();

  virtual void SetRadii( RadiusListType radii );
  virtual RadiusListType GetRadii();

private:
  SeedListType  m_SeedList;
  int            m_DarkLumen;
  RadiusListType m_RadiusList;
  bool           m_Verbose;

protected:
  void operator=( const Self& ); 

  LumenIntensitySimilarityFilter(const Self &);

  /** Constructor 
  * Note that it is declared protected. This makes sure that 
  * users can only instantiate this class with the ::New() method. */
  LumenIntensitySimilarityFilter();

  /** Destructor */
  virtual ~LumenIntensitySimilarityFilter() {}

  //virtual void BeforeThreadedGenerateData();
  //virtual void ThreadedGenerateData(const typename Superclass::OutputImageRegionType& outputRegionForThread, int threadId);
  virtual void GenerateData();
};

} // namespace itk

#include "itkLumenIntensitySimilarityFilter.txx"

#endif //__LumenIntensitySimilarityFilter__H
