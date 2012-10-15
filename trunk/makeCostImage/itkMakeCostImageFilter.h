//
//
#ifndef __MakeCostImageFilter__H
#define __MakeCostImageFilter__H
#include "itkConceptChecking.h"
#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkMedialnessImageFilter.h"
#include "itkLumenIntensitySimilarityFilter.h"
#include <itkMaximumImageFilter.h>
#include "itkMultiplyImageFilter.h"
#include "itkInvertImageFilter.h"


namespace itk {

template < typename TInputImage0, typename TInputImage1, typename TMaskImage, typename TOutputImage > 
class MakeCostImageFilter : public ImageToImageFilter< TInputImage0, TOutputImage >
{
public:

  typedef MakeCostImageFilter Self;
  typedef ImageToImageFilter<TInputImage0,TOutputImage> Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  itkTypeMacro(MakeCostImageFilter, ImageToImageFilter);
  itkNewMacro(Self);

  typedef TInputImage0  BBImageType;
  typedef TInputImage1  PCImageType;
  typedef TMaskImage    MaskImageType;
  typedef TOutputImage  OutputImageType;

  typedef float InternalPrecision;
  typedef Image<InternalPrecision,BBImageType::ImageDimension> PrecisionImageType;
  typedef Point<InternalPrecision,BBImageType::ImageDimension> PrecisionPointType;
  typedef std::vector<PrecisionPointType>   SeedListType;
  typedef std::vector< InternalPrecision >  RadiusListType;


  itkSetMacro(UseMask, bool);
  itkGetMacro(UseMask, bool);
  itkGetMacro(DarkLumen, int);
  itkSetMacro(DarkLumen, int);
  itkSetMacro(Verbose, bool);
  itkGetMacro(Verbose, bool);
  itkSetMacro(UseOnlyBBImage, bool);
  itkGetMacro(UseOnlyBBImage, bool);
  itkSetMacro(SaveIntermediateResults, bool);
  itkGetMacro(SaveIntermediateResults, bool);

  virtual void SetBlackBloodImage( BBImageType * bbImage );
  virtual BBImageType* GetBlackBloodImage();

  virtual void SetPhaseContrastImage( PCImageType * pcImage );
  virtual PCImageType* GetPhaseContrastImage();

  virtual void SetMaskImage( MaskImageType * maskImage );
  virtual MaskImageType* GetMaskImage();

  virtual void SetSeedList( SeedListType seedList );
  virtual SeedListType GetSeedList();

  virtual void SetRadii( RadiusListType radii );
  virtual RadiusListType GetRadii();

  void UpdateParameters();

private:
  bool           m_UseMask;
  SeedListType   m_SeedList;
  int            m_DarkLumen;
  bool           m_UseOnlyBBImage;
  RadiusListType m_RadiusList;
  bool           m_Verbose;
  bool           m_SaveIntermediateResults;


protected:
  void operator=( const Self& ); 

  MakeCostImageFilter(const Self &);

  // Used filters
  typedef MedialnessImageFilter< BBImageType, MaskImageType, PrecisionImageType> BBMedialnessFilterType;
  typedef MedialnessImageFilter< PCImageType, MaskImageType, PrecisionImageType> PCMedialnessFilterType;

  typedef MaximumImageFilter<PrecisionImageType, PrecisionImageType, PrecisionImageType> MaximumFilterType;
  
  typedef LumenIntensitySimilarityFilter<BBImageType, PrecisionImageType> BBLISFilterType;
  typedef LumenIntensitySimilarityFilter<PCImageType, PrecisionImageType> PCLISFilterType;
  
  typedef MultiplyImageFilter<PrecisionImageType,PrecisionImageType,PrecisionImageType> MultiplyImageFilterType;
  
  typedef InvertImageFilter<PrecisionImageType, PrecisionImageType > InverterType;

  typename BBMedialnessFilterType::Pointer m_BBMedialnessFilter;
  typename PCMedialnessFilterType::Pointer m_PCMedialnessFilter;
  typename MaximumFilterType::Pointer m_MedMaxFilter;
  typename BBLISFilterType::Pointer m_BBLisFilter;
  typename PCLISFilterType::Pointer m_PCLisFilter;
  typename MaximumFilterType::Pointer m_LisMaxFilter;
  typename MultiplyImageFilterType::Pointer m_MultiplyFilter;
  typename InverterType::Pointer m_InvertFilter;

  /** Constructor 
  * Note that it is declared protected. This makes sure that 
  * users can only instantiate this class with the ::New() method. */
  MakeCostImageFilter();

  /** Destructor */
  virtual ~MakeCostImageFilter() {}

  //virtual void BeforeThreadedGenerateData();
  //virtual void ThreadedGenerateData(const typename Superclass::OutputImageRegionType& outputRegionForThread, int threadId);
  virtual void GenerateData();
};

} // namespace itk

#include "itkMakeCostImageFilter.txx"

#endif //__MakeCostImageFilter__H
