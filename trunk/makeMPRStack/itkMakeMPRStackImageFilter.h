//
//
#ifndef __MakeMPRStackImageFilter__H
#define __MakeMPRStackImageFilter__H
#include "itkConceptChecking.h"
#include "itkImageToImageFilter.h"
#include "itkImage.h"


namespace itk {

template < typename TInputImage, typename TOutputImage > 
class MakeMPRStackImageFilter : public ImageToImageFilter< TInputImage, TOutputImage >
{
public:

  typedef MakeMPRStackImageFilter Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  itkTypeMacro(MakeMPRStackImageFilter, ImageToImageFilter);
  itkNewMacro(Self);

  typedef TInputImage   InputImageType;
  typedef TOutputImage  OutputImageType;

  typedef float InternalPrecision;
  typedef Point<InternalPrecision,InputImageType::ImageDimension> PrecisionPointType;
  typedef std::vector<PrecisionPointType>   PathType;

  itkSetMacro(Verbose, bool);
  itkGetMacro(Verbose, bool);
  itkSetMacro(SaveIntermediateResults, bool);
  itkGetMacro(SaveIntermediateResults, bool);


  virtual void SetPath( PathType path );
  virtual PathType GetPath();


  void UpdateParameters();

private:
  PathType   m_Path;
  bool       m_Verbose;
  bool       m_SaveIntermediateResults;


protected:
  void operator=( const Self& ); 

  MakeMPRStackImageFilter(const Self &);


  //* Constructor 
  MakeMPRStackImageFilter();

  /** Destructor */
  virtual ~MakeMPRStackImageFilter() {}

  virtual void GenerateData();
};

} // namespace itk

#include "itkMakeMPRStackImageFilter.txx"

#endif //__MakeMPRStackImageFilter__H
