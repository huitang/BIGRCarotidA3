//
//

#ifndef __MakeCostImageFilter__TXX
#define __MakeCostImageFilter__TXX

#include "itkImageFileWriter.h" 

namespace itk {


template < typename TInputImage0, typename TInputImage1, typename TMaskImage, typename TOutputImage > 
MakeCostImageFilter< TInputImage0, TInputImage1, TMaskImage, TOutputImage > 
::MakeCostImageFilter(){
  m_SaveIntermediateResults = false;
  m_BBMedialnessFilter = BBMedialnessFilterType::New();
  m_PCMedialnessFilter = PCMedialnessFilterType::New();
  m_MedMaxFilter       = MaximumFilterType::New();
  m_BBLisFilter        = BBLISFilterType::New();
  m_PCLisFilter        = PCLISFilterType::New();
  m_LisMaxFilter       = MaximumFilterType::New();
  m_MultiplyFilter     = MultiplyImageFilterType::New();
  m_InvertFilter       = InverterType::New();

  m_MedMaxFilter->SetInput(0, m_BBMedialnessFilter->GetOutput() );
  m_MedMaxFilter->SetInput(1, m_PCMedialnessFilter->GetOutput() );
  m_LisMaxFilter->SetInput(0, m_BBLisFilter->GetOutput() );
  m_LisMaxFilter->SetInput(1, m_PCLisFilter->GetOutput() );
  m_MultiplyFilter->SetInput(0, m_MedMaxFilter->GetOutput());
  m_MultiplyFilter->SetInput(1, m_LisMaxFilter->GetOutput());
  m_InvertFilter->SetInput( m_MultiplyFilter->GetOutput() );

  m_BBMedialnessFilter->SetNumberOfAngles( 24 );
  m_BBMedialnessFilter->SetDarkLumen( 1 );
  m_BBMedialnessFilter->SetSigma( 1 );
  m_BBMedialnessFilter->SetMinimalRadius( 0 );
  m_BBMedialnessFilter->SetMaximalRadius( 4 );
  m_BBMedialnessFilter->SetNumberOfRadii( 20 );
  m_BBMedialnessFilter->SetThreshold( 0 );

  m_PCMedialnessFilter->SetNumberOfAngles( 24 );
  m_PCMedialnessFilter->SetDarkLumen( 0 );
  m_PCMedialnessFilter->SetSigma( 1 );
  m_PCMedialnessFilter->SetMinimalRadius( 0 );
  m_PCMedialnessFilter->SetMaximalRadius( 4 );
  m_PCMedialnessFilter->SetNumberOfRadii( 20 );
  m_PCMedialnessFilter->SetThreshold( 0 );

  m_UseOnlyBBImage = false;
  m_RadiusList.push_back( 3.5 );
  m_RadiusList.push_back( 2.5 );
  m_RadiusList.push_back( 2.5 );
  m_BBLisFilter->SetSeedList( m_SeedList );
  m_BBLisFilter->SetDarkLumen( 1 );
  m_BBLisFilter->SetRadii( m_RadiusList );
  m_BBLisFilter->SetVerbose( m_Verbose );

  m_PCLisFilter->SetSeedList( m_SeedList );
  m_PCLisFilter->SetDarkLumen( 0 );
  m_PCLisFilter->SetRadii( m_RadiusList );
  m_PCLisFilter->SetVerbose( m_Verbose );

  m_InvertFilter->SetConstant( 0.00000001 );
  m_InvertFilter->SetExponential( 4 );
}

// Input 0
template < typename TInputImage0, typename TInputImage1, typename TMaskImage, typename TOutputImage > 
void
MakeCostImageFilter< TInputImage0, TInputImage1, TMaskImage, TOutputImage > 
::SetBlackBloodImage( typename MakeCostImageFilter< TInputImage0, TInputImage1, TMaskImage, TOutputImage >::BBImageType * bbImage ){
  this->ProcessObject::SetNthInput( 0, const_cast< BBImageType * >( bbImage ) );
}

template < typename TInputImage0, typename TInputImage1, typename TMaskImage, typename TOutputImage > 
typename MakeCostImageFilter< TInputImage0, TInputImage1, TMaskImage, TOutputImage >::BBImageType* 
MakeCostImageFilter< TInputImage0, TInputImage1, TMaskImage, TOutputImage > 
::GetBlackBloodImage() {
  return ( static_cast< BBImageType* >(this->ProcessObject::GetInput(0)) );
}

// Input 1
template < typename TInputImage0, typename TInputImage1, typename TMaskImage, typename TOutputImage > 
void
MakeCostImageFilter< TInputImage0, TInputImage1, TMaskImage, TOutputImage > 
::SetPhaseContrastImage( typename MakeCostImageFilter< TInputImage0, TInputImage1, TMaskImage, TOutputImage >::PCImageType * pcImage ){
  this->ProcessObject::SetNthInput( 1, const_cast< PCImageType * >( pcImage ) );
}

template < typename TInputImage0, typename TInputImage1, typename TMaskImage, typename TOutputImage > 
typename MakeCostImageFilter< TInputImage0, TInputImage1, TMaskImage, TOutputImage >::PCImageType* 
MakeCostImageFilter< TInputImage0, TInputImage1, TMaskImage, TOutputImage > 
::GetPhaseContrastImage() {
  return ( static_cast< PCImageType* >(this->ProcessObject::GetInput(1)) );
}

// Input 2
template < typename TInputImage0, typename TInputImage1, typename TMaskImage, typename TOutputImage > 
void
MakeCostImageFilter< TInputImage0, TInputImage1, TMaskImage, TOutputImage > 
::SetMaskImage( typename MakeCostImageFilter< TInputImage0, TInputImage1, TMaskImage, TOutputImage >::MaskImageType * maskImage ){
  this->ProcessObject::SetNthInput( 2, const_cast< MaskImageType * >( maskImage ) );
}

template < typename TInputImage0, typename TInputImage1, typename TMaskImage, typename TOutputImage > 
typename MakeCostImageFilter< TInputImage0, TInputImage1, TMaskImage, TOutputImage >::MaskImageType* 
MakeCostImageFilter< TInputImage0, TInputImage1, TMaskImage, TOutputImage > 
::GetMaskImage() {
  return ( static_cast< MaskImageType* >(this->ProcessObject::GetInput(2)) );
}

template < typename TInputImage0, typename TInputImage1, typename TMaskImage, typename TOutputImage > 
void
MakeCostImageFilter< TInputImage0, TInputImage1, TMaskImage, TOutputImage > 
::SetSeedList( typename MakeCostImageFilter< TInputImage0, TInputImage1, TMaskImage, TOutputImage >::SeedListType seedList ){
  m_SeedList = seedList;
}

template < typename TInputImage0, typename TInputImage1, typename TMaskImage, typename TOutputImage > 
typename MakeCostImageFilter< TInputImage0, TInputImage1, TMaskImage, TOutputImage >::SeedListType
MakeCostImageFilter< TInputImage0, TInputImage1, TMaskImage, TOutputImage > 
::GetSeedList() {
  return m_SeedList;
}

template < typename TInputImage0, typename TInputImage1, typename TMaskImage, typename TOutputImage > 
void
MakeCostImageFilter< TInputImage0, TInputImage1, TMaskImage, TOutputImage > 
::SetRadii( typename MakeCostImageFilter< TInputImage0, TInputImage1, TMaskImage, TOutputImage >::RadiusListType radii ){
  m_RadiusList = radii;
}

template < typename TInputImage0, typename TInputImage1, typename TMaskImage, typename TOutputImage > 
typename MakeCostImageFilter< TInputImage0, TInputImage1, TMaskImage, TOutputImage >::RadiusListType
MakeCostImageFilter< TInputImage0, TInputImage1, TMaskImage, TOutputImage > 
::GetRadii() {
  return m_RadiusList;
}

template < typename TInputImage0, typename TInputImage1, typename TMaskImage, typename TOutputImage > 
void
MakeCostImageFilter< TInputImage0, TInputImage1, TMaskImage, TOutputImage > 
::UpdateParameters()
{

  if ( m_UseOnlyBBImage ) {
    m_InvertFilter->SetInput( m_BBMedialnessFilter->GetOutput() );
  } else {
    m_MedMaxFilter->SetInput(0, m_BBMedialnessFilter->GetOutput() );
    m_MedMaxFilter->SetInput(1, m_PCMedialnessFilter->GetOutput() );
    m_LisMaxFilter->SetInput(0, m_BBLisFilter->GetOutput() );
    m_LisMaxFilter->SetInput(1, m_PCLisFilter->GetOutput() );
    m_MultiplyFilter->SetInput(0, m_MedMaxFilter->GetOutput());
    m_MultiplyFilter->SetInput(1, m_LisMaxFilter->GetOutput());
    m_InvertFilter->SetInput( m_MultiplyFilter->GetOutput() );
  }
  m_BBMedialnessFilter->SetVerbose( m_Verbose );
  m_PCMedialnessFilter->SetVerbose( m_Verbose );

  m_BBLisFilter->SetVerbose( m_Verbose );
  m_BBLisFilter->SetSeedList( m_SeedList );
  m_PCLisFilter->SetVerbose( m_Verbose );
  m_PCLisFilter->SetSeedList( m_SeedList );
}

template < typename TInputImage0, typename TInputImage1, typename TMaskImage, typename TOutputImage > 
void
MakeCostImageFilter< TInputImage0, TInputImage1, TMaskImage, TOutputImage > 
::GenerateData()
{

  if ( m_Verbose ) {
    std::cout << "Start calculating cost image" << std::endl;
    std::cout << "number of seed points: " << m_SeedList.size() << std::endl;
    std::cout << "number of radii: " << m_RadiusList.size() << std::endl;
  }

  typename BBImageType::Pointer bbImage = this->GetBlackBloodImage();
  typename PCImageType::Pointer pcImage;
  if ( !m_UseOnlyBBImage ) {
    pcImage = this->GetPhaseContrastImage();
  }
  typename MaskImageType::Pointer   mask;

  m_BBMedialnessFilter->SetUseMask( m_UseMask );
  m_PCMedialnessFilter->SetUseMask( m_UseMask);
  if (m_UseMask) {
    mask = this->GetMaskImage();
    m_BBMedialnessFilter->SetMaskImage( mask );
    m_PCMedialnessFilter->SetMaskImage( mask );
  }

  m_BBMedialnessFilter->SetInput( bbImage );
  typedef itk::ImageFileWriter< PrecisionImageType > WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( "bbMedImage.mhd");
  writer->SetInput( m_BBMedialnessFilter->GetOutput() );
  writer->Update();

  if ( !m_UseOnlyBBImage ) {
    m_PCMedialnessFilter->SetInput( pcImage );
    writer->SetFileName( "pcMedImage.mhd");
    writer->SetInput( m_PCMedialnessFilter->GetOutput() );
    writer->Update();
  }

  m_BBLisFilter->SetInput( bbImage );
  writer->SetFileName( "bbLisImage.mhd");
  writer->SetInput( m_BBLisFilter->GetOutput() );
  writer->Update();

  if ( !m_UseOnlyBBImage ) {
    m_PCLisFilter->SetInput( pcImage );
    writer->SetFileName( "pcLisImage.mhd");
    writer->SetInput( m_PCLisFilter->GetOutput() );
    writer->Update();
  }

  m_InvertFilter->GraftOutput( this->GetOutput() );
  m_InvertFilter->Update();
  this->GraftOutput( m_InvertFilter->GetOutput() );

  if ( m_Verbose ) {
    std::cout << "End calculating cost image" << std::endl;
  }
}

} // namespace itk
#endif
