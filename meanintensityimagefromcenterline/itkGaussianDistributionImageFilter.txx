//
//

#ifndef __GaussianDistributionImageFilter__TXX
#define __GaussianDistributionImageFilter__TXX

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

namespace itk {


template< class TInputImage, class TOutputImage > 
GaussianDistributionImageFilter< TInputImage, TOutputImage >
::GaussianDistributionImageFilter(){
  m_Mean = .0;
  m_Variance = 1.0;
  m_Gaussian = GaussianType::New();
  m_Gaussian->SetMean( m_Mean );
  m_Gaussian->SetVariance( m_Variance );
}

template< class TInputImage, class TOutputImage > 
void GaussianDistributionImageFilter< TInputImage, TOutputImage >
::SetMean( double mean ){
  m_Mean = mean;
  m_Gaussian->SetMean( m_Mean );
}

template< class TInputImage, class TOutputImage > 
double GaussianDistributionImageFilter< TInputImage, TOutputImage >
::GetMean(){
  return m_Mean;
}

template< class TInputImage, class TOutputImage > 
void GaussianDistributionImageFilter< TInputImage, TOutputImage >
::SetVariance( double var ){
  m_Variance = var;
  m_Gaussian->SetVariance( m_Variance );
}

template< class TInputImage, class TOutputImage > 
double GaussianDistributionImageFilter< TInputImage, TOutputImage >
::GetVariance(){
  return m_Mean;
}

template< class TInputImage, class TOutputImage > 
void GaussianDistributionImageFilter< TInputImage, TOutputImage >
::SetSigma( double sig ){
  m_Variance = sig*sig;
  m_Gaussian->SetVariance( m_Variance );
}

template< class TInputImage, class TOutputImage > 
double GaussianDistributionImageFilter< TInputImage, TOutputImage >
::GetSigma(){
  return sqrt( m_Variance );
}

template< class TInputImage, class TOutputImage > 
void GaussianDistributionImageFilter< TInputImage, TOutputImage >
::BeforeThreadedGenerateData(){
}

template< class TInputImage, class TOutputImage > 
void GaussianDistributionImageFilter< TInputImage, TOutputImage >
::ThreadedGenerateData(const typename ImageToImageFilter<TInputImage,TOutputImage>::OutputImageRegionType& outputRegionForThread, int threadId)
{

  //InputImageType::Pointer input = this->GetInput();
  typedef typename itk::ImageRegionConstIterator< InputImageType > ConstIteratorType;
  ConstIteratorType inputIt(   this->GetInput(),  outputRegionForThread  );

  //OutputImageType::Pointer output = this->GetOutput();
  typedef typename itk::ImageRegionIterator< OutputImageType >       IteratorType;
  IteratorType      outputIt(  this->GetOutput(), outputRegionForThread );


  inputIt.GoToBegin();
  outputIt.GoToBegin();

  while( !inputIt.IsAtEnd() )
  {
    outputIt.Set(  static_cast<OutputPixelType>( m_Gaussian->EvaluatePDF( static_cast<double>( inputIt.Get() )))  );
    ++inputIt;
    ++outputIt;
  }
}

} // namespace itk
#endif
