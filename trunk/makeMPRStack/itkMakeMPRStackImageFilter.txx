//
//

#ifndef __MakeMPRStackImageFilter__TXX
#define __MakeMPRStackImageFilter__TXX

#include "boost/program_options.hpp"
//#include "boost/mpl/size_t_fwd.hpp" 
#include "mlAPI.h"
#include "mlUtilsSystem.h"
#include "mlInitSystemML.h"
#include "mlDataTypes.h"
#include "mlModuleIncludes.h"
#include "mlITKSupportToolFunctions.h"
#include "General/Sources/ML/EMCUtilities/mlPortSubImageToOutput.h"
// Include path classes
#include "pathTrackingLibrary.h"
#include "SeedPointFileIO.h"

#include "itkImageFileWriter.h" 

namespace itk {


template < typename TInputImage, typename TOutputImage > 
MakeMPRStackImageFilter< TInputImage, TOutputImage > 
::MakeMPRStackImageFilter(){

}



template < typename TInputImage, typename TOutputImage > 
void
MakeMPRStackImageFilter< TInputImage, TOutputImage > 
::SetPath( typename MakeMPRStackImageFilter< TInputImage, TOutputImage >::PathType path ){
  m_Path = path;
}

template < typename TInputImage, typename TOutputImage > 
typename MakeMPRStackImageFilter< TInputImage, TOutputImage >::PathType
MakeMPRStackImageFilter< TInputImage, TOutputImage > 
::GetPath() {
  return m_Path;
}


template < typename TInputImage, typename TOutputImage > 
void
MakeMPRStackImageFilter< TInputImage, TOutputImage > 
::GenerateData()
{

  if ( m_Verbose ) {
    std::cout << "Start constructing MPR stack" << std::endl;
    std::cout << "number of path points: " << m_Path.size() << std::endl;
  }

  typename InputImageType::ConstPointer  input = ( this->GetInput() );
  typename OutputImageType::Pointer output = this->GetOutput();

  // ROI
  InputImageType::RegionType desiredRegion;
  // Start of ROI
  InputImageType::IndexType start;
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  desiredRegion.SetIndex( start );

  // Create vector of image extents
  std::vector<int> imgExt(3, 0);

  // Intensity
  typedef ExtractImageFilter< InputImageType, InputImageType > FilterType;
  FilterType::Pointer filter0 = FilterType::New();
  InputImageType::SizeType inputSize = input->GetLargestPossibleRegion().GetSize();
  desiredRegion.SetSize( inputSize );

  // Get ROI from image
  filter0->SetExtractionRegion( desiredRegion );
  filter0->SetInput( input );
  filter0->Update();
  filter0->GetOutput()->GetPixelContainer()->SetContainerManageMemory(false);

  // Set pointer to data
  float* data = NULL;
  data = filter0->GetOutput()->GetPixelContainer()->GetImportPointer();

  // Init ML
  std::cout << std::endl << "Initializing ML" << std::endl;
  MLInit(ML_MAJOR_VERSION, ML_MAJOR_CAPI_VERSION, ML_CAPI_REVISION);
  std::cout << "ML initialized" << std::endl;
  char  buffer[4096]="\n";

  // Create ML image and import it into ML
  std::cout << std::endl<< "create ml image" << std::endl;
  ::ml::SubImage *mlImage = new ::ml::SubImage;
  mlImage->setDataType( MLfloatType );
  std::cout << "image size" << inputSize << std::endl;
  ::ml::ImageVector imageSize(inputSize[0],inputSize[1],inputSize[2],1,1,1);
  mlImage->setImgExt( imageSize );
  std::cout << "set data pointer" << std::endl;
  mlImage->setData( data );
  if ( mlImage->isValid() ) {
    std::cout << "mlimage valid:" <<std::endl;
  } else {
    std::cout << "mlimage not valid:" <<std::endl;
  }
  std::cout << "ml image created" << std::endl;

  ::ml::MedicalImageProperties *props = new ::ml::MedicalImageProperties;
  ::ml::setMLWorldFromITKScaleOriginAndOrientation<InputImageType> (filter0->GetOutput(), *props, true);

  // Load dlls
  MLErrorCode libErr = MLLoadLibraryWOError("EMCUtilities");
  libErr = 10*MLLoadLibraryWOError("MLResample1");
  libErr += 100*MLLoadLibraryWOError("MLImageFile");
  libErr += 1000*MLLoadLibraryWOError("MLDicomTree_OFFIS");
  libErr += 10000*MLLoadLibraryWOError("MLParser");
  if ( libErr != ML_RESULT_OK ) {
    std::cout << "Error loading dll " << libErr << std::endl;
    //return;
  }

  // Create modules
  std::cout << "Create modules" << std::endl;
  mlModule* import = MLCreateModuleFromName("PortSubImageToOutput");
  mlModule* ptkf   = MLCreateModuleFromName("PathToKeyFrame");
  mlModule* mpr    = MLCreateModuleFromName("MPRPathVis");
  mlModule* writer = MLCreateModuleFromName("ImgSave");
  mlModule* sbase  = MLCreateModuleFromName("SaveBase");
  if ( !import || !ptkf || !mpr || !writer ) {
    std::cout << "Error creating module" << std::endl;
    return;
  }

  // Get fields
  std::cout << "get fields" << std::endl;
  mlField* importerOutput = MLModuleGetField( import, "output0");
  mlField* ptkfInput      = MLModuleGetField( ptkf,   "inputKeys");
  mlField* ptkfOutput     = MLModuleGetField( ptkf,   "outputKeys");
  mlField* mprnKeys       = MLModuleGetField( mpr,    "maxKeyFrame");
  mlField* mprImageInput  = MLModuleGetField( mpr,    "input0");
  mlField* mprPathInput   = MLModuleGetField( mpr,    "frames");
  mlField* mprStackOutput = MLModuleGetField( mpr,    "output1");
  mlField* writerInput0   = MLModuleGetField( writer, "input0");
  mlField* saveField      = MLModuleGetField( writer, "save");
  mlField* statusField    = MLModuleGetField( writer, "status");
  mlField* writerFileName = MLModuleGetField( writer, "filename");
  if ( !importerOutput || !ptkfInput      || !ptkfOutput   || 
       !mprnKeys       || !mprImageInput  || !mprPathInput || 
       !mprStackOutput || !writerInput0   || !saveField    || 
       !statusField    || !writerFileName ){
    std::cout << "could not find field: " 
              << importerOutput << ptkfInput     << ptkfOutput 
              << mprnKeys       << mprImageInput << mprPathInput 
              << mprStackOutput << writerInput0  << saveField 
              << statusField    << writerFileName 
              << std::endl;
    return;
  }

  // Connect modules
  std::cout << std::endl << "connect modules " << std::endl;
  MLFieldConnectFrom(mprImageInput,importerOutput);
  MLFieldConnectFrom(mprPathInput,ptkfOutput);
  MLFieldConnectFrom(writerInput0,mprStackOutput);
  std::cout << "modules connected" << std::endl;

  std::cout << "import ml image" << std::endl;
  ::ml::PortSubImageToOutput* importModule = static_cast< ::ml::PortSubImageToOutput* >(import);
  importModule->setImageProperties( props );
  importModule->setSubImage( mlImage );


  std::cout << "Create marker list and set to input of PathToKeyFrame module" << std::endl;
  ::ml::XMarkerList markerListPath;
  PathType::iterator pIt = m_Path.begin();
  for (;pIt < m_Path.end();++pIt){
    markerListPath.push_back( ::ml::XMarker(::ml::Vector3( (*pIt)[0],(*pIt)[1],(*pIt)[2] ) ) );
  }
  ::ml::BaseField* markerListInput = static_cast< ::ml::BaseField* >( ptkfInput );
  markerListInput->setBaseValue ( &markerListPath );
  MLFieldTouch(ptkfInput);
  MLFieldTouch( mprImageInput ); 
  MLFieldTouch( importerOutput );

  MLFieldGetValue(mprnKeys, buffer, 4096);
  std::cout << "Number of Keys: " << buffer << std::endl;

  
  ml::PagedImage* outputImage = static_cast<ml::PagedImage*>(MLModuleGetField( mpr, "output0"));
  std::cout << "Image valid: " << outputImage->isValid() << std::endl;

  // Set FoV
  ml::FloatField* fOVField = static_cast< ml::FloatField* >( MLModuleGetField (mpr, "fieldOfView") );
  fOVField->setFloatValue( static_cast< float>(30.0f) );

  // Set fillValue
  ml::DoubleField* fillValue = static_cast< ml::DoubleField* >( MLModuleGetField (mpr, "fillValue") );
  fillValue->setDoubleValue( static_cast< double >(10000.0) );

  // Get the file name field of the writer.
  MLModuleSetFieldValue( writer, "format",   "DICOMTIFF" );
  MLModuleSetFieldValue( writer, "filename", "d:/testimage.dcm" );
  MLFieldGetValue(writerFileName, buffer, 4096);
  std::cout << "Save file name: " << buffer << std::endl;
  std::cout << "save image" << std::endl;
  MLFieldTouch( saveField );
  MLFieldGetValue( statusField, buffer, 4096);
  std::cout << "Write status: " << buffer << std::endl;

  // Free allocated memory
  delete[] data;
  delete mlImage;
  delete props;

  if ( m_Verbose ) {
    std::cout << "End constructing MPR image" << std::endl;
  }
}

} // namespace itk
#endif
