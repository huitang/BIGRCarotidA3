//
//
#ifndef __MEDIALNESSIMAGEFILTER__
#define __MEDIALNESSIMAGEFILTER__

#include "itkImageSource.h"
#include "itkConceptChecking.h"
#include "itkImageToImageFilterDetail.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <windows.h>
#include "vtksys/CommandLineArguments.hxx"
#include "vtkMetaImageReader.h"
#include "vtkImageData.h"
#include "vtkMetaImageWriter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkNthElementImageAdaptor.h"
#include "itkIndex.h"
#include "itkImageAdaptor.h"
#include "itkImageRegionIterator.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkCastImageFilter.h"


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

        typedef TOutputImage                             OutputImageType;
       // typedef typename OutputImageType::Pointer        OutputImagePointer;
        typedef typename OutputImageType::PixelType      OutputPixelType;
        typedef typename OutputImageType::RegionType     OutputRegionType;
        typedef typename OutputImageType::SizeType       OutputSizeType;
        typedef typename OutputImageType::SizeValueType  OutputSizeValueType;
        typedef typename OutputImageType::IndexType      OutputIndexType;
        typedef typename OutputImageType::IndexValueType OutputIndexValueType;
        //typedef typename OutputImageType::Dimension   Dimension;

        void SetNumberOfAngles(int a)
        {
            m_a=a;
        }

        int GetNumberOfAngles()
        {
            return m_a;
        }

        void SetSigma(float s)
        {
            m_sigma=s;
        }

        float GetSigma()
        {
            return m_sigma;
        }

        void SetMinRadius(float r)
        {
            m_r=r;
        }

        float GetMinRadius()
        {
            return m_r;
        }

        void SetMaxRadius(float R)
        {
            m_R=R;
        }

        float GetMaxRadius()
        {
            return m_R;
        }

        void SetNumberOfRadius(int n)
        {
            m_n=n;
        }

        int GetNumberOfRadius()
        {
            return m_n;
        }

        void SetDarkLumen(bool d)
        {
            m_darkLumen=d;
        }

        bool GetDarkLumen()
        {
            return m_darkLumen;
        }

        void SetMaskImage( MaskImageType * maskImage)
        {

        this->ProcessObject::SetNthInput( 1, const_cast< MaskImageType * >( maskImage ) );
        }


      MaskImageType * GetMaskImage()
        {
            return ( static_cast< MaskImageType* >( this->ProcessObject::GetInput(1) ) );
        }

    protected:
         void operator=( const Self& ); 
        int m_a;
        bool m_darkLumen;
        float m_sigma;
        float m_r;
        float m_R;
        int m_n;
        MedialnessImageFilter(const Self &);

    /** Constructor 
     * Note that it is declared protected. This makes sure that 
     * users can only instantiate this class with the ::New() method. */
    MedialnessImageFilter()
    {
        ;
    }

    /** Destructor */
    virtual ~MedialnessImageFilter()
    {
        ;
    }

    /** GenerateData
     * This method does the work. It is automatically invoked
     * when the user calls Update(); */
    virtual void GenerateData(void)
    {
        
        //allocate memory
        typename OutputImageType::Pointer output = this->GetOutput();
        output->SetBufferedRegion( output->GetRequestedRegion() );
        output->Allocate();
        typedef typename InputImageType::SizeType InputSizeType;
        InputSizeType size = this->GetOutput()->GetLargestPossibleRegion().GetSize(); 
        const unsigned int Dimension = 3;
        //calculate gradient of input image
        typedef float ComponentType;
        typedef itk::CovariantVector< ComponentType,Dimension> OutGradientVectorPixelType;
        typedef itk::Image <OutGradientVectorPixelType,Dimension> OutGradientVectorImageType;
        typedef itk::GradientRecursiveGaussianImageFilter<InputImageType,OutGradientVectorImageType> GradientFilterType;
        GradientFilterType::Pointer gradientFilter= GradientFilterType::New();
        gradientFilter->SetSigma(this->m_sigma);
        gradientFilter->SetInput(this->GetInput());
        gradientFilter->UpdateLargestPossibleRegion();
        cout<<"Gradient Calculation Done"<<endl;
        OutGradientVectorImageType::Pointer vectorImage = gradientFilter->GetOutput();
        typedef itk::VectorIndexSelectionCastImageFilter<OutGradientVectorImageType, InputImageType> IndexSelectionFilterType;
        IndexSelectionFilterType::Pointer indexSelectionFilterX = IndexSelectionFilterType::New();
        indexSelectionFilterX->SetInput(vectorImage);
        indexSelectionFilterX->SetIndex(0);
        indexSelectionFilterX->UpdateLargestPossibleRegion();
        InputImageType::Pointer gradientx=indexSelectionFilterX->GetOutput();
        IndexSelectionFilterType::Pointer indexSelectionFilterY = IndexSelectionFilterType::New();
        indexSelectionFilterY->SetInput(vectorImage);
        indexSelectionFilterY->SetIndex(1);
        indexSelectionFilterY->UpdateLargestPossibleRegion();
        InputImageType::Pointer gradienty=indexSelectionFilterY->GetOutput();
              
        //get image spacing and size
        InputImageType::IndexType idx;
        InputImageType::IndexType neibourIdx;
        std::vector<float> imgExt(3, 0);
        imgExt[0] = this->GetInput()->GetBufferedRegion().GetSize()[0];
        imgExt[1] = this->GetInput()->GetBufferedRegion().GetSize()[1];
        imgExt[2] = this->GetInput()->GetBufferedRegion().GetSize()[2];
        const InputImageType::SpacingType& sp = this->GetOutput()->GetSpacing();
        std::vector<float> voxelSizes (3, 1.0f);
        voxelSizes[0] = sp[0];
        voxelSizes[1] = sp[1];
        voxelSizes[2] = sp[2];

        //Region Iterator
        typedef itk::ImageRegionIteratorWithIndex<InputImageType> ImageIterator;
        ImageIterator it0( gradientx, gradientx->GetRequestedRegion() );
        it0.GoToBegin(); 
        ImageIterator it1( gradienty, gradienty->GetRequestedRegion() );
        it1.GoToBegin();  
        ImageIterator it2( this->GetMaskImage(),this->GetMaskImage()->GetRequestedRegion() );
        it2.GoToBegin();
        ImageIterator itOut( this->GetOutput(), this->GetOutput()->GetRequestedRegion() );
        itOut.GoToBegin();
        //calculatemedialness
 
        for (int zi = 0; zi<imgExt[2];  zi++){
            for (int yi = 0;  yi<imgExt[1];  yi++){
                for (int xi = 0; xi<imgExt[0];  ++xi){  
                    idx[0]=xi;
                    idx[1]=yi;
                    idx[2]=zi;
                    it0.SetIndex(idx);
                    it1.SetIndex(idx);
                    it2.SetIndex(idx);
                    if (it2.Get())
                    {
                        std::vector<std::vector<float>> allBoundaryMeasure;
                        std::vector<std::vector<float>> allIntensityProfile;
                        std::vector<float> maxBoundaryMeasure(this->m_a);
                        std::vector<int> maxBoundaryMeasureIndex(this->m_a);
                        float stepOfAngle =2*3.1415926/this->m_a;
                        for (int i=0; i<this->m_a;i++)
                        {
                            std::vector<float> boundaryMeasureAtiAngle;
                            float deltaRadius= (this->m_R-this->m_r)/(this->m_n-1);
                            for (int ri=0;ri<this->m_n;ri++)                                      
                            {    
                                //vec3 currentSV(x,y,z);
                                //vec3 currentSW;
                                //getInImg(0)->transformToWorldCoord(currentSV,currentSW);
                                //float xW = currentSW[0] + (r+r*deltaRadius)*cos(i*stepOfAngle);
                                //float yW = currentSW[1] + (r+r*deltaRadius)*sin(i*stepOfAngle);
                                //float zW = currentSW[2];  

                                //vec3 currentRW(xW,yW,zW);
                                //vec3 currentRV;
                                //getInImg(0)->transformToVoxelCoord(currentRW,currentRV);

                                int x=min(max(xi+((ri+ri*deltaRadius)*cos(i*stepOfAngle)/voxelSizes[0]),(float)0),(float)(imgExt[0]-1));
                                int y=min(max(yi+((ri+ri*deltaRadius)*sin(i*stepOfAngle)/voxelSizes[1]),(float)0),(float)(imgExt[1]-1));
                                int z=zi;
                                neibourIdx[0]=(int)x;
                                neibourIdx[1]=(int)y;
                                neibourIdx[2]=(int)z;
                                it0.SetIndex(neibourIdx);
                                it1.SetIndex(neibourIdx);
                                if (this->m_darkLumen)
                                {
                                    boundaryMeasureAtiAngle.push_back((float)max(it0.Get()*cos(i*stepOfAngle)+it1.Get()*sin(i*stepOfAngle),(float)0));

                                }
                                else
                                { 
                                    boundaryMeasureAtiAngle.push_back((float)min(it0.Get()*cos(i*stepOfAngle)+it1.Get() *sin(i*stepOfAngle),(float)0));
                                }

                                if (abs((float)(it0.Get()*cos(i*stepOfAngle)+it1.Get()*sin(i*stepOfAngle)))>maxBoundaryMeasure[i])
                                {
                                    maxBoundaryMeasure[i]=abs((float)(it0.Get() *cos(i*stepOfAngle)+it1.Get() *sin(i*stepOfAngle)));
                                    maxBoundaryMeasureIndex[i]=ri;
                                }


                            }
                            allBoundaryMeasure.push_back(boundaryMeasureAtiAngle);

                        }

                        std::vector<std::vector<float>> allBoundaryMeasureN;     
                        for (int i=0; i<this->m_a;i++)
                        {
                            std::vector<float> boundaryMeasureAtiAngleN;
                            for (int r=0;r<this->m_n;r++)                                      
                            {  
                                std::vector<float> currentBoundaryMeasureFrom0ToR;
                                float currentMax=maxBoundaryMeasure[i];
                                if (abs(currentMax)>0.005)
                                {
                                    boundaryMeasureAtiAngleN.push_back((float)abs((float)allBoundaryMeasure[i][r])/abs(currentMax));
                                    if ((float)abs((float)allBoundaryMeasure[i][r])/abs(currentMax)>1)
                                    {
                                        std::cout<<abs(currentMax)<<std::endl;
                                    }
                                }
                                else
                                {
                                    boundaryMeasureAtiAngleN.push_back((float)0);
                                }

                            }

                            allBoundaryMeasureN.push_back(boundaryMeasureAtiAngleN);
                        }
                        float maxMedialness=0;
                        for (int ri=0;ri<this->m_n;ri++)                                      
                        {   
                            float currentMedialness=0;
                            for (int i=0; i<this->m_a;i++)
                            {  
                                currentMedialness=currentMedialness+allBoundaryMeasureN[i][ri];
                            }
                            if (currentMedialness>maxMedialness)
                            {
                                maxMedialness=currentMedialness;
                            }
                        }	
                        itOut.SetIndex(idx);
                        itOut.Set((OutputImageType::PixelType)maxMedialness/m_a);
                    }
                    else
                    {
                        itOut.SetIndex(idx);
                        itOut.Set(0);
                    }
                }
            }
        }
    }
};
}

#endif