#ifndef __mlITKTools_H
#define __mlITKTools_H

#include "itkImportImageFilter.h"
#include "mlPagedImage.h"
#include "mlMatrix4.h"
#include "EMCUtilitiesSystem.h"

template <typename ITK_DATATYPE, unsigned int DIM>
static typename itk::ImportImageFilter<ITK_DATATYPE, DIM>::Pointer mlToITKImage(const ml::PagedImage *mlImage, void* data )
{
  typedef itk::Image< ITK_DATATYPE, DIM > ImageType;

  typename ImageType::SizeType  size;
  size[0] = mlImage->getImgExt()[0];  // size along X
  size[1] = mlImage->getImgExt()[1];  // size along Y
  size[2] = mlImage->getImgExt()[2];  // size along Z

  typename ImageType::IndexType start;
  start[0] = 0;  // first index on X
  start[1] = 0;  // first index on Y
  start[2] = 0;  // first index on Z

  typename ImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );

  double origin[ DIM ];
  const ml::Matrix4 *worldMatrix = mlImage->getToWorldMatrix(); 
  origin[0] = (*worldMatrix)[0][3];
  origin[1] = (*worldMatrix)[1][3];
  origin[2] = (*worldMatrix)[2][3];

  double spacing[ DIM ];
  spacing[0] = mlImage->getVoxelSize()[0];
  spacing[1] = mlImage->getVoxelSize()[1];
  spacing[2] = mlImage->getVoxelSize()[2];

  int numberOfVoxels = mlImage->getSize();

  // Set up the itk import filters
  typedef itk::ImportImageFilter
    < ITK_DATATYPE, 
    DIM > ImportFilterType;

  typename ImportFilterType::Pointer imageImporter = ImportFilterType::New();

  imageImporter->SetRegion(  region  );
  imageImporter->SetOrigin(  origin  );
  imageImporter->SetSpacing( spacing );
  imageImporter->SetImportPointer( (ITK_DATATYPE*)data, 
    numberOfVoxels,
    false );
  imageImporter->Update();  
  return imageImporter;
}

template <typename ITK_DATATYPE, unsigned int DIM>
static typename itk::ImportImageFilter<ITK_DATATYPE, DIM>::Pointer itkToMLImage(itk::Image<ITK_DATATYPE, DIM>::Pointer itkImage, ml::PagedImage *mlImage)
{
  typedef itk::Image< ITK_DATATYPE, DIM > ImageType;

  typename ImageType::SizeType  size;
  size[0] = mlImage->getImgExt()[0];  // size along X
  size[1] = mlImage->getImgExt()[1];  // size along Y
  size[2] = mlImage->getImgExt()[2];  // size along Z

  typename ImageType::IndexType start;
  start[0] = 0;  // first index on X
  start[1] = 0;  // first index on Y
  start[2] = 0;  // first index on Z

  typename ImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );

  double origin[ DIM ];
  const ml::Matrix4 *worldMatrix = mlImage->getToWorldMatrix(); 
  origin[0] = (*worldMatrix)[0][3];
  origin[1] = (*worldMatrix)[1][3];
  origin[2] = (*worldMatrix)[2][3];

  double spacing[ DIM ];
  spacing[0] = mlImage->getVoxelSize()[0];
  spacing[1] = mlImage->getVoxelSize()[1];
  spacing[2] = mlImage->getVoxelSize()[2];

  int numberOfVoxels = mlImage->getSize();

  // Set up the itk import filters
  typedef itk::ImportImageFilter
    < ITK_DATATYPE, 
    DIM > ImportFilterType;

  typename ImportFilterType::Pointer imageImporter = ImportFilterType::New();

  imageImporter->SetRegion(  region  );
  imageImporter->SetOrigin(  origin  );
  imageImporter->SetSpacing( spacing );
  imageImporter->SetImportPointer( (ITK_DATATYPE*)data, 
    numberOfVoxels,
    false );
  imageImporter->Update();  
  return imageImporter;
}

#endif
