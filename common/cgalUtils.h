#ifndef _cgalutils_h
#define _cgalutils_h

#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"

// Include conversion options
#include "VTKPolyDataToCgalPolyhedron.h"
#include "ObjToCgalPolyhedron.h"

#include "boost/filesystem.hpp"
namespace fs = boost::filesystem;

void loadPoly( const std::string & filename, Polyhedron & poly )
{
  // Read surface from input file
  fs::path p ( filename );
  if ( p.extension() == ".vtk" ) 
  {
    // Read VTK polydata object
    vtkPolyDataReader * reader;
    vtkPolyData * polyVTK;
    reader = vtkPolyDataReader::New();
    reader->SetFileName( filename.c_str() );
    reader->Update();
    polyVTK = reader->GetOutput();
    polyVTK->Update();
    std::cout << "vtkPolyData [ " << polyVTK->GetNumberOfPoints() << " ] " << std::endl;

    // Convert VTK polydata to polyhedron
    VTKPolyDataToOFF( polyVTK, poly );

    // Delete reader
    reader->Delete();
  }
  else if ( p.extension() == ".off" )
  {
    std::ifstream is ( filename.c_str() );
    is >> poly;
  }
  else if ( p.extension() == ".obj" )
  {
    ObjToOFF( filename, poly );
  }
  else
  {
    std::cerr << "Unknown file extension " << p.extension() << "!" << std::endl;
    exit( -1 );
  }
}

#endif // _cgalutils_h

