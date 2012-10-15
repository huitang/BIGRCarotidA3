#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>

//#include <CGAL/Polyhedron_3.h>

double to_double( const float a )
{
  return static_cast< double >( a );
}

//template< class Kernel >
void cgalPolyhedronToVTKPolyData( Polyhedron polyIn, vtkPolyData * polyOut ) {
  //typedef typename CGAL::Polyhedron_3<Kernel> PolyType;

  // Reset points and polydata but do not reallocate.
  polyOut->DeleteCells();
  if (polyOut->GetPointData()){ polyOut->GetPointData()->Reset(); }
  
  vtkPoints * points = vtkPoints::New();
 
  // Table which stores node ids for entry index.
  std::vector<unsigned int> idsMap ( polyIn.size_of_vertices() );
  std::vector< vtkIdList * > facetsmap;
  std::map< Kernel::Point_3 *, vtkIdType > vertexmap; 
  std::vector< vtkIdType > ids;

  for ( Polyhedron::Facet_iterator it = polyIn.facets_begin(); it != polyIn.facets_end(); ++it) {
    facetsmap.push_back( vtkIdList::New() );

    typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_around_facet_circulator;
    Halfedge_around_facet_circulator he, end;
    he = end = it->facet_begin();
    CGAL_For_all(he,end)
    {
      Kernel::Point_3 point = he->vertex()->point();
      
      float pc[3];
      pc[0] = to_double( point.x() );
      pc[1] = to_double( point.y() );
      pc[2] = to_double( point.z() );

      std::map< Kernel::Point_3 *, vtkIdType >::iterator pointInMap = vertexmap.find( &he->vertex()->point() );
      
      vtkIdType pid;
      if ( pointInMap == vertexmap.end() ) {
        // Add coordinate to point list.
        pid = points->InsertNextPoint( pc );
        vertexmap.insert( std::make_pair( &he->vertex()->point(), pid ) );
      } else {
        pid = pointInMap->second;
      }
      
      facetsmap[ facetsmap.size()-1 ]->InsertNextId( pid );
    }
  }
  
  // Allocate internal data for poly data.
  polyOut->Allocate();

  // Specify points for the polydata mesh.
  polyOut->SetPoints(points);

  // Build mesh consistency.
  polyOut->BuildLinks();
  polyOut->BuildCells();

  for ( std::vector< vtkIdList * >::iterator it = facetsmap.begin(); it != facetsmap.end(); ++it) {
    polyOut->InsertNextCell(VTK_POLYGON, *it);
  }

  // Build mesh consistency.
  polyOut->BuildLinks();
  polyOut->BuildCells();
}
