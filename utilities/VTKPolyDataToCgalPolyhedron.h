#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkCell.h>

#include <CGAL/Polyhedron_3.h>

typedef Polyhedron::HalfedgeDS HalfedgeDS;

// A modifier creating a Polyhedron from a VTK polydata object.
template <class HDS>
class Build_polygon : public CGAL::Modifier_base<HDS> {
public:
  Build_polygon( vtkPolyData * poly  ) {
    _poly = poly;
  }

  void operator()( HDS& hds) {
    CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
    B.begin_surface( _poly->GetNumberOfVerts(), _poly->GetNumberOfPolys() );
    typedef typename HDS::Vertex Vertex;
    typedef typename Vertex::Point Point;
        
    vtkPoints * allpoints = _poly->GetPoints();

    int nodeid = 0;
    if (allpoints != NULL) {
      const vtkIdType numPolyDataPoints = allpoints->GetNumberOfPoints();
      
      // Vector to map vtkIds to corresponding ids in WEM.
      std::vector<int> idsMap;
      idsMap.resize(numPolyDataPoints, -1);

      // Traverse all cells.
      const vtkIdType numCells = _poly->GetNumberOfCells();
      for (vtkIdType cellId = 0; cellId < numCells; ++cellId) {
        // Get cell and its points.
        vtkCell *cell = _poly->GetCell(cellId);
        if (cell != NULL){
          // Get all points of the cell.
          vtkPoints *points = cell->GetPoints();
          const vtkIdType numCellPoints = points ? points->GetNumberOfPoints() : 0;
          
          // If cell exists and if it has more than 2 points then add a face to the WEM.
          if (numCellPoints > 2){
            std::vector<int> faceids;
            // Parse all points of the vtkCell.
            for (vtkIdType pntId=0; pntId < numCellPoints; ++pntId) {
              double pnt[3]={0,0,0};
              points->GetPoint(pntId, pnt);
                
              // Get actual id of point in mesh and write it into position
              vtkIdType actualPointId = cell->GetPointId(pntId);
                
              // Add new node to Polyhedron or reuse an existing one if there
              // is an existing one corresponding to the the actualPointId.
              int node = nodeid;
              if (idsMap[actualPointId] != -1){
                // Node already exists, get it from node array.
                node = idsMap[actualPointId];
              } else {
                B.add_vertex( Point( static_cast<float>(pnt[0]),
                                     static_cast<float>(pnt[1]),
                                     static_cast<float>(pnt[2]) ) );
                    
 
                // OFF id corresponding to vtkIds.
                idsMap[actualPointId] = nodeid;
                ++nodeid;
              }
                 
              faceids.push_back( node );
            } // for(...)
            B.begin_facet();
            for ( unsigned int i=0; i<faceids.size(); ++i ) {
              B.add_vertex_to_facet( faceids[ i ] );
            }
            B.end_facet();
          }  // if (newFace)
        } // if (numCellPoints > 2)
      }  // if (cell)
    }
     
    B.end_surface();
  }
private:
  vtkPolyData * _poly;
};

template< class Kernel >
void VTKPolyDataToOFF( vtkPolyData * polyIn, CGAL::Polyhedron_3<Kernel> & polyOut ) {
  typedef typename CGAL::Polyhedron_3<Kernel> PolyType;

  Build_polygon<HalfedgeDS> poly( polyIn );
  polyOut.delegate( poly );
}
