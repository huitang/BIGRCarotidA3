#ifndef _objtocgalpolyhedron_h 
#define _objtocgalpolyhedron_h 

#include <CGAL/Polyhedron_3.h>

typedef Polyhedron::HalfedgeDS HalfedgeDS;

// A modifier creating a Polyhedron from a OBJ file.
template <class HDS>
class Build_polygon_obj : public CGAL::Modifier_base<HDS> {
public:
  Build_polygon_obj( const std::string & filename ) {
    _filename = filename;
  }

  void operator()( HDS& hds) {
    std::ifstream file;
    file.open( _filename.c_str(), ios::in );

    std::vector< std::vector< float > > points;
    std::vector< std::vector< int > > faces;
    
    // Read points and vertices from OBJ file
    while ( true )
    {
      std::string type, c0, c1, c2;
      if ( !( file >> type ) ) break;
      if ( !( file >> c0 ) ) break;
      if ( !( file >> c1 ) ) break;
      if ( !( file >> c2 ) ) break;
      if ( type == "v" )
      {
        std::vector< float > point;
        point.push_back( atof( c0.c_str() ) );
        point.push_back( atof( c1.c_str() ) );
        point.push_back( atof( c2.c_str() ) );
        points.push_back( point );
      }
      else if ( type == "f" )
      {
        std::vector< int > face;
        face.push_back( atoi( c0.c_str() ) );
        face.push_back( atoi( c1.c_str() ) );
        face.push_back( atoi( c2.c_str() ) );
        faces.push_back( face );
      }
    }

    file.close();

    // Start polygon
    CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
    B.begin_surface( points.size(), faces.size() );
    typedef typename HDS::Vertex Vertex;
    typedef typename Vertex::Point Point;

    for ( size_t p = 0; p < points.size(); ++p )
    {
      B.add_vertex( Point( points[ p ][ 0 ],
                           points[ p ][ 1 ],
                           points[ p ][ 2 ] ) );
    }

    for ( size_t f = 0; f < faces.size(); ++f )
    {
      B.begin_facet();
      for ( unsigned int i = 0; i < 3; ++i ) {
        B.add_vertex_to_facet( faces[ f ][ i ] - 1 );
      }
      B.end_facet();
    }
    
    B.end_surface();
  }
private:
  std::string _filename;
};

template< class Kernel >
void ObjToOFF( const std::string & filename, CGAL::Polyhedron_3< Kernel > & polyOut ) 
{
  typedef typename CGAL::Polyhedron_3< Kernel > PolyType;

  Build_polygon_obj< HalfedgeDS > poly( filename );
  polyOut.delegate( poly );
}

#endif

