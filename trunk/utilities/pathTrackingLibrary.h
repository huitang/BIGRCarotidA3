#ifndef __pathTracking_H
#define __pathTracking_H

#include <set>
#include "path.h"
#include "itkImageBase.h"


class PathTracking {
public:
  // Constructor.
  PathTracking (const std::vector<int> imgExt,
                   float* &data,
                   const std::vector<float> &voxelSizes);

  typedef position PositionType;
  typedef Path::PathType PathType;

  // Start path tracking using start(point) and end(point)
  void start(const PositionType &start,
             const PositionType &end);

  // Return encapsulated path
  PathType * getPath();

  // Get reference to enclosed path
  PathType & getPathPointer();

  void findPathBetweenTwoMarkers (const PositionType &begin, 
                                  const PositionType &end, 
                                  std::vector<int> &summmedData);

  static void writePathToFile ( const std::string & filename, 
                                PathType          & path, 
                                bool                outputVec, 
                                bool                world, 
                                float               x, 
                                float               y, 
                                float               z );

  static void writePathToFile ( const std::string & filename, 
                                PathType          & path, 
                                itk::ImageBase<3>::Pointer image );

  static void writePathToFile4D ( const std::string & filename, 
                                  const PathType    & path, 
                                  const bool          world, 
                                  const float         x, 
                                  const float         y, 
                                  const float         z, 
                                  const float         r );

  static bool loadPathFromFile ( const std::string & inputfile, 
                                 PathType          & path  );

private:
  // Minimum cost returned by frangi3D (* sigmoid)
  static const float MINCOSTS;

  // Data pointer
  float* _data;

  // Function to compute distances between center point and its 26 neighbours
  // for the voxel sizes of the input volumes
  void fillNeighborDistances();

  // Trace back the path from e to s using referenceData image
  // In this image, the value of a voxel references to its parent voxel
  PathType traceBack(const int s,
                     const int e,
                     std::vector<int> &referenceData);

  // Process all 26 candidates of the currently visited voxel
  void process26Candidates(
                const int candidate,
                std::set< std::pair<float, int> > &candidates,
                int &c1,
                int &c2,
                bool &finished, 
                std::vector<int> &referenceData,
                double currentSum);

  // Process one of the 26 neighbors
  void processNeighbor( 
                const int candidate,
                const int nb,
                const int candsum,
                const double *d,
                bool &finished,
                int &c1,
                int &c2,
                std::set< std::pair<float, int> > &candidates,
                std::vector<int> &summmedData,
                double currentSum);

  // Find minimal cost path between begin and end using summeddata to store
  // references in

  // Array to store neighbour distances in
  double _neighborDist[27];

  // Strides input images
  const int _imageStrideY;
  const int _imageStrideZ;

  // Max index
  const int _maxIndex;

  // Increment when moving in z and y from last x(, y) voxel
  const int _zIncrement;
  const int _yIncrement;

  // Size of input images (assumed to be equal for all input images)
  const std::vector<int> _imgExt;
  
  // Voxel sizes of input images (assumed to be equal for all input images)
  const std::vector<float> & _voxelSizes;
  
  // The resulting path
  PathType _path;

  long _num;
};

#endif // __mlPathTracking_H
