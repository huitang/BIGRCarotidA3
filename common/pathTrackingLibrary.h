#ifndef __pathTracking_H
#define __pathTracking_H

#include <set>
#include "../common/path.h"

class PathTracking {
public:
  // Constructor.
  PathTracking (const std::vector<int> imgExt,
                   float* &data,
                   const std::vector<float> &voxelSizes);

  // Start path tracking using start(point) and end(point)
  void start(const position &start,
             const position &end);

  // Return encapsulated path
  pathtype * getPath();

  // Get reference to enclosed path
  pathtype & getPathPointer();

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
  pathtype traceBack(const int s,
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
  void findPathBetweenTwoMarkers (const position &begin, 
                                  const position &end, 
                                  std::vector<int> &summmedData);

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
  pathtype _path;

  long _num;
};

#endif // __mlPathTracking_H
