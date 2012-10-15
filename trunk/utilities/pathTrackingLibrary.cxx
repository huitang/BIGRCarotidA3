// Local includes
#include "pathTrackingLibrary.h"

#include <cmath>
#include <iostream>
#include <time.h>
#include <algorithm>
#include <stdlib.h>
#include <assert.h>

#include "utils.h"

#undef throw

const float PathTracking::MINCOSTS = 0.001f;

//------------------------------------------------------------------------------
//! Constructor
//------------------------------------------------------------------------------
PathTracking::PathTracking (const std::vector<int> imgExt,
                               float* &data,
                               const std::vector<float> &voxelSizes):
  _data (data),
  _imageStrideY (imgExt[0]),
  _imageStrideZ (imgExt[0]*imgExt[1]),
  _maxIndex (imgExt[0]*imgExt[1]*imgExt[2]-1),
  _zIncrement (_imageStrideZ - 3*_imageStrideY),
  _yIncrement (_imageStrideY - 3),
  _imgExt (imgExt),
  _voxelSizes (voxelSizes),
  _num (0) {
}

// Start path tracking with startpoint 'start' and endpoint 'end'
void PathTracking::start (const PositionType & start, const PositionType & end) 
{
  // referenceData contains ints that point to the neighbour a voxel is 
  // visited from, for a fast traceback
  std::vector<int> referenceData (_imgExt[2]*_imgExt[1]*_imgExt[0], 0);

  // Compute voxel distances between neighbors
  fillNeighborDistances();

  // Find minimal cost path between markers
  findPathBetweenTwoMarkers(start, end, referenceData);
}

// Precompute distances between voxels
void PathTracking::fillNeighborDistances() 
{
  // Voxel dimensions, initialized to unity.
  const float &vx = _voxelSizes[0];
  const float &vy = _voxelSizes[1];
  const float &vz = _voxelSizes[2];
  
  // Traverse all neighbors
  double *d = _neighborDist;
  
  // 26-connected: determine euclidean distance.
  for (int z = -1; z <= 1; ++z)
    for (int y = -1; y <= 1; ++y)
      for (int x = -1; x <= 1; ++x) {
        *d++ = sqrt ( x*x*vx*vx +
                      y*y*vy*vy +
                      z*z*vz*vz );
      }
}

// Given a reference-image, trace back the optimal (mincost) path from
// starting point s to end point e.
PathTracking::PathType PathTracking::traceBack(const int s, const int e, 
                                            std::vector<int> &referenceData) {
  // Variable holding path traced
  PathType result;
  
  // Start at index s
  int index = s;
  // While the next position is not equal to the end position
  while (index != e) {
    // Add voxel to the path
    int indexModStrides2 = index%_imageStrideZ;
    result.push_back(position(float(indexModStrides2%_imageStrideY),
                              float(indexModStrides2/_imageStrideY),
                              float(index/_imageStrideZ)));
    // Get next voxel
    index = abs (referenceData[index]);
  }

  // Add last voxel
  int eModStrides2 = e%_imageStrideZ;
  result.push_back(position(float(eModStrides2%_imageStrideY),
                            float(eModStrides2/_imageStrideY),
                            float(e/_imageStrideZ)));

  // Return path
  return result;
}

// Find the minimal cost path between begin and end. Store parent references.
// If no path is found, an empty path is returned.
void PathTracking::findPathBetweenTwoMarkers (const PositionType &begin, 
                                              const PositionType &end, 
                                              std::vector<int> &referenceData) {
  // Start clock for timing
  clock_t start = clock();

  // Clear the output path
  _path.clear();

  // Check bounds for seedpoints
  if (begin[0] < 0 || begin[1] < 0 || begin[2] < 0 ||
      end[0] < 0 || end[1] < 0 || end[2] < 0 ||
      begin[0] >= _imgExt[0] || begin[1] >= _imgExt[1] || 
      begin[2] >= _imgExt[2] || end[0] >= _imgExt[0] || 
      end[1] >= _imgExt[1] || end[2] >= _imgExt[2]) {
    std::cerr << "Begin or end voxel is outside valid region" << std::endl;
    return;
  }

  // Compute the array index of the start and end position
  int v1ImgIndex = _imageStrideZ * int (begin[2]) + _imageStrideY * int (begin[1]) + int (begin[0]);
  int v2ImgIndex = _imageStrideZ * int (end[2]) + _imageStrideY * int (end[1]) + int (end[0]);

  // Create candidates queues
  std::set< std::pair<float, int> > candidates;

  // Compute start and endpoint costs
  float v1cost, v2cost;
  v1cost = _data[v1ImgIndex];
  v2cost = _data[v2ImgIndex];

  // Initialize candidates and sum image with begin and end marker
  candidates.insert(std::make_pair(v1cost,v1ImgIndex));
  candidates.insert(std::make_pair(v2cost,v2ImgIndex));
  referenceData[v1ImgIndex] = 1;
  referenceData[v2ImgIndex] = -1;

  bool pathFound = false;

  // Voxels for the connection. Path is from s to c1, and from c2
  // to end, where c1 and c2 are neighbors.
  int c1=0, c2=0;

  while (!candidates.empty() && !pathFound) {
    // Get first candidate
    std::pair<float,int> cand = *candidates.begin();
    // Remove first element
    candidates.erase(candidates.begin());
    _num++;

    // Get cost and voxel
    const float candcost  = cand.first;
    const int   candvoxel = cand.second;

    // Process neighbours
    process26Candidates(candvoxel, candidates, c1, c2, pathFound,
                                   referenceData, double (candcost));
  }

  if (pathFound) {
    // Successfully found a connection, now use the sum image to
    // trace back the path from the connection to the begin and
    // end point.
    PathType p1 = traceBack(c1,v1ImgIndex,referenceData);
    PathType p2 = traceBack(c2,v2ImgIndex,referenceData);

    if (!p1.empty() || !p2.empty()) {
      _path.insert(_path.end(), p1.rbegin(), p1.rend());
      _path.insert(_path.end(), p2.begin(), p2.end());
    }
  }

  // Stop timer and output running time
  clock_t ready = clock();
  std::cout << "Running time of path tracking: " << float(ready-start)/CLOCKS_PER_SEC << std::endl;
  std::cout << "Num voxels processed: " << _num << std::endl;
}

//----------------------------------------------------------------------------
// Processes all 26-connected (8-connected if 2D) neighbors of a
// candidate. 
//----------------------------------------------------------------------------
void PathTracking::process26Candidates(const int candidate,
                   std::set< std::pair<float, int> > &candidates,
                   int &c1, int &c2, bool &finished, 
                   std::vector<int> & referenceData, double currentSum) {
  // Compute index value of top left front neighbour
  int indexNeighbour = candidate - _imageStrideZ - _imageStrideY - 1;

  // Get reference to node where we came from before visiting this one.
  // Sign of reference indicates the front it belongs to
  int candsum = referenceData[candidate];

  double *d = _neighborDist;
  // Process all neigbours
  for (int z = -1; z < 2; ++z, indexNeighbour+=_zIncrement) {
    for (int y = -1; y < 2; ++y, indexNeighbour+=_yIncrement) {
      for (int x = -1; x < 2; ++x, ++d, ++indexNeighbour) {
        if (!x && !y && !z) continue; // Current point, skip it

        if (indexNeighbour<0) continue;
        if (indexNeighbour>_maxIndex) continue;

        processNeighbor(candidate, indexNeighbour, candsum, d,
                        finished, c1, c2, candidates, referenceData, currentSum);

        if (finished) return;
      }
    }
  }
}

void PathTracking::processNeighbor( const int candidate,
                     const int nb,
                     const int candsum,
                     const double *d,
                     bool &finished,
                     int &c1, int &c2,
                     std::set< std::pair<float, int> > &candidates,
                     std::vector<int> &referenceData,
                     double currentSum) {
  // Get reference to parent of this neighbor, should be zero if not already
  // visited. Sign of value determines front voxel belongs to
  int referenceNb = referenceData[nb];
  if (((referenceNb<0 && candsum>0) || (referenceNb>0 && candsum<0)) && !finished) {
    // If signs of references are not equal, the connection has been found.
    c1 = candsum > 0 ? candidate:nb;
    c2 = candsum > 0 ? nb:candidate;
    finished = true;
    return;
  }
  // If the referenceNb = 0, the voxel has not been visited before
  if (referenceNb == 0) {
    // Get cost for neighbor from cost image
    double nbcost = currentSum+_data[nb];

    // Determine front, neighbor belongs to and create reference to parent
    const int nbindex = candsum > 0 ? candidate : -candidate;
    referenceData[nb] = nbindex;

    candidates.insert(std::make_pair(float(nbcost),nb));
  }
}

// Return pointer to encapsulated path
PathTracking::PathType * PathTracking::getPath () {
  return &_path;
}

// Return reference to encapsulated path
PathTracking::PathType & PathTracking::getPathPointer() {
  return _path;
}


// Write a path to a file
// By setting world to true, coordinates are transformed to world coordinates (+0.5)*voxelsize
void PathTracking::writePathToFile (const std::string & filename, PathType & path, bool outputVec=false, bool world=false, float x=1.0f, float y=1.0f, float z=1.0f) {
  // Open file
  std::fstream file_op(filename.c_str(), std::ios::out);

  // Check if file is open
  if (file_op.is_open()) {
    // Write path to file
    for (PathType::iterator it = path.begin(); it != path.end(); ++it) {
      position v = *it;
      float vx = v[0];
      float vy = v[1];
      float vz = v[2];
      if (world) {
        vx = (vx + 0.5f) * x;
        vy = (vy + 0.5f) * y;
        vz = (vz + 0.5f) * z;
      }
      file_op << vx << " " << vy << " " << vz;
      if (outputVec) {
        file_op << " " << v._idv[0] << " " << v._idv[1] << " " << v._idv[2];
      }
      file_op << std::endl;
    }
    file_op.close();
    std::cout << "Path written to " << filename.c_str() << std::endl;
  } else {
    std::cout << "Cannot write path to " << filename.c_str() << std::endl;
  }  
}


// Write a path to a file
// By setting world to true, coordinates are transformed to world coordinates (+0.5)*voxelsize
void PathTracking::writePathToFile (const std::string & filename, PathType & path, itk::ImageBase<3>::Pointer image) {
  // Open file
  std::fstream file_op(filename.c_str(), std::ios::out);

  // Check if file is open
  if (file_op.is_open()) {
    // Write path to file
    for (PathType::iterator it = path.begin(); it != path.end(); ++it) {
      position v = *it;
      itk::ImageBase<3>::IndexType idx;
      idx[0] = v[0];
      idx[1] = v[1];
      idx[2] = v[2];
      itk::ImageBase<3>::PointType vw;
      image->TransformIndexToPhysicalPoint(idx,vw);
      file_op << vw[0] << " " << vw[1] << " " << vw[2];
      file_op << std::endl;
    }
    file_op.close();
    std::cout << "Path written to " << filename.c_str() << std::endl;
  } else {
    std::cout << "Cannot write path to " << filename.c_str() << std::endl;
  }  
}


// Write a path to a file
// By setting world to true, coordinates are transformed to world coordinates (+0.5)*voxelsize
void PathTracking::writePathToFile4D (const std::string & filename, const PathType & path, const bool world=false, const float x=1.0f, const float y=1.0f, const float z=1.0f, const float r=1.0f) {
  // Open file
  std::fstream file_op(filename.c_str(), std::ios::out);

  // Check if file is open
  if (file_op.is_open()) {
    // Write path to file
    for (PathType::const_iterator it = path.begin(); it != path.end(); ++it) {
      position v = *it;
      float vx = v[0];
      float vy = v[1];
      float vz = v[2];
      float vr = v[3];
      if (world) {
        vx = (vx + 0.5f) * x;
        vy = (vy + 0.5f) * y;
        vz = (vz + 0.5f) * z;
        //vr = vr; //(vr + 0.5f) * r;
      }
      file_op << vx << " " << vy << " " << vz << " " << vr << std::endl;
    }
    file_op.close();
    std::cout << "Path written to " << filename.c_str() << std::endl;
  } else {
    std::cout << "Cannot write path to " << filename.c_str() << std::endl;
  }  
}


// Load a path from file
bool PathTracking::loadPathFromFile (const std::string & inputfile, PathType & path) {
  path.clear();

  // Check filename
  if (inputfile != "") {
    // Open file
    std::fstream file_op(inputfile.c_str(), std::ios::in);
    // Check if file is open
    if (file_op.is_open()) {
      std::string line;
      while (std::getline(file_op, line)) {
        std::vector<std::string> tokens;
        tokenize(line, tokens);
        position v (float(atof (tokens[0].c_str())), float(atof(tokens[1].c_str())), float(atof(tokens[2].c_str())));
        if (tokens.size()>5) {
          v._idv[0] = float(atof(tokens[3].c_str()));
          v._idv[1] = float(atof(tokens[4].c_str()));
          v._idv[2] = float(atof(tokens[5].c_str()));
        }
        if (tokens.size()==4) {
          v._idv[0] = float(atof(tokens[3].c_str()));
          v._idv[1] = 0;
          v._idv[2] = 0;
        }
        path.push_back (v);
      }
      file_op.close();
    } else {
      return false;
    }
  } else {
    return false;
  }

  return true;
}
