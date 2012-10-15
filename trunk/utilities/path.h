#ifndef __path_H
#define __path_H

#include "position.h"
#include <vector>


#include <set>
#include <list>
#include <algorithm>
#include <iostream>

// Borrowed from Theo, order markers on length
//-----------------------------------------------------------------
// Given the input markers, output them reordered, such that the form
// a "linear path". E.g. if, of a list of three markers, the last
// markers is physically in between the other two, we want to search a
// path between the first and third marker, and between the third and
// second marker.  We want to find the ordering that minimizes the
// length of the lines between neighboring markers. Straightforward
// implementation by testing all permutations is not feasible, thus we
// implement a heuristic: the two closest markers probably are
// neighbors, and than we iteratively add the marker that is closest
// to either one of the markers at the front or end of the current
// list.

class Path {

public:
  // Define type for paths
  typedef std::vector< position > PathType;

  static PathType orderVoxels (const PathType &input);

  // Compute the sum of euclidian lengths between consecutive XMarkers
  static double computePathLength (const PathType & list);

  // Smooth path with gaussian kernel and resample
  static void smoothAndResamplePath( PathType    & path, 
                                     PathType    & outputpath, 
                                     const float sigma, 
                                     const float precision, 
                                     const float sampleDistance, 
                                     const bool  orderMarkers);

  static void smoothPath( const PathType & inputpath, 
                          PathType & outputpath,
                          const float sigma,
                          const float precision);
};

#endif
