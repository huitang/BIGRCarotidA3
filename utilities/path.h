#ifndef __path_H
#define __path_H

#include "../utilities/position.h"
#include <vector>

// Define type for paths
typedef std::vector< position > pathtype;

#include "../utilities/splinefitter.h"
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
pathtype orderVoxels (const pathtype &input) {
  if (input.size() < 3) return input;

  // Find two closest markers
  double minlen = (input[0]-input[1]).length2();
  std::pair<int, int> minpair(0,1);
  for (int i = 0; i < int(input.size()); ++i)
    for (int j = i+1; j < int(input.size()); ++j) {
      double len = (input[i]-input[j]).length2();
      if (len < minlen) {
        minpair = std::make_pair<int,int>(i,j);
        minlen = len;
      }
    }

  // Remaining markers to process: all except the pair of markers that
  // are closest to each other.
  std::set<int> remaining;
  for (int x = 0; x < (int)input.size(); ++x)
    if (x != minpair.first && x != minpair.second) remaining.insert(x);

  // Series contains the indices of the markers in the correct order.
  std::list<int> series;
  series.push_back(minpair.first);
  series.push_back(minpair.second);
  int curfront = minpair.first;
  int curback = minpair.second;

  // While not all markers processed.
  while (!remaining.empty()) {
    // Find closest to front and back of current list
    double frontdistmin = 1E30;
    double backdistmin = 1E30;
    int idfrontmin = -1;
    int idbackmin = -1;
    // Test all remaining markers ...
    for (std::set<int>::iterator s = remaining.begin(); s != remaining.end(); ++s) {
      // If closer to front ...
      double frontdist = (input[curfront]-input[*s]).length2();
      if (frontdist < frontdistmin) {
        frontdistmin = frontdist;
        idfrontmin = *s;
      }
      // If closer to back...
      double backdist = (input[curback]-input[*s]).length2();
      if (backdist < backdistmin) {
        backdistmin = backdist;
        idbackmin = *s;
      }
    }
    // Pick marker that is closest to either the front or the back
    if (frontdistmin < backdistmin) {
      // marker to front is closest
      series.push_front(idfrontmin);
      curfront = idfrontmin;
      remaining.erase(idfrontmin);
    } else {
      // marker to back is closest
      series.push_back(idbackmin);
      curback = idbackmin;
      remaining.erase(idbackmin);
    }
  }

  // We prefer to have the first marker at the front.
  if (series.back() == 0)
    std::reverse(series.begin(), series.end());

  // Fill output marker list.
  pathtype result;
  for (std::list<int>::iterator j = series.begin(); j != series.end(); ++j) {
    result.push_back(input[*j]);
  }

  return result;
}

// Compute the sum of euclidian lengths between consecutive XMarkers
double computePathLength (const pathtype & list) {
  double length = 0.0;
  for (pathtype::const_iterator it=list.begin(); it<list.end()-1; ++it) {
    position current = *it;
    position next = *(it+1);
    // Add length of line from current to next marker
    length += (next-current).length();
  }
  return length;
}

// Smooth path with gaussian kernel and resample
void smoothAndResamplePath (pathtype & path, pathtype & outputpath, const float sigma=1.0f, const float precision=5.0f, const float sampleDistance=1.0f, const bool orderMarkers=true) {
  const int nrVoxels = path.size();

  outputpath.clear();
  if (nrVoxels > 3) {
    // Sort voxels according to euclidean distance
    pathtype sortedMarkers;
    if (orderMarkers) {
      sortedMarkers = orderVoxels (path);
    } else {
      std::copy (path.begin(), path.end(), back_inserter(sortedMarkers));
    }

    // Smooth markers with Gaussian kernel
    pathtype smoothedMarkers;
    smoothPath (sortedMarkers, smoothedMarkers, sigma, precision);

    // Fit 3th order spline through sorted markers and return nrSplinePoints 
    // line segments per spline segment (a spline segment is a section between 
    // input markers)
    std::vector<lineParam> centerSplineLines;
    pathToSpline (smoothedMarkers, sampleDistance, centerSplineLines);

    // Sample with equidistance along the line segments on the spline
    centerlineSpline equiDistSpline;
    sampleEquiDistant (centerSplineLines, sampleDistance, equiDistSpline);

    // Output a path at the start of the equidistance sampled line segments
    splineToPath (equiDistSpline, sampleDistance, outputpath);
  }
}

void smoothPath( const pathtype & inputpath, 
                 pathtype & outputpath,
                 const float sigma,
                 const float precision) {
  const float sqrt2pi = float (sqrt(2 * 3.14159265358));

  // Clear output path
  outputpath.clear();

  int nrVoxels = inputpath.size();
  if (nrVoxels > 3) {
    // Compute cumulative distances across path, from startpoint
    // Cumulative distance from endpoint = cumulative distance last
    // point from startpoint - cumulative distance current point
    // from startpoint
    std::vector<float> cumulativeDistance (inputpath.size(), 0.0f);
    for (int i=1; i<int(inputpath.size()); ++i) {
      position pos = inputpath[i];
      position lastPos = inputpath[i-1];
      cumulativeDistance[i] = (pos-lastPos).length() + cumulativeDistance[i-1];
    }

    // Distance threshold depends on sigma kernel and precision
    float distThreshold = sigma * precision;

    //Put markers in list
    for(int i=0; i<nrVoxels; ++i) {
      // Create output marker
      position voxelOut;
      if (sigma < 0.0001) {
        // Sigma too small, output unmodified position
        voxelOut = position (inputpath[i][0], inputpath[i][1], inputpath[i][2]);
      } else {
        // Check distanceThreshold
        float curDistThreshold = distThreshold;
        if (curDistThreshold > cumulativeDistance[i]) {
          // Distance to start marker is smaller than kernel size
          curDistThreshold = cumulativeDistance[i];
        }
        float cumDistToEnd = *(cumulativeDistance.end()-1)-cumulativeDistance[i];
        if (curDistThreshold > cumDistToEnd) {
          // Distance to end marker is smaller than kernel size
          curDistThreshold = cumDistToEnd;
        }

        // Now we know the real kernelsize (curDistThreshold)
        // Walk over positions

        // Maintain sum of values (for every dimension) and sum of weights
        std::vector<float> sumValues (3, 0.0f);
        float sumWeights = 0.0f;

        // First walk to the left
        float cumDist  = 0.0f;
        int   curIndex = i;
        while (cumDist <= curDistThreshold) {
          // Compute gaussian weight
          float weight = 1.0f/(sqrt2pi * sigma) * exp(-(cumDist*cumDist)/(2.0f*sigma*sigma));

          // Add values
          sumValues[0] += weight*inputpath[curIndex][0];
          sumValues[1] += weight*inputpath[curIndex][1];
          sumValues[2] += weight*inputpath[curIndex][2];
          sumWeights += weight;

          // Compute distance to left marker
          cumDist += (inputpath[curIndex]-inputpath[curIndex-1]).length();
          curIndex--;
        }

        if (i+1 < nrVoxels) {
          curIndex = i+1;
          position diff = inputpath[i]-inputpath[curIndex];
          cumDist = sqrt (pow (diff[0], 2.0f) + pow (diff[1], 2.0f) + pow (diff[2], 2.0f));
          while (cumDist <= curDistThreshold) {
            // Compute gaussian weight
            float weight = 1.0f/(sqrt2pi * sigma) * exp(-(cumDist*cumDist)/(2.0f*sigma*sigma));

            // Add values
            sumValues[0] += weight*inputpath[curIndex][0];
            sumValues[1] += weight*inputpath[curIndex][1];
            sumValues[2] += weight*inputpath[curIndex][2];
            sumWeights += weight;

            // Compute distance to right marker
            cumDist += (inputpath[curIndex]-inputpath[curIndex+1]).length();
            curIndex++;
          }
        }

        voxelOut[0] = sumValues[0] / sumWeights;
        voxelOut[1] = sumValues[1] / sumWeights;
        voxelOut[2] = sumValues[2] / sumWeights;
      }
      // Add smoothed voxel to output
      outputpath.push_back(voxelOut);
    }
  }
}

#endif
