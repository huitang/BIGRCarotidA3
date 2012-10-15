#ifndef __splinefitter_H
#define __splinefitter_H

#include <vector>
#include "../utilities/position.h"

// Struct for linesegments
struct lineParam {
    float ax,bx,ay,by,az,bz;
    float r1,r2;
    float d1,d2;
};

typedef std::vector< lineParam > centerlineSpline;
typedef std::vector< position > centerlinePoints;
typedef std::pair< centerlinePoints, centerlineSpline > centerline;

float distanceToLine (const lineParam& lp, 
                      const position & pos,
                      const bool outSideSegMeas);

float distanceBetweenTwoLines (const lineParam& lp1, 
                               const lineParam& lp2);

float intersectionPlaneAndLine (const lineParam& plane, 
                                const lineParam& lp2,
                                position point);

void  addLines (std::vector<lineParam> & lineParams, 
                const float a[3], 
                const float b[3], 
                const float c[3], 
                const float d[3], 
                const float sampleDistance);

void sampleEquiDistant (const centerlineSpline & inputLines,
                        const float sampleDistance,
                        centerlineSpline & outputLines);

void splineToPath (const centerlineSpline & inputSpline,
                   const float sampleDistance,
                   const float markerType,
                   pathtype & outputMarkers,
                   const bool setVec=false);

void pathToSpline (const centerlinePoints & inputMarkers,
                   float sampleDistance,
                   centerlineSpline & outputSpline);

void smoothPath (const pathtype & inputpath, 
                 pathtype & outputpath,
                 const float sigma,
                 const float precision);

#include "splinefitter.cpp"

#endif
