#ifndef __splinefitter_cpp
#define __splinefitter_cpp

#include <iostream>
#include "splinefitter.h"

float SplineFitter::distanceToLine(const lineParam& lp, const position & pos, const bool outSideSegMeas = true) {
  const float x = pos[0];
  const float y = pos[1];
  const float z = pos[2];

  float t = (lp.ax * (x - lp.bx) + lp.ay*(y - lp.by) + lp.az*(z - lp.bz)) /
            (lp.ax * lp.ax + lp.ay * lp.ay + lp.az * lp.az);

  if(t < 0) {
    if(!outSideSegMeas)
      return 10E5;
    t = 0;
  }

  if(t > 1) {
    if(!outSideSegMeas)
      return 10E5;
    t = 1;
  }

  const float lx = lp.bx + t * lp.ax;
  const float ly = lp.by + t * lp.ay;
  const float lz = lp.bz + t * lp.az;

  const float dist =  sqrt(  (lx - x)*(lx - x) + 
                            (ly - y)*(ly - y) +
                            (lz - z)*(lz - z));

  return dist;
}

float SplineFitter::intersectionPlaneAndLine (const lineParam& plane, const lineParam& lp2) {
  const position planePos(plane.bx, plane.by, plane.bz);
  position N (plane.ax, plane.ay, plane.az);
  N.normalize();

  const position P1 (lp2.bx, lp2.by, lp2.bz);
  const position lineDir (lp2.ax, lp2.ay, lp2.az);

  const position P2 = P1 + lineDir;

  const float divValue  = N.dot(P2 - P1);
  if(divValue == 0) {
    return distanceToLine(lp2, planePos, false);
  } else {
    float mu = N.dot(planePos - P1) / divValue;
    if(mu < 0)
      return 10E5;

    if(mu > 1)
      return 10E5;
    position muLineDir (mu*lineDir[0], mu*lineDir[1], mu*lineDir[2]);
    const position posOnLine2 = P1 + muLineDir;

    return (planePos - posOnLine2).length();
  }
}

float SplineFitter::distanceBetweenTwoLines(const lineParam & lp1, const lineParam & lp2) {
  position r1(lp1.bx, lp1.by, lp1.bz);
  position r2(lp2.bx, lp2.by, lp2.bz);
  position a1(lp1.ax, lp1.ay, lp1.az);
  a1.normalize();
  position e1 = a1;
  position a2(lp2.ax, lp2.ay, lp2.az);
  a2.normalize();
  position e2 = a2;
  position r12 = r2 - r1;

  float divValue = (1 - (e1.dot(e2)) * (e1.dot(e2)));

  if(divValue != 0) {
    float lambda0 =  (r12.dot(e1) - r12.dot(e2)*(e1.dot(e2)) ) / divValue;
    float tOnLine1 = lambda0 / a1.length();
    if(tOnLine1 < 0)
      tOnLine1 = 0;
    if(tOnLine1 > 1)
      tOnLine1 = 1;
    position tOnLine1A1 (tOnLine1 * a1[0], tOnLine1 * a1[1], tOnLine1 * a1[2]);
    position dPos = r1 + tOnLine1A1;

    position pos(dPos[0],dPos[1],dPos[2]);

    return distanceToLine(lp2, pos, false);

  } else { //Lines parallel
    position p11 = r1;
    position p12 = r1 + a1;

    position p21 = r2;
    position p22 = r2 + a2;

    float dist1 = (p11 - p21).length();
    float dist2 = (p11 - p22).length();
    float dist3 = (p12 - p21).length();
    float dist4 = (p12 - p22).length();

    float minDist = dist1;
    if(dist2 < minDist)
      minDist = dist2;
    if(dist3 < minDist)
      minDist = dist3;
    if(dist4 < minDist)
      minDist = dist4;

    return minDist;
  }
}

void SplineFitter::addLines(std::vector<lineParam>& lineParams, const float a[3], const float b[3], const float c[3], const float d[3], const float sampleDistance) {
  float distStartEnd = 0;
  for(int i = 0; i < 3; i++) {
    const float compDif = a[i] + b[i] + c[i];
    distStartEnd += compDif*compDif;
  }
  distStartEnd  = sqrt(distStartEnd);
  int nrLines = int (floor (float (distStartEnd) / sampleDistance)) + 1;

  for(int p = 0; p < nrLines; p ++) {
    lineParam lp;
    const float t1 = (float)p / nrLines;
    const float t2 = (float)(p+1) / nrLines;

    for(int component = 0; component < 3; component++) {
      const float cx1 = a[component] * t1 * t1 * t1 + b[component] * t1 * t1 + c[component] * t1 + d[component];
      const float cx2 = a[component] * t2 * t2 * t2 + b[component] * t2 * t2 + c[component] * t2 + d[component];

      switch(component) {
        case 0:
          lp.bx = cx1;
          lp.ax = cx2-cx1;
          break;
        case 1:
          lp.by = cx1;
          lp.ay = cx2-cx1;
          break;
        case 2:
          lp.bz = cx1;
          lp.az = cx2-cx1;
          break;
      }
    }
    lineParams.push_back(lp);
  }
}

void SplineFitter::sampleEquiDistant(const centerlineSpline& inputLines,
                       const float sampleDistance,
                       centerlineSpline& outputLines) {
  if(!inputLines.empty()) {
    centerlineSpline copyInputLines(inputLines);
    std::copy(inputLines.begin(), inputLines.end(), copyInputLines.begin());

    //Last point
    lineParam lastLp;
    lastLp.bx = (inputLines.end() - 1)->bx + (inputLines.end() - 1)->ax;
    lastLp.by = (inputLines.end() - 1)->by + (inputLines.end() - 1)->ay;
    lastLp.bz = (inputLines.end() - 1)->bz + (inputLines.end() - 1)->az;

    copyInputLines.push_back(lastLp);
    outputLines.clear();
    
    //Calculate distance from each point to each spline in the reference 
    centerlineSpline::const_iterator cpIt = copyInputLines.begin();

    position currentPoint(cpIt->bx, cpIt->by, cpIt->bz);
    position prevPoint = currentPoint;
    position prevDir;

    for(; cpIt != copyInputLines.end(); ) {
      float distToPrev = 0;
      //Calculate position of next point
      position prevPointInSeg = prevPoint;
      while((distToPrev < sampleDistance) && (cpIt != copyInputLines.end())) {
        currentPoint = position (cpIt->bx, cpIt->by, cpIt->bz);
        distToPrev  += (currentPoint - prevPointInSeg).length();

        if(distToPrev > sampleDistance) {
          const float tooFar = (distToPrev - sampleDistance);
          currentPoint[0] = currentPoint[0] - prevDir[0] * tooFar / prevDir.length();
          currentPoint[1] = currentPoint[1] - prevDir[1] * tooFar / prevDir.length();
          currentPoint[2] = currentPoint[2] - prevDir[2] * tooFar / prevDir.length();
          distToPrev += sampleDistance;
        } else {
          prevDir = position(cpIt->ax, cpIt->ay, cpIt->az);
          ++cpIt;
        }
        prevPointInSeg = currentPoint;
      }
      lineParam lp;
      lp.ax = (currentPoint[0] - prevPoint[0]);
      lp.ay = (currentPoint[1] - prevPoint[1]);
      lp.az = (currentPoint[2] - prevPoint[2]);
      lp.bx = prevPoint[0];
      lp.by = prevPoint[1];
      lp.bz = prevPoint[2];
      outputLines.push_back(lp);
      prevPoint = currentPoint;
    }
  }
}

void SplineFitter::splineToPath(  const centerlineSpline& inputSpline,
                    const float sampleDistance,
                    PathType & outputpath) {
  outputpath.clear();
  if(!inputSpline.empty()) {
    //Output voxels
    std::vector<lineParam>::const_iterator splineIt = inputSpline.begin();
    for(;splineIt != inputSpline.end(); ++splineIt) {
      position v (splineIt->bx, splineIt->by, splineIt->bz);
      outputpath.push_back(v);
    }
    --splineIt;
    position vOut;
    vOut[0] = splineIt->bx + splineIt->ax;
    vOut[1] = splineIt->by + splineIt->ay;
    vOut[2] = splineIt->bz + splineIt->az;
    outputpath.push_back(vOut);
  }
}

void SplineFitter::pathToSpline(  const centerlinePoints & inputpath,
                    float sampleDistance,
                    centerlineSpline& outputSpline) {
  outputSpline.clear();
  for(int i = 0; i < int (inputpath.size()) - 1; i++) {
    position pos0 = inputpath[i];
    const position pos1 = inputpath[i];
    const position pos2 = inputpath[i+1];
    position pos3 = inputpath[i+1];

    if(i > 0)
      pos0 = inputpath[i-1];

    if(i < int (inputpath.size()) - 2)
      pos3 = inputpath[i+2];

    float a[3], b[3], c[3], d[3];
    for(int component=0; component < 3; component++) {
      const float x0 = pos0[component];
      const float x1 = pos1[component];
      const float x2 = pos2[component];
      const float x3 = pos3[component];

      c[component] = (x2 - x0) / 2.0f;
      d[component] = x1;
      a[component] = ((x3- x1)/2.0f + c[component] - 2.0f * (x2 - x1));
      //a[component] = ((x3- x1)/2 - c[component] - 2 * (x2 - x1) + 2*c[component]) / 5;
      b[component] = (x2 - x1) - c[component] - a[component];
    }

    addLines(outputSpline, a, b, c, d, sampleDistance);
  }
}

#endif
