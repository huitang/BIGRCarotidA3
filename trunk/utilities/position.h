#ifndef __voxel_H
#define __voxel_H

#include <cmath>

//-----------------------------------------------------------------
// Simple voxel type, used to store the voxel indices and vector
// elements
struct position {
  // Coordinates
  float _idx[4];

  // Vector
  float _idv[4];

  position() {
    _idx[0] = 0;
    _idx[1] = 0;
    _idx[2] = 0;
    _idx[3] = 0;
    _idv[0] = 0;
    _idv[1] = 0;
    _idv[2] = 0;
    _idv[3] = 0;
  }

  // Constructor
  position(float x, float y, float z, float r=0, float vx=0, float vy=0, float vz=0, float vr=0) {
    _idx[0] = x;
    _idx[1] = y;
    _idx[2] = z;
    _idx[3] = r;
    _idv[0] = vx;
    _idv[1] = vy;
    _idv[2] = vz;
    _idv[3] = vr;
  }

  // Operators
  bool operator== (const position &rhs) const {
    return 
      _idx[0]==rhs._idx[0] &&
      _idx[1]==rhs._idx[1] &&
      _idx[2]==rhs._idx[2] &&
      _idx[3]==rhs._idx[3] &&
      _idv[0]==rhs._idv[0] &&
      _idv[1]==rhs._idv[1] &&
      _idv[2]==rhs._idv[2] &&
      _idv[3]==rhs._idv[3];
  }
  
  bool operator< (const position &rhs) const {
    return ( _idx[0] < rhs._idx[0] ||
             ( _idx[0] == rhs._idx[0] &&
               ( _idx[1] < rhs._idx[1] ||
                 ( _idx[1] == rhs._idx[1] &&
                   ( _idx[2] < rhs._idx[2] ||
                     ( _idx[2] == rhs._idx[2] &&
                       ( _idx[3] < rhs._idx[3]) ) ) ) ) ) );
  }

  position operator- (const position &rhs) const {
    position temp = (*this);
    temp._idx[0] -= rhs._idx[0];
    temp._idx[1] -= rhs._idx[1];
    temp._idx[2] -= rhs._idx[2];
    temp._idx[3] -= rhs._idx[3];
    return temp;
  }

  position operator+ (const position &rhs) const {
    position temp = (*this);
    temp._idx[0] += rhs._idx[0];
    temp._idx[1] += rhs._idx[1];
    temp._idx[2] += rhs._idx[2];
    temp._idx[3] += rhs._idx[3];
    return temp;
  }

  position operator* (float f) const {
    position temp (*this);
    temp._idx[0] *= f;
    temp._idx[1] *= f;
    temp._idx[2] *= f;
    temp._idx[3] *= f;
    return temp;
  }

  position operator* (double f) const {
    position temp (*this);
    temp._idx[0] *= float(f);
    temp._idx[1] *= float(f);
    temp._idx[2] *= float(f);
    temp._idx[3] *= float(f);
    return temp;
  }

  position operator/ (float f) const {
    position temp (*this);
    temp._idx[0] /= f;
    temp._idx[1] /= f;
    temp._idx[2] /= f;
    temp._idx[3] /= f;
    return temp;
  }

  position operator/ (double f) const {
    position temp (*this);
    temp._idx[0] /= float(f);
    temp._idx[1] /= float(f);
    temp._idx[2] /= float(f);
    temp._idx[3] /= float(f);
    return temp;
  }

  // Dot product with r
  float dot (position r) const {
    return _idx[0]*r[0] + _idx[1]*r[1] + _idx[2]*r[2] + _idx[3]*r[3];
  }

  // Squared length
  float length2 () const {
    return _idx[0]*_idx[0]+_idx[1]*_idx[1]+_idx[2]*_idx[2]+_idx[3]*_idx[3];
  }

  // Length
  float length () const {
    return sqrt (length2());
  }

  // Vector length
  float vecLength() const {
    return sqrt (_idv[0]*_idv[0]+_idv[1]*_idv[1]+_idv[2]*_idv[2]+_idv[3]*_idv[3]);
  }

  // Normalize coordinates
  void normalize() {
    float l = length();
    _idx[0] /= l;
    _idx[1] /= l;
    _idx[2] /= l;
    _idx[3] /= l;
  }

  // Index operators
  float operator[](int i) const { return _idx[i]; }
  float &operator[](int i)  { return _idx[i]; }

  position operator* (float lhs) {
    position temp;
    temp[0] = _idx[0]*lhs;
    temp[1] = _idx[1]*lhs;
    temp[2] = _idx[2]*lhs;
    temp[3] = _idx[3]*lhs;
    return temp;
  }

  position operator* (double lhs) {
    position temp;
    temp[0] = _idx[0]*float(lhs);
    temp[1] = _idx[1]*float(lhs);
    temp[2] = _idx[2]*float(lhs);
    temp[3] = _idx[3]*float(lhs);
    return temp;
  }

  position operator/ (float lhs) {
    position temp;
    temp[0] = _idx[0]/lhs;
    temp[1] = _idx[1]/lhs;
    temp[2] = _idx[2]/lhs;
    temp[3] = _idx[3]/lhs;
    return temp;
  }

  position operator/ (double lhs) {
    position temp;
    temp[0] = _idx[0]/float(lhs);
    temp[1] = _idx[1]/float(lhs);
    temp[2] = _idx[2]/float(lhs);
    temp[3] = _idx[3]/float(lhs);
    return temp;
  }

};



//-----------------------------------------------------------------
// Simple voxel type, used to store the voxel indices and vector
// elements
struct positionShort {
  // Coordinates
  short _idx[4];

  positionShort() {
    _idx[0] = 0;
    _idx[1] = 0;
    _idx[2] = 0;
    _idx[3] = 0;
  }

  // Constructor
  positionShort(short x, short y, short z, short r=0) {
    _idx[0] = x;
    _idx[1] = y;
    _idx[2] = z;
    _idx[3] = r;
  }

  // Operators
  bool operator== (const positionShort &rhs) const {
    return 
      _idx[0]==rhs._idx[0] &&
      _idx[1]==rhs._idx[1] &&
      _idx[2]==rhs._idx[2] &&
      _idx[3]==rhs._idx[3];
  }

  bool operator< (const positionShort &rhs) const {
    return ( _idx[0] < rhs._idx[0] ||
      ( _idx[0] == rhs._idx[0] &&
      ( _idx[1] < rhs._idx[1] ||
      ( _idx[1] == rhs._idx[1] &&
      ( _idx[2] < rhs._idx[2] ||
      ( _idx[2] == rhs._idx[2] &&
      ( _idx[3] < rhs._idx[3]) ) ) ) ) ) );
  }

  positionShort operator- (const positionShort &rhs) const {
    positionShort temp = (*this);
    temp._idx[0] -= rhs._idx[0];
    temp._idx[1] -= rhs._idx[1];
    temp._idx[2] -= rhs._idx[2];
    temp._idx[3] -= rhs._idx[3];
    return temp;
  }

  positionShort operator+ (const positionShort &rhs) const {
    positionShort temp = (*this);
    temp._idx[0] += rhs._idx[0];
    temp._idx[1] += rhs._idx[1];
    temp._idx[2] += rhs._idx[2];
    temp._idx[3] += rhs._idx[3];
    return temp;
  }

  // Index operators
  short operator[](int i) const { return _idx[i]; }
  short &operator[](int i)  { return _idx[i]; }
};

#endif
