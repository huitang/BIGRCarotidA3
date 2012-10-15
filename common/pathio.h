#ifndef __pathio_H
#define __pathio_H

// Standard includes
#include <string>
#include <fstream>
#include <iostream>

// Utils include (for tokenize)
#include "../common/utils.h"

// Write a path to a file
// By setting world to true, coordinates are transformed to world coordinates (+0.5)*voxelsize
void writePathToFile (const std::string & filename, pathtype & path, bool outputVec=false, bool world=false, float x=1.0f, float y=1.0f, float z=1.0f) {
  // Open file
  std::fstream file_op(filename.c_str(), std::ios::out);

  // Check if file is open
  if (file_op.is_open()) {
    // Write path to file
    for (pathtype::iterator it = path.begin(); it != path.end(); ++it) {
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
void writePathToFile4D (const std::string & filename, const pathtype & path, const bool world=false, const float x=1.0f, const float y=1.0f, const float z=1.0f, const float r=1.0f) {
  // Open file
  std::fstream file_op(filename.c_str(), std::ios::out);

  // Check if file is open
  if (file_op.is_open()) {
    // Write path to file
    for (pathtype::const_iterator it = path.begin(); it != path.end(); ++it) {
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
bool loadPathFromFile (const std::string & inputfile, pathtype & path) {
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

#endif
