#ifndef __utils_H
#define __utils_H

#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <set>

// List type to store all vessels in
typedef std::vector< std::string > vesselListType;

// Define column numbers for different results
// These are the positions in which they are defined in the
// comparison files
#define COLUMN_NAMEGOLD 0
#define COLUMN_NAMEPATH 1
#define COLUMN_BIN 2
#define COLUMN_PERCENTAGE 3
#define COLUMN_ACCURACY 4
#define COLUMN_NUMPOINTSWITHIN 5
#define COLUMN_NUMPOINTSALL 6
#define COLUMN_AREAWITHIN 7
#define COLUMN_LENGTHGOLDWITHIN 8
#define COLUMN_NUMGOLDWITHIN 9
#define COLUMN_DICE 10

// Column offset for multiple bins
#define COLUMN_BINOFFSET 9

// Parameter positions in the file
// This will be the index in the file array after splitting the file on spaces
#define FILEPOS_ALPHA 1
#define FILEPOS_BETA 2
#define FILEPOS_C 3
#define FILEPOS_G1 4
#define FILEPOS_G2 5
#define FILEPOS_G3 6
#define FILEPOS_SCALES 7
#define FILEPOS_SIGMOIDA 8
#define FILEPOS_SIGMOIDB 9

// Parameter positions in the filename for Li evaluation
#define FILEPOS_LAMBDA 1
#define FILEPOS_RSTEP 2
#define FILEPOS_OMEGA 3
#define FILEPOS_LAMBDA1 4
#define FILEPOS_LAMBDA2 5

// This is the minimal amount of tokens we expect in every row of the
// comparison files
#define MIN_NUM_TOKENS 10

// Copyright year
#define COPYRIGHT_YEAR 2009

// List type to store randomly chosen vessels in
typedef std::set< std::string > randomVesselListType;
// Histogram type for percentages
typedef std::map< float, int > histogramType;

// Struct to store specific parameter settings and results in
struct parameterResultList {
  // Frangi alpha setting
  std::string frangiAlpha;
  // Frangi beta setting
  std::string frangiBeta;
  // Frangi c setting
  std::string frangiC;
  // Scales used
  std::string frangiScales;
  // Normalization factor scale 1
  std::string g1;
  // Normalization factor scale 2
  std::string g2;
  // Normalization factor scale 3
  std::string g3;
  
  // Intensity sigmoid position setting
  std::string intensitySigmoidA;
  // Intensity sigmoid steepness setting
  std::string intensitySigmoidB;

  // Total percentage in vessel
  double totalPercentage;

  // Sum, squared sum, mean and stddev of evaluated percentages in vessel
  double sumPercentage;
  double sumPercentageSquared;
  double meanPercentage;
  double stdDevPercentage;

  // Accuracy for correct parts
  double accuracy;

  // Sum of all areas corresponding to correct parts
  double sumArea;
  // Length of gold standard for correct parts
  double sumLengthGoldWithin;

  // Sum of num points within vessel
  double sumNumInTube;
  // Sum of all points of automatically found paths
  double sumNumAll;

  // Total number of vessels evaluated
  int num;

  // Histogram of percentages
  histogramType percentagesHistogram;

  // Constructor for initialization of values
  parameterResultList():    
    frangiAlpha(""),    
    frangiBeta(""),    
    frangiC(""),
    frangiScales(""),
    g1(""),
    g2(""),
    g3(""),
    intensitySigmoidA(""),
    intensitySigmoidB(""),
    totalPercentage(0.0),
    sumPercentage(0.0),
    sumPercentageSquared(0.0),
    meanPercentage(0.0),
    stdDevPercentage(0.0),

    accuracy(0.0),

    sumArea(0.0),
    sumLengthGoldWithin(0.0),

    sumNumInTube(0.0),
    sumNumAll(0.0),

    num(0) {}

  // Comparison operator. Structs are compared on totalPercentage
  bool operator< (const parameterResultList &rhs) const {
    return ( totalPercentage < rhs.totalPercentage );
  }
};

// Struct to store specific parameter settings and results in
struct parameterResultListLi {
  // Lambda setting (weighting term in r-dimension)
  std::string lambda;
  // Radius step setting
  std::string radiusStep;
  // Omage setting
  std::string omega;
  // Lambda 1 setting
  std::string lambda1;
  // Lambda 2 setting
  std::string lambda2;
  
  // Total percentage in vessel
  double totalPercentage;

  // Sum, squared sum, mean and stddev of evaluated percentages in vessel
  double sumPercentage;
  double sumPercentageSquared;
  double meanPercentage;
  double stdDevPercentage;

  // Accuracy for correct parts
  double accuracy;

  // Sum of all areas corresponding to correct parts
  double sumArea;
  // Length of gold standard for correct parts
  double sumLengthGoldWithin;

  // Sum of num points within vessel
  double sumNumInTube;
  // Sum of all points of automatically found paths
  double sumNumAll;

  // Total number of vessels evaluated
  int num;

  // Histogram of percentages
  histogramType percentagesHistogram;

  // Constructor for initialization of values
  parameterResultListLi():    
    lambda(""),    
    radiusStep(""),    
    omega(""),
    lambda1(""),
    lambda2(""),
    
    totalPercentage(0.0),
    sumPercentage(0.0),
    sumPercentageSquared(0.0),
    meanPercentage(0.0),
    stdDevPercentage(0.0),

    accuracy(0.0),

    sumArea(0.0),
    sumLengthGoldWithin(0.0),

    sumNumInTube(0.0),
    sumNumAll(0.0),

    num(0) {}

  // Comparison operator. Structs are compared on totalPercentage
  bool operator< (const parameterResultListLi &rhs) const {
    return ( totalPercentage < rhs.totalPercentage );
  }
};

// Type to store the n top results in
typedef std::multiset< parameterResultList > parameterTopType;
// Type to store the n top results for every bin in
typedef std::vector< parameterTopType > parameterTopForBinType;
// Type to store the n top results for every bin and number of vessels in
typedef std::vector< parameterTopForBinType > parameterTopForBinAndExpType;
// Type to store the n top results for every bin, number of vessels and 
// experiment number in
typedef std::vector< parameterTopForBinAndExpType > 
                                          parameterTopForBinExpAndVesselType;

// Li comparison types
// ============================
// Type to store the n top results in
typedef std::multiset< parameterResultListLi > parameterTopTypeLi;
// Type to store the n top results for every bin in
typedef std::vector< parameterTopTypeLi > parameterTopForBinTypeLi;
// Type to store the n top results for every bin and number of vessels in
typedef std::vector< parameterTopForBinTypeLi > parameterTopForBinAndExpTypeLi;
// Type to store the n top results for every bin, number of vessels and 
// experiment number in
typedef std::vector< parameterTopForBinAndExpTypeLi > 
parameterTopForBinExpAndVesselTypeLi;

// Print a program header
void programHeader (const std::string & name, const std::string author, 
                    const int copyright=2009) {
  std::cout << std::endl << "_________________________________________________" 
            << std::endl
            << std::endl << " " << name << " (c) " << copyright << ", " 
            << author 
            << std::endl << "_________________________________________________" 
            << std::endl << std::endl;
}

void tokenize(const std::string & str, std::vector<std::string> & tokens, 
              const std::string& delimiters = " ") {
  // Skip delimiters at beginning.
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

  while (std::string::npos != pos || std::string::npos != lastPos) {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}

void replace1st (std::string & str, const std::string & find, const std::string & replace) {
  std::string::size_type loc = str.find (find, 0);
  if (loc != std::string::npos) {
    str.erase(loc, find.length());
    str.insert(loc, replace);
  }
}

// Gaussian function
inline float gaussian1D (const float x, const float sigma) {
  return exp (-(x*x)/(2*sigma*sigma));
}

// Sigmoid function
inline float sigmoid (const float x, const float a, const float b) {
  return 1.0f / (1.0f + exp((a-x)/b));
}

// Load all vessels from a comparison file
void getVesselsFromFile (const std::string file, vesselListType & vesselList) {
  vesselList.clear();
  // Open file
  std::fstream file_op(file.c_str(), std::ios::in);
  // Check if file is open
  if (file_op.is_open()) {
    // Get filename from line
    std::string line;
    while (std::getline(file_op, line)) {
      std::vector<std::string> tokens;
      tokenize(line, tokens);
      if (tokens.size()>1) {
        // Find filename
        std::vector<std::string> tokensPath;
        tokenize(tokens[1], tokensPath, "/");
        // Find datasetnumber and vessel
        std::vector<std::string> tokensFile;
        tokenize(tokensPath[tokensPath.size()-1], tokensFile, "_");
        std::string dataset = tokensFile[tokensFile.size()-2];
        std::vector<std::string> tokensVessel;
        // Remove extension
        tokenize(tokensFile[tokensFile.size()-1], tokensVessel, ".");
        std::string vessel = tokensVessel[0];
        std::stringstream datasetVessel;
        datasetVessel << dataset << "_" << vessel;
        vesselList.push_back(datasetVessel.str());
      }
    }
  }
}

// Load all vessels from a comparison file
void getVesselsFromFileLi (const std::string file, vesselListType & vesselList) {
  vesselList.clear();
  // Open file
  std::fstream file_op(file.c_str(), std::ios::in);
  // Check if file is open
  if (file_op.is_open()) {
    // Get filename from line
    std::string line;
    while (std::getline(file_op, line)) {
      std::vector<std::string> tokens;
      tokenize(line, tokens);
      if (tokens.size()>1) {
        // Find filename
        std::vector<std::string> tokensPath;
        tokenize(tokens[1], tokensPath, "/");
        std::string dataset = tokensPath[tokensPath.size()-3];
        std::string vessel = tokensPath[tokensPath.size()-2];
        std::stringstream datasetVessel;
        datasetVessel << dataset << "_" << vessel;
        vesselList.push_back(datasetVessel.str());
      }
    }
  }
}

#endif
