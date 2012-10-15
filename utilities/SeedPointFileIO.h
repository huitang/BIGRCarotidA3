
#ifndef __SeedPointFileIO_H
#define __SeedPointFileIO_H

#include "itkMacro.h"
#include "itkPoint.h"
#include "itkImage.h"
#include <string>

namespace itk {

#define Dim 3

class SeedPointFileIO : public Object {

public:
  typedef float InternalPrecision;
  typedef Point<InternalPrecision,Dim>      PrecisionPointType;
  typedef std::vector<PrecisionPointType> PointListType;
  typedef Image<unsigned char, 3>         BinaryImageType;

  typedef SeedPointFileIO Self;
  typedef Object    Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  itkTypeMacro(SeedPointFileIO, Object);
  itkNewMacro(Self);

  itkSetMacro(FileName, std::string );
  itkGetMacro(FileName, std::string );
  itkSetMacro(Verbose, bool );
  itkGetMacro(Verbose, bool );
  PointListType GetPoints();

protected:
  SeedPointFileIO();
  ~SeedPointFileIO() {}

private:
  std::string m_FileName;
  PointListType m_Points;
  bool m_Verbose;
  bool m_ListValid;

  virtual void ReadPointsFromFile();

};

}
#endif //__SeedPointFileIO_H
