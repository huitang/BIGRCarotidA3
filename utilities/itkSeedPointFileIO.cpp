#ifndef __itkSeedPointFileIO_CPP
#define __itkSeedPointFileIO_CPP

#include "itkSeedPointFileIO.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include <iostream>
#include <fstream>
#include <math.h>

#include <itksys/SystemTools.hxx>
#include <boost/algorithm/string.hpp>

namespace itk {
SeedPointFileIO::SeedPointFileIO()
{
  m_FileName = "";
  m_Points.clear();
  m_ListValid = false;
  m_Verbose = true;
}

SeedPointFileIO::PointListType SeedPointFileIO::GetPoints()
{
  if ( !m_ListValid && m_FileName != "" ){
    this->ReadPointsFromFile();
  }
  return m_Points;
}

void SeedPointFileIO::ReadPointsFromFile()
{
  // Different possibilities to give the seedpoints, added for compatibility
  // with mevis lab. All in world coordinates
  // 1 as plain string, separated by whitespaces
  // 2 as itkimage file, nonzeros are converted
  // 3 as txtfile (check on extension, possibly with six values on each line
  //   but only first three are taken)
  // 4 as ml xml file
  //

  m_Points.clear();

  if (m_FileName.empty())
  {
    std::cout << "missing filename for seedpoints" << std::endl;
  }
  if (m_Verbose)
  {
    std::cout << "reading points from  \t" << m_FileName << std::endl;
  }


  // 1 string
  // 2 itkfile
  // 3 txtfile
  // 4 xmlfile

  bool succeeded(false);
  int  guess(0);

  // determine extension
  const std::string ext = itksys::SystemTools::GetFilenameLastExtension(m_FileName);

  if (ext == ".vtk" || ext == ".mhd" || ext == ".dcm" ||
      ext == ".VTK" || ext == ".MHD" || ext == ".DCM" ) guess = 2;
  if (ext == ".txt" || ext == ".TXT" ) guess = 3;
  if (ext == ".xml" || ext == ".XML" ) guess = 4;
  if (guess == 0) guess = 1;

  // try as string
  if (guess == 1)
  {
    InternalPrecision buf;
    std::istringstream ss(m_FileName);

    std::vector<InternalPrecision> tokens;
    while (ss >> buf) 
    {
      tokens.push_back(buf);
    }
    const unsigned int n = tokens.size()/3;
    if (n>=2)
    {
      for (unsigned int i=0; i<n; ++i)
      {
        PrecisionPointType q;

        q[0] = tokens[i*3];
        q[1] = tokens[i*3+1];
        q[2] = tokens[i*3+2];

        m_Points.push_back(q);
      }

      succeeded = true;
    }
    else
    {
      if (m_Verbose)
      {
        std::cout << "error interpreting as string " << std::endl;
      }
    }
  }

  // try itkfile
  if (guess == 2)
  {
    typedef ImageFileReader<BinaryImageType> ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(m_FileName);

    try
    {
      reader->Update();

      const BinaryImageType::SpacingType     ispacing = reader->GetOutput()->GetSpacing();
      const BinaryImageType::PointType       iorigin = reader->GetOutput()->GetOrigin();

      ImageRegionConstIteratorWithIndex<BinaryImageType> 
        it (reader->GetOutput(), reader->GetOutput()->GetLargestPossibleRegion());

      for (it.GoToBegin(); !it.IsAtEnd(); ++it)
      {
        if (it.Value() != NumericTraits< BinaryImageType::PixelType>::Zero)
        {
          PrecisionPointType q;
          for (unsigned int i=0; i<Dim; ++i)
          {
            q[i] = iorigin[i] + it.GetIndex()[i] * ispacing[i];
          }
          m_Points.push_back(q);
        }
      }

      succeeded = true;
    }
    catch(ExceptionObject &)
    {
      if (m_Verbose)
      {
        std::cout << "error interpreting as itk file " << std::endl;
      }
    } // end catch reader
  }

  // try as txtfile
  if (guess == 3)
  {
    std::ifstream seedfile (m_FileName.c_str());
    if (!seedfile.is_open())
    {
      if (m_Verbose)
      {
        std::cout << "error opening as txt file " << std::endl;
      }
    }
    else
    {
      succeeded = true;
      while (!seedfile.eof())
      {
        std::string line;
        std::getline(seedfile,line);
        if (!line.empty())
        {
          InternalPrecision buf;
          std::istringstream ss(line);
          std::vector<InternalPrecision> tokens;
          while (ss>>buf)
          {
            tokens.push_back(buf);
          }
          // world coordinates incl offset!
          if (tokens.size() == 3 || tokens.size() == 6)
          {
            PrecisionPointType q;
            q[0] = tokens[0];
            q[1] = tokens[1];
            q[2] = tokens[2];
            m_Points.push_back(q);
          }
        }
        else
        {
          std::getline(seedfile,line);
        }
      }
    }
  }

  // xml
  if (guess == 4)
  {

    std::ifstream seedfile(m_FileName.c_str());

    if (!seedfile.is_open())
    {
      if (m_Verbose)
      {
        std::cout << "error opening as xml file " << std::endl;
      }
    }
    else
    {
      // mevis        <MeVis-XML ...>
      // xmarker      <XMarkerList ... >
      // size         <ListSize>
      // coord        <pos>

      // convert to vector
      std::vector <std::string> text;
      while (!seedfile.eof())
      {
        std::string line;
        std::getline(seedfile,line);
        boost::to_lower(line);

        if (!line.empty())
        {
          text.push_back(line);
        }
      }
      // check for validity and store line with points
      int n1(0), n2(-1), n3(-2), n4(-3);
      std::vector<std::string> pl;

      for (unsigned int i=0; i<text.size(); ++i)
      {
        if (boost::find_first(text[i], "<mevis-xml"))
          n1 = i;
        if (boost::find_first(text[i], "<xmarkerlist"))
          n2 = i;
        if (boost::find_first(text[i], "</xmarkerlist>"))
          n3 = i;
        if (boost::find_first(text[i], "</mevis-xml"))
          n4 = i;
        if (boost::find_first(text[i], "<pos>"))
          pl.push_back(text[i]);
      }
      if (n4 >= n3 && n3 >= n2 && n2 >= n1)
      {
        for (unsigned int i=0; i<pl.size(); ++i)
        {
          boost::erase_all(pl[i], "<pos>");
          boost::erase_all(pl[i], "</pos>");
          boost::trim(pl[i]);

          InternalPrecision buf;
          std::istringstream ss(pl[i]);
          std::vector<InternalPrecision> tokens;
          while (ss>>buf)
          {
            tokens.push_back(buf);
          }
          // world coordinates 
          if (tokens.size() >= 3)
          {
            PrecisionPointType q;
            q[0] = tokens[0];
            q[1] = tokens[1];
            q[2] = tokens[2];
            m_Points.push_back(q);
          }
        }
        succeeded = true;
      }
      else
      {
        if (m_Verbose)
        {
          std::cout << "error interpreting as mevis-xml file " << std::endl;
        }
        succeeded = false;
      }
    } // end else opening xml
  } // end guess 4

  if (!succeeded)
  {
    std::cout << "ERROR: retrieving seed points ... " << std::endl; 
    m_ListValid = false;
  }
  else 
  {
    m_ListValid = true;
  }


  if (m_Verbose)
  {
    std::cout << "number of points    \t" << m_Points.size() << std::endl; 
    for (unsigned int i=0; i<m_Points.size(); ++i)
    {
      std::cout << "point " << i << "              \t" << m_Points[i] << std::endl;
    }
  }

}


}
#endif //__SeedPointFileIO_CPP
