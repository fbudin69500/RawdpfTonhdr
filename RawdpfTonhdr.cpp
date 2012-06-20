#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include "RawdpfTonhdrCLP.h"
#include <itkImageFileReader.h>
#include <itkImage.h>
#include <cmath>

#define ZERO 0.0000001

int main( int argc , char* argv[] )
{
  PARSE_ARGS ;
  if( index.size() != 3 )
  {
    std::cerr << "Index has to have the following structure: first,last,step" << std::endl ;
  }
  //open dpf file
  std::ifstream dpffile( dpfFileName.c_str() ) ;
  if( dpffile.fail() )
  {
    std::cerr << "Could not open " << dpfFileName << std::endl ;
    return 1 ;
  }
  //Check nhdr extension
  if( outputnhdrFileName.compare( outputnhdrFileName.size() - 5, 5 ,".nhdr") )
  {
    outputnhdrFileName = outputnhdrFileName + ".nhdr" ;
  }
  //Check if raw is gzip
  bool gzip = false ;
  if( !rawFileName.compare(rawFileName.size() - 3 , 3 ,".gz") )
  {
    gzip = true ;
  }
  
  //create nhdr file
  std::ofstream nhdr( outputnhdrFileName.c_str() ) ;
  if( nhdr.fail() )
  {
    std::cerr << "Could not create " << outputnhdrFileName << std::endl ;
    return 1 ;
  }
  //Creates variables
  int size[ 3 ] ;
  double spacing[ 3 ] ;
  int bvalue = 0 ;
  std::vector< std::vector< double > > gradient ;
  //Reads dpf and gets info
  std::string line ;
  std::string trash ;
  size_t pos = std::string::npos ;
  do
  {
    std::getline( dpffile , line ) ;
    std::istringstream iline( line ) ;
    std::string buffer ;
    iline >> buffer ;
    if( !buffer.compare( 0 , 11 , "ImageWidth:" ) )
    {
      std::istringstream width( line ) ;
      width >> trash >> size[ 0 ] ;
    }
    else if( !buffer.compare( 0 , 12 , "ImageHeight:" ) )
    {
      std::istringstream height( line ) ;
      height >> trash >> size[ 1 ] ;
    }
    else if( !buffer.compare( 0 , 12 , "ImageSlices:" ) )
    {
      std::istringstream slices( line ) ;
      slices >> trash >> size[ 2 ] ;
    }
    else if( line.find( "PixelSize(X):" ) != std::string::npos )
    {
      pos = line.find( "PixelSize(X):" ) ;
      std::string sub = line.substr( pos , line.size() - 1 ) ;
      std::istringstream spacingx( sub ) ;
      spacingx >> trash >> spacing[ 0 ] ;
    }
    else if( line.find( "PixelSize(Y):" ) != std::string::npos )
    {
      pos = line.find( "PixelSize(Y):" ) ;
      std::string sub = line.substr( pos , line.size() - 1 ) ;
      std::istringstream spacingy( sub ) ;
      spacingy >> trash >> spacing[ 1 ] ;
    }
    else if( !buffer.compare( 0 , 15 , "SliceThickness:" ) )
    {
      std::istringstream slices( line ) ;
      slices >> trash >> spacing[ 2 ] ;
    }
    else if( !buffer.compare( 0 , 8 , "B-Value:" ) )
    {
      std::istringstream bval( line ) ;
      bval >> trash >> bvalue ;
    }
    else if( !buffer.compare( 0 , 23 , "Begin_Of_Gradient_Table" ) )
    {
      gradient.clear() ;
      do
      {
        std::getline( dpffile , line ) ;
        std::istringstream iline2( line ) ;
        iline2 >> buffer ;
        std::istringstream grad( line ) ;
        std::vector< double > temp_gradient( 3 ) ;
        grad >> trash >> temp_gradient[ 0 ] >> trash >> temp_gradient[ 1 ] >> trash >> temp_gradient[ 2 ] ;
        if( !notnormalized )
        {
          double norm = sqrt( temp_gradient[ 0 ] * temp_gradient[ 0 ]
                            + temp_gradient[ 1 ] * temp_gradient[ 1 ]
                            + temp_gradient[ 2 ] * temp_gradient[ 2 ]
                            ) ;
          if( norm > ZERO )
          {
            for( int j = 0 ; j < 3 ; j++ )
            {
              temp_gradient[ j ] /= norm ;
            }
          }
        }
        gradient.push_back( temp_gradient ) ;
      }while( buffer.compare( 0 , 21 , "End_Of_Gradient_Table" )
              && !dpffile.bad() && !dpffile.eof() && !dpffile.fail()
            ) ;
    }
  }while( !dpffile.bad() && !dpffile.eof() && !dpffile.fail() ) ;
  gradient.pop_back() ;
  
  //Prints info
  if( verbose )
  {
     std::cout << "Image Width: " << size[ 0 ] << std::endl ;
     std::cout << "Image Height: " << size[ 1 ] << std::endl ;
     std::cout << "Image Slices: " << size[ 2 ] << std::endl ;
     std::cout << "Pixel Size (X): " << spacing[ 0 ] << std::endl ;
     std::cout << "Pixel Size (Y): " << spacing[ 1 ] << std::endl ;
     std::cout << "Slice Thickness: " << spacing[ 2 ] << std::endl ;
     if( gzip )
     {
        std::cout << "encoding: gzip" << std::endl ;
     }
     else
     {
        std::cout << "encoding: raw" << std::endl ;
     }
     std::cout << "Data File: " << rawFileName << std::endl ;
     if( !b0dwi.compare( "DWI" ) )
     {
        std::cout << "B-Value: " << bvalue << std::endl ;
        for( unsigned int i = 0 ; i < gradient.size() ; i++ )
        {
           std::cout << "Gradient " << i 
                 << " : " << gradient[ i ][ 0 ]
                 << " " << gradient[ i ][ 1 ]
                 << " " << gradient[ i ][ 2 ] << std::endl ;
        }
     }
  }
  
  //writes dwi nhdr
  nhdr << "NRRD0005" << std::endl ;
  nhdr << "# Complete NRRD file format specification at:" << std::endl ;
  nhdr << "# http://teem.sourceforge.net/nrrd/format.html" << std::endl ;
  nhdr << "type: " << type << std::endl ;
  if( !b0dwi.compare( "DWI" ) )
  {
    nhdr << "dimension: 4" << std::endl ;
  }
  else
  {
    nhdr << "dimension: 3" << std::endl ;
  }
  nhdr << "space: " << space << std::endl ;
  nhdr << "sizes: " << size[ 0 ] << " "
       << size[ 1 ] << " " << size[ 2 ] ;
  if( !b0dwi.compare( "DWI" ) )
  {
       nhdr << " " << gradient.size() ;
  }
  nhdr << std::endl ;
  nhdr << "space directions: (" << spacing[ 0 ]
       << ",0,0) (0," << spacing[ 1 ]
       <<",0) (0,0," << spacing[ 2 ]
        << ")" ;
  if( !b0dwi.compare( "DWI" ) )
  {
     nhdr << " none" ;
  }
  nhdr << std::endl ;
  nhdr << "kinds: space space space" ;
  if( !b0dwi.compare( "DWI" ) )
  {
    nhdr << " list" ;
  }
  nhdr << std::endl ;
  nhdr << "endian: " << endian << std::endl ;
  if( gzip )
  {
    nhdr << "encoding: gzip" << std::endl ;
  }
  else
  {
    nhdr << "encoding: raw" << std::endl ;
  }
  nhdr << "space origin: (0,0,0)" << std::endl ;
  nhdr << "data file: " << rawFileName ;
  if( !b0dwi.compare( "DWI" ) )
  {
    nhdr << " " << index[ 0 ] << " " << index[ 1 ] << " " << index[ 2 ] ;
  }
  nhdr << std::endl ;
  if( !b0dwi.compare( "DWI" ) )
  {
    nhdr << "modality:=DWMRI" << std::endl ;
    nhdr << "measurement frame: (1,0,0) (0,1,0) (0,0,1)" << std::endl ;
    nhdr << "DWMRI_b-value:=" << bvalue << std::endl ;
    for( unsigned int i = 0 ; i < gradient.size() ; i++ )
    {
      nhdr << "DWMRI_gradient_" <<  std::setfill('0') << std::setw( 4 ) << i << ":="
           << gradient[ i ][ 0 ] << "     " << gradient[ i ][ 1 ] << "      "
           << gradient[ i ][ 2 ] << std::endl ;
    }
  }
  nhdr.close() ;
  dpffile.close() ;
  
  if( check )
  {
     if( verbose )
     {
        std::cout << "Trying to open the newly created file" << std::endl ;
     }
     typedef itk::Image< unsigned char , 3 > ImageType ;
     typedef itk::ImageFileReader< ImageType > ReaderType ;
     ReaderType::Pointer reader = ReaderType::New() ;
     reader->SetFileName( outputnhdrFileName.c_str() ) ;
     try
     {
       reader->Update() ;
     }
     catch ( itk::ExceptionObject & err )
     {
        std::cerr << "  Caught an ITK exception: " << std::endl ;
        std::cerr << err << " " << __FILE__ << " " << __LINE__ << std::endl ;
        return 1 ;
     }
     std::cout << "New file opened successfully" << std::endl ;
  }
  return 0 ;
}
