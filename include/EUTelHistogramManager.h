/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
// #ifdef EXPERIMENTAL
#ifndef EUTELHISTOGRAMMANAGER_H
#define EUTELHISTOGRAMMANAGER_H

// personal includes ".h"

// marlin includes ".h"
#include "marlin/Exceptions.h"

// lcio includes <.h>

// system includes <>
#include <string>
#include <exception>

namespace eutelescope {

  //! Histogram information helper class
  /*! This class is used to store the information about histograms
   *  read from the XML information file.
   *  
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$
   */
  class EUTelHistogramInfo {
    
  public:
    //! Histogram name
    std::string _name;
    
    //! Histogram title
    std::string _title;

    //! Histogram type
    /*! I don't think it will be used very often, but it may be
     *  useful in the future. Possible values are:
     * 
     *  @li H1D, H2D, H3D, for histograms 1, 2 and 3 dimensional
     *
     *  @li P1D, P2D  for profiles 1 and 2 dimensional
     *
     *  @li C1D, C2D, C3D, for clouds 1, 2 and 3 dimensional
     */
    std::string _type;

    //! Number of bin along x
    int _xBin;

    //! Minimum value along x
    float _xMin;
    
    //! Maximum value along x
    float _xMax;

    //! Number of bin along y
    int _yBin;
    
    //! Minimum value along y
    float _yMin;
    
    //! Maximum value along y
    float _yMax;

    //! Number of bin along z
    int _zBin;
    
    //! Minimum value along z
    float _zMin;
    
    //! Maximum value along z
    float _zMax;
    
    //! Streamer for the EUTelHistogramInfo class
    /*! Mostly used for debug purposes, this operator can stream out
     *  all histogram information read from the XML input
     *  file. Depending on the histogram type, only the needed
     *  information are displayed.
     *
     *  @return The output stream
     *  @param os The input stream
     *  @param histoInfo The object to stream out.
     *  
     */ 
    friend std::ostream& operator<<(std::ostream& os, const EUTelHistogramInfo& histoInfo ) {
      os << "===============================================\n" 
	 << "| Name         " << histoInfo._name << std::endl;
      if ( histoInfo._title != "" ) os << "| Title        " << histoInfo._title << std::endl;
      os << "| Type         " << histoInfo._type << std::endl;
      
      if (  ( histoInfo._type == "H1D" ) || 
	    ( histoInfo._type == "H2D" ) ||
	    ( histoInfo._type == "H3D" ) || 
	    ( histoInfo._type == "P1D" ) ||
	    ( histoInfo._type == "P2D" ) )  {
	os << "| xBin         " << histoInfo._xBin << std::endl
	   << "| xMin         " << histoInfo._xMin << std::endl
	   << "| xMax         " << histoInfo._xMax << std::endl;
      }
      
      if ( ( histoInfo._type == "H2D" ) ||
	   ( histoInfo._type == "H3D" ) ||
	   ( histoInfo._type == "P2D" ) )  {
	os << "| yBin         " << histoInfo._yBin << std::endl
	   << "| yMin         " << histoInfo._yMin << std::endl
	   << "| yMax         " << histoInfo._yMax << std::endl;
      }
      
      if ( ( histoInfo._type == "H3D" ) ) {
	os << "| zBin         " << histoInfo._zBin << std::endl
	   << "| zMin         " << histoInfo._zMin << std::endl
	   << "| zMax         " << histoInfo._zMax << std::endl;
      }
      
      if ( ( histoInfo._type == "P1D" ) ) {
	os << "| yMin         " << histoInfo._yMin << std::endl
	   << "| yMax         " << histoInfo._yMax << std::endl;
      }
      
      if ( ( histoInfo._type == "P2D" ) ) {
	os << "| zMin         " << histoInfo._zMin << std::endl
	   << "| zMax         " << histoInfo._zMax << std::endl;
      } 
      
      return os;
    }

  };




  //! Histogram manager class
  /*! This is the first attempt to have a histogram manager. 
   *
   *  The basic ideas behind this manager class are the following:
   *
   *  @li Easiness to rebin: For the time being histograms are booked
   *  according to hardcoded values and to modify them the user needs
   *  to open the code, find the right place, change the values and
   *  recompile the project. An elegant way out but probably too
   *  simplistic would be to book only clouds instead of
   *  histograms. The main disadvantage with that is that ROOT
   *  does not support clouds. A smarter solution would be to have a
   *  file containing the histogram booking info the user can modify.
   *
   *  <b>Histogram booking info</b>
   *  These kinds of information are contained into a XML file that is
   *  parsed when the EUTelHistogramManager is initialized and for
   *  each entry found a EUTelHistoInfo is added to the list of
   *  available histograms.
   * 
   *  @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @version $Id$
   */ 
    
  class EUTelHistogramManager  {

  public:
    //! Default constructor with histogram info file name
    /*! This is the default constructor
     *  
     *  @param histoInfoFileName The histogram information file name
     */ 
    EUTelHistogramManager(std::string histoInfoFileName) : _histoInfoFileName(histoInfoFileName), _histoInfoMap() {;}
   
    //! Destructor
    /*! Deletes all the entries of the map since they all have been
     *  created with new
     */
    virtual ~EUTelHistogramManager();
    
    //! Initialize the histogram manager
    /*! This method initializes the histogram manager, loading and
     *  parsing the XML information file. For each histogram block
     *  identified in the XML file, a new EUTelHistoInfo entry is
     *  added to the local map using as key the histogram name.
     *  
     *  @return True if the initialization has been successful 
     *  @throw exception in case of problem reading the info file
     *  @throw ParseException in case of problem parsing the XML file
     */
    bool init() throw( std::exception, marlin::ParseException ) ;


    //! Get the histogram information corresponding to histoName
    /*! This method is used to return a pointer to a
     *  EUTelHistogramInfo object containing the information
     *  concerning the histogram named @c histoName
     *
     *  @param histoName The name of the histogram
     *  @return A pointer to a EUTelHistogramInfo object. Returns 0x0
     *  in case @c histoName is not in the map
     */
    EUTelHistogramInfo * getHistogramInfo(std::string histoName) const ;

  private:
    //! Histogram information file name
    /*! This is the name of the file containing the histogram booking
     *  information
     */ 
    std::string _histoInfoFileName;

    //! Histogram information map
    /*! This is the key point of the histogram manager. For each
     *  histogram block found in the XML file there is a corresponding
     *  entries in this map keyed after the histogram name. 
     * 
     *  When booking the histogram, the user should ask to the
     *  EUTelHistogramManager is any histogram information is
     *  available for that particular histogram.
     */
    std::map< std::string , EUTelHistogramInfo *> _histoInfoMap;


  };
 
}
#endif
// #endif
