// -*- C++ -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELESCOPE_H
#define EUTELESCOPE_H


namespace eutelescope
{

  //! Global constants used in the Eutelescope package
  /*!
   * This class has only static data members used only to define global
   * constant to be used within the Eutelescope package. Please add here
   * whatever constant you want to use.  A typical useful of this class
   * is to define name of collection to be retrieved/saved from/to
   * files.
   *
   * @Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   * @Version $Id: EUTELESCOPE.h,v 1.6 2007-02-22 08:09:36 bulgheroni Exp $
   */

  class EUTELESCOPE
  {

  public:
    //! Default destructor. 
    /*! This is the default destructor, but it is actually a NO-OP
     *  since there is nothing to be destroyed.
     */
    virtual ~ EUTELESCOPE ()  {;  }
    
  public:

    // PARAMETER NAMES USED IN THE HEADER IMPLEMENTATION

    //! Parameter key to store/recall the header version number
    static const char * HEADERVERSION;

    //! Parameter key to store/recall the data type 
    /*! The value of DATATYPE can be set using one of the static const
     *  defined in this class. Use DAQDATA for data coming from a real
     *  acquisition; SIMULDATA for data produced by a simulation job;
     *  CONVDATA for data converted from another data format.
     *  
     *  @see EUTELESCOPE::DAQDATA, EUTELESCOPE::SIMULDATA, EUTELESCOPE::CONVDATA
     */
    static const char * DATATYPE;

    //! Parameter key to store/recall the number of events in the file
    static const char * NOOFEVENT;

    //! Parameter key to store/recall the date time string in human readable format
    static const char * DATETIME;

    //! Parameter key to store/recall the DAQ hardware name
    /*! The value of the DAQHWNNAME stored in the header of LCIO file
     *  can be one of the DAQ name provided in this class.  @see
     *  EUTELESCOPE::EUDRB, EUTELESCOPE::IPHCIMAGER,
     *  EUTELESCOPE::SUCIMAIMAGER
     */
    static const char * DAQHWNAME;

    //! Parameter key to store/recall the DAQ hardware version 
    static const char * DAQHWVERSION;

    //! Parameter key to store/recall the DAQ software name
    /*! This parameter can be set using one of the const string
     *  defined as well in this class (e.g. EUDAQ). 
     *  
     *  @see EUTELESCOPE::EUDAQ
     */
    static const char * DAQSWNAME;

    //! Parameter key to store/recall the DAQ software version
    static const char * DAQSWVERSION;

    //! Parameter key to store/recall the simulation software name
    static const char * SIMULSWNAME;

    //! Parameter key to store/recall the simulation software version
    static const char * SIMULSWVERSION;

    //! Parameter key to store/recall the geometry identification number
    static const char * GEOID;

    //! Parameter key to store/recall the number of detector in the file
    static const char * NOOFDETECTOR;

    //! Parameter key to store/recall the minimum X pixel
    static const char * MINX;

    //! Parameter key to store/recall the maximum X pixel
    static const char * MAXX;

    //! Parameter key to store/recall the minimum Y pixel
    static const char * MINY;

    //! Parameter key to store/recall the maximum Y pixel
    static const char * MAXY;

    //! Parameter key to store/recall the list of applied processors
    static const char * APPLIEDPROCESSOR;

    //! Parameter key to store/recall the intermediate file names
    static const char * INTERMEDIATEFILE;

    // standard name to be saved into the header

    //! Constant used to identify really acquired data
    static const char * DAQDATA;

    //! Constant used to identify simulated data
    static const char * SIMULDATA;

    //! Constant used to identify data converted from another data format
    static const char * CONVDATA;

    //! Constant used to identify the DAQ board designed by INFN
    static const char * EUDRB;

    //! Constant used to identify the DAQ board developed in IPHC/LEPSI
    static const char * IPHCIMAGER;

    //! Constant used to identify the DAQ board developed by the SUCIMA collaboration
    static const char * SUCIMAIMAGER;

    //! Constant used to identify the DAQ software developed by Geneva team
    static const char * EUDAQ;

    // pixel flags

    //! Constant to identify good pixels
    static const int GOODPIXEL;

    //! Constant to identify bad pixels
    static const int BADPIXEL;

    //! Constant to identify hit pixels
    static const int HITPIXEL;

    // algorithm names

    //! Bad pixel masking algorithm identifier
    /*! @see EUTelPedestalNoiseProcessor::maskBadPixel() for a
     *  detailed description of the algorithm
     */
    static const char * NOISEDISTRIBUTION;
    
    //! Bad pixel masking algorithm identifier
    /*! @see EUTelPedestalNoiseProcessor::maskBadPixel() for a
     * detailed description of the algorithm.
     */
    static const char * ABSOLUTENOISEVALUE;


    //! Pedestal calculation algorithm identifier
    /*! This string is used to identify a pedestal calculation
     *  algorithm. @a MEANRMS means that the pedestal value is defined
     *  as the mean value of each single pixel signal distribution,
     *  while the noise is given by the distribution RMS. The
     *  algorithm is based on the use of std::vector's
     */ 
    static const char * MEANRMS;

    //! Pedestal calculation algorithm identifier
    /*! This string is used to identify a pedestal calculation
     *  algorithm. @a AIDAPROFILE means that pedestal and noise are
     *  determined taking advantages from the use of an
     *  AIDA::IProfile2D object. For each event, all pixel signals are
     *  filled into a 2d profile properly booked in the init()
     *  phase. When the loop is complete, the mean value and the noise
     *  of each pixel is dumped from the IProfile2D to pedestal and
     *  noise std::vector's. Of course, the use of it is limited to
     *  the cases in which MARLIN_USE_AIDA is defined. Otherwise, the
     *  algorithm will fall back to something, may be to so elegant,
     *  but definitely more standard (EUTELESCOPE::MEANRMS).
     */
    static const char * AIDAPROFILE;

    //! Fixed frame clustering algorithm
    /*! For a detailed description @see
     *  EUTelClusteringProcessor::_clusteringAlgo
     */
    static const char * FIXEDFRAME;

    //! Fixed weight algorithm for the pedestal / noise update
    /*! The name for the pedestal and noise update algorithm. @see
     *  EUTelUpdatePedestalNoiseProcessor
     */
    static const char * FIXEDWEIGHT;

    //! Default Tracker(Raw)Data encoding for full matrix
    /*! This constant string is used with CellIDEncoder to define the
     *  default encoding used for describe cells into a
     *  Tracker(Raw)Data object
     *
     *  "sensorID:5,xMin:12,xMax:12,yMin:12,yMax:12"
     */
    static const char * MATRIXDEFAULTENCODING;

    //! Default TrackerData encoding for cluster
    /*! This constant string is used with CellIDEncoder to define the
     *  default encoding used for describe cells into a clusters. This
     *  encoding is different from the one for complete matrices.
     *
     *  "sensorID:5,clusterID:8,xSeed:12,ySeed:12,xCluSize:5,yCluSize:5"
     */
    static const char * CLUSTERDEFAULTENCODING;
  };

}


#endif
