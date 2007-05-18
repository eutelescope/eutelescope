// -*- mode: c++; mode: auto-fill; mode: flyspell-prog -*-
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
   * @Version $Id: EUTELESCOPE.h,v 1.9 2007-05-18 07:57:43 bulgheroni Exp $
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

    //! Parameter key to store/recall the beam location
    static const char * BEAMLOCATION;

    //! Parameter key to store/recall the beam type
    static const char * BEAMTYPE;

    //! Parameter key to store/recall the beam energy
    static const char * BEAMENERGY;

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

    //! Parameter key to store/recall the user defined run header comment
    static const char * USERCOMMENT;
    
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

    // PARAMETER NAMES USED IN THE EVENT IMPLEMENTATION
    
    //! Constant used to get/set the event type 
    static const char * EVENTTYPE;


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

    //! Cluster separation algorithm with only flagging capability
    /*! The name for the cluster separation algorithm that is not
     *  really dividing merging clusters, but only flagging their
     *  quality. @see EUTelClusterSeparationProcessor
     */
    static const char * FLAGONLY;

    // Encoding strings

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
     *  "sensorID:5,clusterID:8,xSeed:12,ySeed:12,xCluSize:5,yCluSize:5:quality:5"
     *
     *  Note about cluster quality: this is a three bit flag to be
     *  used with the cluster quality enum.
     *
     *  @see ClusterQuality
     */
    static const char * CLUSTERDEFAULTENCODING;
  };

  //! Event type enum
  /*! This enumeration type has been introduced in order to
   *  distinguish the three main different kind of events available in
   *  a EUTelEventImpl. Those are the following:
   *
   *  \li <b>kBORE</b>: Beginning Of Run Event. This is a dummy event
   *  containing none of the real data collection. For the time being
   *  it is not used since the first event of the run can be
   *  identified precisely being the first after the run header
   *  record.
   *
   *  \li <b>kEORE</b>: End of Run Event. This is the last event of
   *  the run. If the DAQ software was able to terminate properly the
   *  run this event should be the last, otherwise it has to be
   *  appended "manually" using a sort of LCIO file fix.
   *
   *  \li <b>kDE</b>: Data Event. This is flag is used for all the
   *  other events in the run, the ones containing the data collection
   *  to be used for the analysis.
   *
   *  \li <b>kUNKNOWN</b>: type unknown or not set. This value has
   *  been introduced in order to extend the compatibility to non
   *  EUTelEventImpl event. In such case, in fact, asking for a not
   *  existing parameter will return 0.
   *
   *  @Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @Version $Id: EUTELESCOPE.h,v 1.9 2007-05-18 07:57:43 bulgheroni Exp $
   */
  enum EventType {
    kUNKNOWN  = 0,
    kBORE     = 1,
    kDE       = 2,
    kEORE     = 3
  };



  //! Cluster quality enum
  /*! This enum can be attached to a LCIO class describing a cluster
   *  or it can be inserted into the CellID describing the cluster
   *  collection. It is a five bit flag, that can be used to
   *  discriminate among different cluster qualities. This is because
   *  not all clusters passing the required cuts can be considered to
   *  be at the same quality levels. For example there are clusters
   *  centered on or so close to the detector edge that are
   *  incomplete. Those, even if they are passing the threshold for
   *  cluster identification they cannot be used for resolution
   *  studies, because their calculated position if biased by the
   *  missing pixels. The same apply for pixels with one missing pixel
   *  because was flagged as bad during the clustering processor.
   *
   *  Here a description of all allowed value of cluster quality and
   *  their meaning:
   *
   *  \li <b>kGoodCluster</b>: this flag is used to identify clusters
   *  having no problem at all.
   *
   *  \li <b>kIncompleteCluster</b>: this flag is used to identify
   *  clusters in which at least one pixel is missing because it was
   *  flagged as bad during the previous analysis processors.
   *
   *  \li <b>kBorderCluster</b>: this flag is used with clusters found
   *  to close to the detector edge and for that reason, their
   *  completeness cannot be garanted.
   *
   *  \li <b>kMergedCluster</b>: this flag is used to label clusters
   *  that have been tracked in a position so close to another cluster
   *  that charge sharing between the two cannot be excluded.
   *
   *  There are still two "not assigned" bits that can be used in the
   *  future to mark other different kind of bad quality clusters.
   *
   *  @Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   *  @Version $Id: EUTELESCOPE.h,v 1.9 2007-05-18 07:57:43 bulgheroni Exp $
   */ 
  
  enum ClusterQuality {
    kGoodCluster       = 0,
    kIncompleteCluster = 1L << 0,
    kBorderCluster     = 1L << 1,
    kMergedCluster     = 1L << 2
  };

  //! Cluster quality bit-wise AND operator
  /*! This is a convenience operator used to identify the reason of a
   *  non good quality clusters. Bad quality clusters may be so for more
   *  than one reason simultaneously. This operator is used in the
   *  identification of such reasons.
   *  
   *  @param a A cluster quality value
   *  @param b Another cluster quality value
   *  @return the bit wise and among @a a and @a b
   */
  ClusterQuality operator&(ClusterQuality a, ClusterQuality b);
    


  //! Cluster quality bit-wise OR operator
  /*! This is a crucial operator for ClusterQuality since, during the
   *  cluster search processor, a cluster maybe flagged with one or
   *  more than one "bad" qualities. For this reason, using this
   *  operator can allow to flag the same cluster with more than one
   *  bad qualities.
   *
   *  @param a A cluster quality value
   *  @param b Another cluster quality value
   *  @return the bit wise or among @a a and @a b
   */ 
  ClusterQuality operator|(ClusterQuality a, ClusterQuality b);


  //! Cluster quality self bit-wise OR operator 
  /*! @see operator|(ClusterQuality a, ClusterQuality b)
   *  @bug Not working!!!! Have a look on the web...
   */ 
  ClusterQuality operator|=(ClusterQuality a, ClusterQuality b);

}


#endif
