// -*- C++ -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELRUNHEADERIMPL_H
#define EUTELRUNHEADERIMPL_H

// personal includes ".h"
#include "EUTELESCOPE.h"

// marlin includes ".h"

// lcio includes <.h>
#include <lcio.h>
#include <IMPL/LCRunHeaderImpl.h>

// system includes <>


namespace eutelescope
{

  //! Implementation of the Run Header for the EUDET telescope.
  /*! This is used within the EUDET JRA1 collaboration to store the
   *  run header into the LCIO files produced both by the system DAQ
   *  and by simulation software. This class is inheriting from
   *  LCRunHeaderImpl and in order to maintain compatibility with the
   *  rest of the LCIO framework, it does not add any other data
   *  members, but it simply provide methods to set other specific
   *  name parameters inside the LCParameter object. The following
   *  other parameters have been defined: 
   *
   *  \li <b>HeaderVersion</b>: a float number representing the
   *  version of this header class. Standard version number v01-23-04
   *  are converted in a float number like 1.2304.
   *
   *  \li <b>DataType</b>: this string is used to distinguish real
   *  data run from simulation. In the EUTELESCOPE class few static
   *  constant string are available for this purpose:
   *  EUTELESCOPE::DAQDATA, EUTELESCOPE::SIMULDATA and
   *  EUTELESCOPE::CONVDATA. The last one has to be used for real data
   *  or simulated data files that were generated and saved with
   *  another data file.
   *
   *  \li <b>EventNumber</b>: this is the total number of events
   *  saved into the file. Due to I/O data format, it is possible to
   *  know a priori how many events are stored into the file. For
   *  some specific analysis processes like pedestal calculation and
   *  eta function estimation, just quote two of them, it is
   *  compulsory to have a processor performing multiple loops over
   *  all events. This has been made possible by the use of a
   *  RewindDataFileException that can only be thrown within the
   *  processEvent(LCEvent*) callback. So it means that you need to
   *  call it before going into end(): if you know how many events
   *  are stored in your file then when isLastEvent() becomes true
   *  throw the exception. If you don't know it, you are doom to fall
   *  in the end() without any possibility to go back at the
   *  beginning of the loop.
   * 
   *  \li <b>DateTime</b>: this is a string showing in a human
   *  readable format the current date and time. Calling the set
   *  method, the current data is automatically saved.
   *
   *  \li <b>DAQHWName</b>: this string is used to identify the
   *  hardware part of the DAQ. In the EUTELESCOPE class few static
   *  constant string are available for this purpose:
   *  EUTELESCOPE::EUDRB, EUTELESCOPE::IPHCIMAGER and
   *  EUTELESCOPE::SUCIMAIMAGER.
   *
   *  \li <b>DAQHWVersion</b>: this is a float number representing
   *  the version of the DAQ hardware system. This value is
   *  meaningfull only when the run is not a simulation. The
   *  conversion from version to float is the same of HeaderVersion
   *
   *  \li <b>DAQSWName</b>: this is a string used to identify the
   *  software part of the used DAQ. In the EUTELESCOPE class there is
   *  (for the time being) just one static string for this purpose:
   *  EUTELESCOPE::EUDAQ.
   *
   *  \li <b>DAQSWVersion</b>: this is a float number representing
   *  the version of the DAQ software. This value is meaningfull only
   *  when the run is not a simulation. The conversion from version to
   *  float is the same of HeaderVersion.
   *
   *  \li <b>SimulSWName</b>: this string parameter is used to
   *  record the name of the simulation program that was used to
   *  generate this file.
   *  
   *  \li <b>SimulSWVer</b>: this float number is used to eventually
   *  store also the version of the simulation software that has been
   *  used for data generation.
   *
   *  \li <b>GeoID</b>: this is an integer number used to link the
   *  current telescope geometrical configuration with an entry in the
   *  geometry database. This is used during the reconstruction phase
   *  into a DB query to retrieve precise information about the
   *  detector positioning and alignement.
   *
   *  \li <b>NoOfDetector</b>: this int number represents the number
   *  of pixel detectors into the current telescope configuration. It
   *  might correspond to the number of planes but not necessarily,
   *  since more than one detector can also be used to make a
   *  plane. As a consequence, it also represent the number of
   *  Tracker(Raw)Data objects are saved into each collection. When a
   *  detector is sub-divided into channels (as Mimo*2 as two
   *  channels) the user may want to save each of them separately into
   *  a Tracker(Raw)Data object, because single detector
   *  characterization is done at the level of a single
   *  Tracker(Raw)Data object. 
   *  <b>Note</b>: This information is of course repeated
   *  into the geometry DB, but, at least, for only monitoring,
   *  detector debug and single plane characterization is better
   *  having this information also localy.
   *
   *  \li <b>MinX</b>: this is an IntVec containing NoOfDetector
   *  numbers. Each of this represents the pixel number of the minimum
   *  pixel along the horizontal direction starting from 0. When a
   *  TrackerRawData contains all the pixels of a sensor, this number
   *  is bound to be zero. But when the TrackerRawData contains just
   *  one channel/submatrix of a sensor, it can be different from
   *  zero. See also the <b>Note</b> paragraph of the NoOfDetector
   *  description.
   *
   *  \li <b>MaxX</b>: this IntVec is the counterpart of MinX. Each
   *  component represents the pixel number of the last pixel along
   *  the horizontal direction belonging to the sensor. To obtain the
   *  number of pixel along this direction one should make (MaxX -
   *  MinX) + 1 to take into account also the first (0) pixel. See
   *  also the <b>Note</b> paragraph of the NoOfDetector
   *  description.
   *
   *  \li <b>MinY</b> and <b>MaxY</b> are playing the same role of
   *  MinX and MaxX but along the vertical direction.
   *  
   *  \li <b>ProcessorList</b>: this is a StringVec used to store
   *  which processors have been applied to this lcio file. If this is
   *  a "just saved" file from the DAQ, the list would be empty.
   *
   *  \li <b>IntermediateFile</b>: during the analysis procedure that
   *  can involve several steps/processors, it may be usefull to save
   *  the intermediate results into a file. This list of string
   *  contains the intermediate file used to produce this lcio
   *  file. Again this is usefull to reconstruct the file history.
   *
   * @Author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   * @Version $Id: EUTelRunHeaderImpl.h,v 1.2 2007-02-09 20:34:26 bulgheroni Exp $
   * 
   */

  class EUTelRunHeaderImpl:public IMPL::LCRunHeaderImpl
  {

  public:

    //! Default constructor
    EUTelRunHeaderImpl ();

    //! Destructor
    virtual ~ EUTelRunHeaderImpl ()
    { /* NO-OP */ ;
    }

    //! Set the header version number
    /*! this is a float number representing the version of this header
     *  class.Standard version number v01 - 23 - 04 are converted in a
     *  float number like 1.2304.
     *
     *  @param ver The version number the user wants to set
     */
    virtual void setHeaderVersion (float ver);
     
    //! Set the number of events in the file
    /*! this is an integer number equal to the number of events in the
     *  file.
     * 
     *  @par num The number of events.
     */
    virtual void setNoOfEvent(int num);

    //! Set the data type
    /*! this string is used to distinguish real
     *  data run from simulation. In the EUTELESCOPE namespace few
     *  static constant string are available for this purpose:
     *  EUTELESCOPE::DAQDATA and EUTELESCOPE::SIMULDATA.
     *
     *  @param type The type of data contained into the file
     */
    virtual void setDataType (std::string type);

    //! Set the current date and time of the day
    /*! this string is used to distinguish real
     *  data run from simulation. In the EUTELESCOPE namespace few
     *  static constant string are available for this purpose:
     *  EUTELESCOPE::DAQDATA and EUTELESCOPE::SIMULDATA.
     *
     *  @todo Consider the possibility to add another method having as
     *  an input parameter a UTIL::LCTime pointer. This will allow the
     *  user to set the date and time of the current file to any
     *  value.
     */
    virtual void setDateTime ();

    //! Set the DAQ hardware name
    /*! this string is used to identify the hardware part of the
     *  DAQ. In the EUTELESCOPE class few static constant string are
     *  available for this purpose: EUTELESCOPE::EUDRB and
     *  EUTELESCOPE::IPHCIMAGER.
     *
     *  @param name The DAQ hardware name
     */
    virtual void setDAQHWName (std::string name);

    //! Set the DAQ hardware version
    /*! this is a float number representing
     *  the version of the DAQ software. This value is meaningfull only
     *  when the run is not a simulation. The conversion from version to
     *  float is the same of HeaderVersion.
     *
     *  @param ver The DAQ hardware version 
     */
    virtual void setDAQHWVersion (float ver);

    //! Set the DAQ software name
    /*! this is a string used to identify the
     *  software part of the used DAQ. In the EUTELESCOPE class there is
     *  (for the time being) just one static string for this purpose:
     *  EUTELESCOPE::EUDAQ.
     *
     *  @param name The name of the DAQ software
     */
    virtual void setDAQSWName (std::string name);

    //! Set the DAQ version number
    /*! this is a float number representing
     *  the version of the DAQ software. This value is meaningfull only
     *  when the run is not a simulation. The conversion from version to
     *  float is the same of HeaderVersion.
     *
     *  @param ver The version of the DAQ software
     */
    virtual void setDAQSWVersion (float ver);

    //! Set the simulation software name
    /*! this string parameter is used to
     *  record the name of the simulation program that was used to
     *  generate this file.
     *
     *  @param name The name of the simulation program used for data generation
     */
    virtual void setSimulSWName (std::string name);

    //! Set the simulation software version
    /*! this float number is used to eventually
     *  store also the version of the simulation software that has been
     *  used for data generation.
     *
     *  @param ver The software version number
     */
    virtual void setSimulSWVersion (float ver);

    //! Set the geometry identification number
    /*! this is an integer number used to link the
     *  current telescope geometrical configuration with an entry in the
     *  geometry database. This is used during the reconstruction phase
     *  into a DB query to retrieve precise information about the
     *  detector positioning and alignement.
     *
     *  @param id The identification number of the current geometry
     */
    virtual void setGeoID (int id);

    //! Set the number of detector in this run
    /*! this int number represents the number of pixel detectors into
     *  the current telescope configuration. It might correspond to
     *  the number of planes but not necessarily, since more than one
     *  detector can also be used to make a plane. As a consequence,
     *  it also represent the number of Tracker(Raw)Data objects are
     *  saved into each collection. When a detector is sub-divided
     *  into channels (as Mimo*2 as two channels) the user may want to
     *  save each of them separately into a Tracker(Raw)Data object,
     *  because single detector characterization is done at the level
     *  of a single Tracker(Raw)Data object.  <b>Note</b>: This
     *  information is of course repeated into the geometry DB, but,
     *  at least, for only monitoring, detector debug and single plane
     *  characterization is better having this information also
     *  localy.
     *
     *  @param num The number of detectors in this file
     */
    virtual void setNoOfDetector (int num);

    //! Set the MinX vector
    /*! this is an IntVec containing NoOfDetector numbers. Each of
     *  this represents the pixel number of the minimum pixel along the
     *  horizontal direction starting from 0. When a TrackerRawData
     *  contains all the pixels of a sensor, this number is bound to be
     *  zero. But when the TrackerRawData contains just one
     *  channel/submatrix of a sensor, it can be different from
     *  zero. See also the <b>Note</b> paragraph of the NoOfDetector
     *  description.
     *
     *  @param xMin The std::vector containing the minimum pixel along X
     */
    virtual void setMinX (lcio::IntVec xMin);

    //! Set the MaxX vector
    /*! this IntVec is the counterpart of MinX. Each component
     *  represents the pixel number of the last pixel along the
     *  horizontal direction belonging to the sensor. To obtain the
     *  number of pixel along this direction one should make (MaxX -
     *  MinX) + 1 to take into account also the first (0) pixel. See
     *  also the <b>Note</b> paragraph of the NoOfDetector
     *  description. 
     *
     *  @param xMax The std::vector containing the maximum pixel along X
     */
    virtual void setMaxX (lcio::IntVec xMax);

    //! Set the MinY vector
    /*! @see EUTelRunHeaderImpl::setMinX
     *  @param yMin  The std::vector containing the minimum pixel along Y
     */
    virtual void setMinY (lcio::IntVec yMin);

    //! Set the MaxY vector
    /*! @see EUTelRunHeaderImpl::setMaxX
     *  @param yMax  The std::vector containing the maximum pixel along Y
     */
    virtual void setMaxY (lcio::IntVec yMax);

    //! Add a processor to the applied processor list
    /*  The analysis procedure of some input data usually requires
     *  that many processors have been applied sequentially. Saving
     *  this list it may be of interest when tring to reconstruct the
     *  history of a file.
     *  
     *  @param processor The name of the process to add to the list
     */
    virtual void addProcessor (std::string processor);

    //! Add an intermediate file to the list of intermediate file
    /*  The analysis procedure may be split in several intermediate
     *  steps with results saved into other files. Use this method to
     *  add a filename to the list of intermediate file of this one.
     * 
     *  @param file The intermediate file name to add to the list
     */
    virtual void addIntermediateFile (std::string file);

    //! return the header version
    inline float getHeaderVersion () const
    {
      return _params.getFloatVal (EUTELESCOPE::HEADERVERSION);
    }

    //! return the data type 
    inline std::string getDataType () const
    {
      return _params.getStringVal (EUTELESCOPE::DATATYPE);
    }

    //! return the number of events
    inline int getNoOfEvent() const 
    {
      return _params.getIntVal(EUTELESCOPE::NOOFEVENT);
    }

    //! return the date and time in a human readable format
    inline std::string getDateTime () const
    {
      return _params.getStringVal (EUTELESCOPE::DATETIME);
    }
    //! return the DAQ hardware name 
    inline std::string getDAQHWName () const
    {
      return _params.getStringVal (EUTELESCOPE::DAQHWNAME);
    }

    //! return the DAQ hardware version 
    inline float getDAQHWVersion () const
    {
      return _params.getFloatVal (EUTELESCOPE::DAQHWVERSION);
    }

    //! return the DAQ software name 
    inline std::string getDAQSWName () const
    {
      return _params.getStringVal (EUTELESCOPE::DAQSWNAME);
    }

    //! return the DAQ software version 
    inline float getDAQSWVersion () const
    {
      return _params.getFloatVal (EUTELESCOPE::DAQSWVERSION);
    }

    //! return the simulation software name
    inline std::string getSimulSWName () const
    {
      return _params.getStringVal (EUTELESCOPE::SIMULSWNAME);
    }

    //! return the simulation software version
    inline float getSimulSWVersion () const
    {
      return _params.getFloatVal (EUTELESCOPE::SIMULSWVERSION);
    }
    //! return the geometry identification number 
    inline int getGeoID () const
    {
      return _params.getIntVal (EUTELESCOPE::GEOID);
    }

    //! return the number of pixel detectors in the file
    inline int getNoOfDetector () const
    {
      return _params.getIntVal (EUTELESCOPE::NOOFDETECTOR);
    }

    //! return the vector of minimum pixel along X
    inline lcio::IntVec getMinX () const
    {
      lcio::IntVec v;
      return _params.getIntVals (EUTELESCOPE::MINX, v);
    }

    //! return the vector of maximum pixel along X
    inline lcio::IntVec getMaxX () const
    {
      lcio::IntVec v;
      return _params.getIntVals (EUTELESCOPE::MAXX, v);
    }

    //! return the vector of minimum pixel along Y
    inline lcio::IntVec getMinY () const
    {
      lcio::IntVec v;
      return _params.getIntVals (EUTELESCOPE::MINY, v);
    }

    //! return the vector of minimum pixel along Y
    inline lcio::IntVec getMaxY () const
    {
      lcio::IntVec v;
      return _params.getIntVals (EUTELESCOPE::MAXY, v);
    }
  };                           // end of EUTelRunHeaderImpl
}                               // eutelescope namespace

#endif // EUTELRUNHEADERIMPL
