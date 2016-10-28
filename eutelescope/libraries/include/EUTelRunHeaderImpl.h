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
#include <EVENT/LCRunHeader.h>
#include <IMPL/LCRunHeaderImpl.h>

// system includes <>


namespace eutelescope {

  //! Implementation of the Run Header for the EUDET telescope.
  /*! This is used within the EUDET JRA1 collaboration to store the
   *  run header into the LCIO files produced both by the system DAQ
   *  and by simulation software. This class is using the decorator
   *  pattern around the LCRunHeaderImpl. The following parameters
   *  have been defined:
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
   *  meaningful only when the run is not a simulation. The
   *  conversion from version to float is the same of HeaderVersion
   *
   *  \li <b>DAQSWName</b>: this is a string used to identify the
   *  software part of the used DAQ. In the EUTELESCOPE class there is
   *  (for the time being) just one static string for this purpose:
   *  EUTELESCOPE::EUDAQ.
   *
   *  \li <b>DAQSWVersion</b>: this is a float number representing
   *  the version of the DAQ software. This value is meaningful only
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
   *  detector positioning and alignment.
   *
   *  \li <b>BeamLocation</b>: this is a string field used to store
   *  where the test beam is taking place. It can be something like
   *  DESY, or CERN or it might be more detailed like CERN-H8. 
   *
   *  \li <b>BeamType</b>: this is another string field used to store
   *  the beam type, i.e. the kind of particle used (muons, electrons,
   *  hadron)
   *
   *  \li <b>BeamEnergy</b>: this is float number representing the
   *  beam enegery in GeV.
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
   *  having this information also locally.
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
   *  can involve several steps/processors, it may be useful to save
   *  the intermediate results into a file. This list of string
   *  contains the intermediate file used to produce this lcio
   *  file. Again this is useful to reconstruct the file history.
   *
   *  \li <b>UserComment</b>: a user defined string field. piipo
   *
   * @author Antonio Bulgheroni, INFN <mailto:antonio.bulgheroni@gmail.com>
   * @version $Id$
   * 
   */

  class EUTelRunHeaderImpl   {

  public:
    
    //! Default constructor
    EUTelRunHeaderImpl ( lcio::LCRunHeader * lcHeader ) : _lcHeader(NULL) { _lcHeader = dynamic_cast<IMPL::LCRunHeaderImpl*> (lcHeader) ; }

    //!Copy Constructor
    EUTelRunHeaderImpl(const EUTelRunHeaderImpl &z);
    
    //!Assignment Operator
    EUTelRunHeaderImpl& operator = (const EUTelRunHeaderImpl &z);
    
    //! Destructor
    virtual ~ EUTelRunHeaderImpl ()  { /* NO-OP */ ;  }

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
     *  @param num The number of events.
     */
    virtual void setNoOfEvent(int num);

    //! Set the data type
    /*! this string is used to distinguish real data run from
     *  simulation. In the EUTELESCOPE class few static constant
     *  string are available for this purpose: EUTELESCOPE::DAQDATA
     *  and EUTELESCOPE::SIMULDATA.
     *
     *  @param type The type of data contained into the file
     */
    virtual void setDataType (std::string type);

    //! Set the current date and time of the day
    /*! this string is used to distinguish real data run from
     *  simulation. In the EUTELESCOPE class few static constant
     *  string are available for this purpose: EUTELESCOPE::DAQDATA
     *  and EUTELESCOPE::SIMULDATA.
     *
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
     *  the version of the DAQ software. This value is meaningful only
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
     *  the version of the DAQ software. This value is meaningful only
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
     *  detector positioning and alignments.
     *
     *  @param id The identification number of the current geometry
     */
    virtual void setGeoID (int id);
    
    //! Set the beam location 
    /*! This is a string field representing where the beam test is
     *  taking place. In case of simulated data is meaningless. It can
     *  be something like "DESY" or "CERN-H8"
     *
     *  @param location The place where data were collected.
     */
    virtual void setBeamLocation (std::string location);

    //! Set beam type
    /*! This string is used to shortly describe the kind of beam used:
     *  like electrons, hadrons, ....
     *
     *  @param type The type of beam
     */
    virtual void setBeamType (std::string type);

    //! Set beam energy 
    /*! This float number represent the beam energy measured in GeV.
     * 
     *  @param energy Beam energy in GeV
     */ 
    virtual void setBeamEnergy (float energy);
    
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
     *  locally.
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

    //! Set the detector modality
    /*! This is the tag taken from the DAQ software to identify how
     *  the EUDRBs are configured. For example if they are working in
     *  RAW mode or in ZS.
     * 
     *  @param mode A string representing the mode of operation
     */
    virtual void setEUDRBMode(std::string mode);

    //! Set which detectors are in the telescope
    /*! This is the tag takend from the DAQ software to identify which
     *  sensors are in the telescope. In case of a mixed
     *  configuration, the processor will discover the correct order
     *
     *  @param det A string naming the sensors in the telescope
     */
    virtual void setEUDRBDet(std::string det);

    //! Add a processor to the applied processor list
    /*! The analysis procedure of some input data usually requires
     *  that many processors have been applied sequentially. Saving
     *  this list it may be of interest when trying to reconstruct the
     *  history of a file.
     *  
     *  @param processor The name of the process to add to the list
     */
    virtual void addProcessor (std::string processor);

    //! Add an intermediate file to the list of intermediate file
    /*! The analysis procedure may be split in several intermediate
     *  steps with results saved into other files. Use this method to
     *  add a filename to the list of intermediate file of this one.
     * 
     *  @param file The intermediate file name to add to the list
     */
    virtual void addIntermediateFile (std::string file);

    //! User defined string field
    /*! This field can be used to store a comment from the user or
     *  anything else the user wants to write in the run header
     *
     *  @param note The user note about the run
     */ 
    virtual void setUserComment(std::string note);

    //! return the header version
    inline float getHeaderVersion () const    {
      return _lcHeader->parameters().getFloatVal (EUTELESCOPE::HEADERVERSION);
    }

    //! return the data type 
    inline std::string getDataType () const     {
      return _lcHeader->parameters().getStringVal (EUTELESCOPE::DATATYPE);
    }

    //! return the number of events
    inline int getNoOfEvent() const     {
      return _lcHeader->parameters().getIntVal(EUTELESCOPE::NOOFEVENT);
    }

    //! return the date and time in a human readable format
    inline std::string getDateTime () const    {
      return _lcHeader->parameters().getStringVal (EUTELESCOPE::DATETIME);
    }

    //! return the DAQ hardware name 
    inline std::string getDAQHWName () const    {
      return _lcHeader->parameters().getStringVal (EUTELESCOPE::DAQHWNAME);
    }

    //! return the DAQ hardware version 
    inline float getDAQHWVersion () const    {
      return _lcHeader->parameters().getFloatVal (EUTELESCOPE::DAQHWVERSION);
    }

    //! return the DAQ software name 
    inline std::string getDAQSWName () const    {
      return _lcHeader->parameters().getStringVal (EUTELESCOPE::DAQSWNAME);
    }

    //! return the DAQ software version 
    inline float getDAQSWVersion () const    {
      return _lcHeader->parameters().getFloatVal (EUTELESCOPE::DAQSWVERSION);
    }

    //! return the simulation software name
    inline std::string getSimulSWName () const    {
      return _lcHeader->parameters().getStringVal (EUTELESCOPE::SIMULSWNAME);
    }

    //! return the simulation software version
    inline float getSimulSWVersion () const    {
      return _lcHeader->parameters().getFloatVal (EUTELESCOPE::SIMULSWVERSION);
    }
    
    //! return the geometry identification number 
    inline int getGeoID () const    {
      return _lcHeader->parameters().getIntVal (EUTELESCOPE::GEOID);
    }

    //! return the beam location
    inline std::string getBeamLocation() const {
      return _lcHeader->parameters().getStringVal(EUTELESCOPE::BEAMLOCATION);
    }

    //! return the number of pixel detectors in the file
    inline int getNoOfDetector () const    {
      return _lcHeader->parameters().getIntVal (EUTELESCOPE::NOOFDETECTOR);
    }

    //! return the vector of minimum pixel along X
    inline lcio::IntVec getMinX () const    {
      lcio::IntVec v;
      return _lcHeader->parameters().getIntVals (EUTELESCOPE::MINX, v);
    }

    //! return the vector of maximum pixel along X
    inline lcio::IntVec getMaxX () const    {
      lcio::IntVec v;
      return _lcHeader->parameters().getIntVals (EUTELESCOPE::MAXX, v);
    }

    //! return the vector of minimum pixel along Y
    inline lcio::IntVec getMinY () const    {
      lcio::IntVec v;
      return _lcHeader->parameters().getIntVals (EUTELESCOPE::MINY, v);
    }

    //! return the vector of minimum pixel along Y
    inline lcio::IntVec getMaxY () const    {
      lcio::IntVec v;
      return _lcHeader->parameters().getIntVals (EUTELESCOPE::MAXY, v);
    }

    //! return the user comment 
    inline std::string getUserComment() const {
      return _lcHeader->parameters().getStringVal(EUTELESCOPE::USERCOMMENT);
    }

    //! return the EUDRB operation mode
    inline std::string getEUDRBMode() const {
      if ( _lcHeader->parameters().getStringVal(EUTELESCOPE::EUDRBMODE) == "" ) return std::string("RAW3");
      return _lcHeader->parameters().getStringVal(EUTELESCOPE::EUDRBMODE);
    }

    //! return the EUDRB global detector
    inline std::string getEUDRBDet() const {
      if ( _lcHeader->parameters().getStringVal(EUTELESCOPE::EUDRBDET) == "" ) return std::string("MIMOTEL");
      return _lcHeader->parameters().getStringVal(EUTELESCOPE::EUDRBDET);
    }

    //! returns the LCRunHeaderImpl underlying object
    inline IMPL::LCRunHeaderImpl * lcRunHeader() { return  _lcHeader ; }

  private:

    IMPL::LCRunHeaderImpl * _lcHeader;

  };                           // end of EUTelRunHeaderImpl
}                               // eutelescope namespace

#endif // EUTELRUNHEADERIMPL
