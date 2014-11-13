/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 */


#ifndef ALIBAVACONVERTER_H
#define ALIBAVACONVERTER_H 1

// personal includes ".h"
#include "ALIBAVA.h"
#include "AlibavaRunHeaderImpl.h"

// marlin includes ".h"
#include "marlin/DataSourceProcessor.h"

// lcio includes <.h>

// system includes <>
#include <string>
#include <vector>
#include <ctime>

namespace alibava {


class AlibavaConverter : public marlin::DataSourceProcessor    {

  public:

    //! Default constructor
    AlibavaConverter ();

    //! New processor
    /*! Return a new instance of a AlibavaConverter. It is
     *  called by the Marlin execution framework and shouldn't be used
     *  by the final user.
     */
    virtual AlibavaConverter * newProcessor ();

    //! Creates events from the eudaq software
    virtual void readDataSource (int numEvents);

    //! Init method
    /*! It is called at the beginning of the cycle and it prints out
     *  the parameters.
     */
    virtual void init ();

    //! End method
    /*! It prints out a good bye message
     */
    virtual void end ();

	double tdc_time(unsigned int atime);
	double get_temperature(unsigned short temp);


  protected:

    //! The input file name
    /*! It is set as a Marlin parameter in the constructor
     */
    std::string _fileName;


    //! Geometry ID
    /*! This is the unique identification number of the telescope
     *  geometry. This identification number is saved in the run
     *  header and then crosscheck against the XML geometry
     *  description during the reconstruction phase.
     *
     *  In the future, this ID can be used to browse a geometry database.
     */
    int _geoID;

	//! The run number
	/*! It is set as a Marlin parameter in the constructor
	 */
	std::string _formattedRunNumber;

	//! The run number
	/*! The _formattedRunNumber changed to integer and saved as _runNumber
	 */
	int _runNumber;

	//! The name of the collection
	/*! The _rawDataCollectionName will be the base of the collection name that 
	 *  raw data will be saved.
	 *  For the collections created for each chip separately chip number information
	 *  will be appended to this name 
	 *
	 *  Example: if ALIBAVA::NOOFCHIPS=2 and _rawDataCollectionName="rawdata"
	 *  these collections will be created:
	 *  rawdata            // containing all chips
	 *  rawdata_chip0      // containing only first chip
	 *  rawdata_chip1      // containing only second chip
	 */
	std::string _rawDataCollectionName;
	
	/////////////////////////
	// optional parameters //
	/////////////////////////

	//! The selected chip numbers
	/*! _chipSelection vector stores the chip numbers that are selected
	 *  if not set, chip 0 and 1 will be stored
	 */
	EVENT::IntVec _chipSelection;
	
	//! The tilt angle
	/*! The tilt angle of the detector can be stored here
	 */
	float _tiltAngle;

	//! The sensor temperature
	/*! The temperature of the sensor can be stored here
	 */
	float _sensorTemperature;

	//! The start event number
	/*! The event number that AlibavaConverter should start storing.
	 *  Default value is -1, in this case it will store every event
	 */
	int _startEventNum;

	//! The stop event number
	/*! The event number that AlibavaConverter should stop storing.
	 *  Default value is -1, in this case it will store every event
	 */
	int _stopEventNum;

	//! The read channels reverse option
	/*! Alibava read channels from right to left
	 *  if you want to revert this i.e. make it from left to right
	 *  set this parameter to true
	 */
	bool _readChannelsReverse;
	
	//! An option to store pedestal and noise values stored in header of alibava data file
	bool _storeHeaderPedestalNoise;

	
  private:
	//! To check if the chip selection is valid
	void checkIfChipSelectionIsValid();
	
  };

  //! A global instance of the processor
  AlibavaConverter gAlibavaConverter;

} // end of alibava namespace
#endif

