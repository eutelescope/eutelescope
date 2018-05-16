/*
 * Created by Eda Yildirim
 *  (2014 DESY)
 *
 *  email:eda.yildirim@cern.ch
 *
 *  modified by: Thomas Eichhorn thomas.eichhorn@desy.de
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

namespace alibava
{

    class AlibavaConverter : public marlin::DataSourceProcessor
    {

	public:

	    AlibavaConverter ( );

	    virtual AlibavaConverter * newProcessor ( );

	    virtual void readDataSource ( int numEvents );

	    virtual void init ( );

	    virtual void end ( );

	    double tdc_time ( unsigned int atime );

	    double get_temperature ( unsigned short temp );

	protected:

	    // The input file name
	    std::string _fileName;

	    // Geometry ID
	    int _geoID;

	    // The run number
	    std::string _formattedRunNumber;

	    // The run number
	    int _runNumber;

	    // The name of the collection
	    std::string _rawDataCollectionName;

	    // The name of the collection
	    std::string _rawChipHeaderCollectionName;

	    // optional parameters

	    // The selected chip numbers
	    EVENT::IntVec _chipSelection;

	    // The start event number
	    int _startEventNum;

	    // The stop event number
	    int _stopEventNum;

	    // An option to store pedestal and noise values stored in header of alibava data file
	    bool _storeHeaderPedestalNoise;

	private:

	    // To check if the chip selection is valid
	    void checkIfChipSelectionIsValid();

    };

    // A global instance of the processor
    AlibavaConverter gAlibavaConverter;

}
#endif
