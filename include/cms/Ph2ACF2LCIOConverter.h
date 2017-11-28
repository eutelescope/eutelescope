/*
 * Created by Thomas Eichhorn
 *  (2017 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

#ifndef PH2ACF2LCIOCONVERTER_H
#define PH2ACF2LCIOCONVERTER_H 1

// marlin includes ".h"
#include "marlin/DataSourceProcessor.h"

// system includes <>
#include <string>
#include <vector>
#include <ctime>

namespace eutelescope
{

    class Ph2ACF2LCIOConverter : public marlin::DataSourceProcessor
    {
	public:

	    Ph2ACF2LCIOConverter ();

	    virtual Ph2ACF2LCIOConverter * newProcessor ( );

	    virtual void readDataSource ( int numEvents );

	    virtual void init ( );

	    virtual void end ( );

	    virtual unsigned createMask ( unsigned a, unsigned b );

	protected:

	    int _runNumber;

	    int _maxRecordNumber;

	    std::string _dataformat;

	    std::string _fileName;

	    std::string _formattedRunNumber;

	    std::string _rawDataCollectionNameBottom;

	    std::string _rawDataCollectionNameTop;

	private:

	    int _nFE;

	    int _nChips;

    };

    Ph2ACF2LCIOConverter gPh2ACF2LCIOConverter;

}
#endif

