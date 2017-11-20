/*
 * Created by Thomas Eichhorn
 *  (2017 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

#ifndef CBCCLUSTERING_H
#define CBCCLUSTERING_H 1

// marlin includes ".h"
#include "marlin/Processor.h"

// system includes <>
#include <string>

namespace eutelescope
{

    //! CMS CBC clustering processor for Marlin.

    class CBCClustering : public marlin::Processor
    {

	public:

	    virtual Processor * newProcessor ( )
	    {
		return new CBCClustering;
	    }

	    //! Default constructor
	    CBCClustering ();

	    virtual void init ();

	    virtual void processRunHeader (LCRunHeader * run);

	    virtual void processEvent (LCEvent * evt);

	    virtual void check (LCEvent * evt);

	    void bookHistos();

	    void fillHistos();

	    virtual void end();

	    std::string _cbcInputCollectionName;

	    std::string _cbcDataOutputCollectionName;

	    std::string _cbcPulseOutputCollectionName;

	    std::string _nonsensitiveaxis;

	    int _chancount;

	    int _maxclusters;

	    int _maxclustersize;

	    int _outputSensorID;

	    int _zsmode;

	protected:

	    std::map < std::string, AIDA::IBaseHistogram * > _aidaHistoMap;

    };

    //! A global instance of the processor
    CBCClustering gCBCClustering;

}

#endif
