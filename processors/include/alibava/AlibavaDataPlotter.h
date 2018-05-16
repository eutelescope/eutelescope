/*
 * Created by Thomas Eichhorn
 *  (2018 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

#ifndef ALIBAVADATAPLOTTER_H
#define ALIBAVADATAPLOTTER_H 1

// alibava includes ".h"
#include "AlibavaBaseProcessor.h"

// marlin includes ".h"
#include "marlin/Processor.h"

// lcio includes <.h>
#include <IMPL/LCRunHeaderImpl.h>
#include <IMPL/TrackerDataImpl.h>

// ROOT includes <>
#include "TObject.h"

// system includes <>
#include <string>
#include <list>

namespace alibava
{

    class AlibavaDataPlotter : public alibava::AlibavaBaseProcessor
    {
	public:

	    virtual Processor * newProcessor ( )
	    {
		return new AlibavaDataPlotter;
	    }

	    AlibavaDataPlotter ( );

	    virtual void init ( );

	    virtual void processRunHeader ( LCRunHeader * run );

	    virtual void processEvent ( LCEvent * evt );

	    virtual void check ( LCEvent * evt );

	    virtual void bookHistos ( );

	    virtual void bookEventHisto ( int eventnum );

	    virtual void fillHistos ( TrackerDataImpl * trkdata );

	    virtual void end ( );

	    void fillEventHisto ( int eventnum, TrackerDataImpl * trkdata );

	    void fillOtherHistos ( TrackerDataImpl * trkdata, float tdctime, float temperature );

	    bool isEventToBePlotted ( int eventnum );

	    std::string getEventHistoName ( int eventnum,  int ichip );

	    std::string _inputCollectionName;

	    IntVec _plotEvents;

	    float _multiplySignalby;

	    float _plotXPercentOfEvents;

	protected:

	};

	//! A global instance of the processor
	AlibavaDataPlotter gAlibavaDataPlotter;
}

#endif
