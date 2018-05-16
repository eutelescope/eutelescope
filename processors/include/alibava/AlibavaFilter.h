/*
 * Created by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

#ifndef AlibavaFilter_H
#define AlibavaFilter_H 1

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

    class AlibavaFilter : public alibava::AlibavaBaseProcessor
    {
	public:

	    virtual Processor * newProcessor ( )
	    {
		return new AlibavaFilter;
	    }

	    AlibavaFilter ( );

	    virtual void init ( );

	    virtual void processRunHeader ( LCRunHeader * run );

	    virtual void processEvent ( LCEvent * evt );

	    virtual void check ( LCEvent * evt );

	    void bookHistos ( );

	    void fillHistos ( TrackerDataImpl * trkdata );

	    virtual void end ( );

	    bool _readcoefficients;
	    bool _rghcorrection;
	    bool _simplemethod;

	    double _readcoefficient1;
	    double _readcoefficient2;

	    float _initcoefficient1;
	    float _initcoefficient2;
	    float _maxrghadc;
	    float _minrghnoise;
	    float _maxsignalfactor;
	    float _rghnegnoisecut;
	    float _seedcut;

	    int _dropsuspectcount;
	    int _maxsuspects;

	    std::string _filteredCollectionName;
	    std::string _filterFileName;

	protected:

    };

    AlibavaFilter gAlibavaFilter;
}

#endif
