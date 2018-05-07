/*
 * Created by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

#ifndef AlibavaSimConverter_H
#define AlibavaSimConverter_H 1

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

    class AlibavaSimConverter:public alibava::AlibavaBaseProcessor
    {

	public:

	    virtual Processor * newProcessor ( )
	    {
		return new AlibavaSimConverter;
	    }

	    AlibavaSimConverter ( );

	    virtual void init ( );

	    virtual void processRunHeader ( LCRunHeader * run );

	    virtual void processEvent ( LCEvent * evt );

	    virtual void check ( LCEvent * evt );

	    virtual void end ( );

	    std::string _outputcollectionname;

	    std::string _nonsensitiveaxis;

	    double _commonmodemean;

	    double _commonmodesigma;

	    double _noisemean;

	    double _noisesigma;

	    double _pedestaldb[256] = { 0.0 };

	    double _pedestalmean;

	    double _pedestalsigma;

	    double _scalefactor;

	    double _tdcmean;

	    double _tdcsigma;

	    EVENT::IntVec _chipSelection;

	protected:

	private:

	    void checkIfChipSelectionIsValid ( );

	};

    AlibavaSimConverter gAlibavaSimConverter;

}

#endif
