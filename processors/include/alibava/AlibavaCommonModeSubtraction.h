/*
 * Created by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@desy.de
 */

#ifndef ALIBAVACOMMONMODESUBTRACTION_H
#define ALIBAVACOMMONMODESUBTRACTION_H 1

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
	class AlibavaCommonModeSubtraction : public alibava::AlibavaBaseProcessor
	{
	    public:
		virtual Processor * newProcessor ( )
		{
		    return new AlibavaCommonModeSubtraction;
		}

		AlibavaCommonModeSubtraction ( );

		virtual void init ( );

		virtual void processRunHeader ( LCRunHeader * run );

		virtual void processEvent ( LCEvent * evt );

		virtual void check ( LCEvent * evt );

		void bookHistos ( );

		void fillHistos ( TrackerDataImpl * trkdata );

		virtual void end ( );

		std::string _commonmodeCollectionName;

		std::string _commonmodeerrorCollectionName;

	    protected:

		std::string _chanDataHistoName;

		std::string getChanDataHistoName ( int chipnum, int ichan );

		std::string getSignalCorrectionName ( );

	};

	//! A global instance of the processor
	AlibavaCommonModeSubtraction gAlibavaCommonModeSubtraction;
}

#endif
