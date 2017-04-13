/*
 * Created by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@cern.ch
 */

#ifndef AlibavaHeader_H
#define AlibavaHeader_H 1

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


namespace alibava {

	//! Example Alibava processor for Marlin.

	class AlibavaHeader:public alibava::AlibavaBaseProcessor   {

	public:

		virtual Processor * newProcessor () { return new AlibavaHeader; }

		AlibavaHeader ();

		virtual void init ();

		virtual void processRunHeader (LCRunHeader * run);

		virtual void processEvent (LCEvent * evt);

		virtual void check (LCEvent * evt);

		void bookHistos();

		void fillHistos(TrackerDataImpl * trkdata,TrackerDataImpl * trkdata2);

		virtual void end();

		std::string getChanDataHistoName(unsigned int ichan);

		std::string getChanDataFitName(unsigned int ichan);

		void calculateHeaderNoise();

		void correlateLastHeader(TrackerDataImpl * trkdata, TrackerDataImpl * trkdata2);

		std::string _rawdatacollection;

		float _pedestal14;

		float _pedestal15;

		float _firstpedestal;

		std::string _filterFileName;

	protected:

	};

	AlibavaHeader gAlibavaHeader;

}

#endif
