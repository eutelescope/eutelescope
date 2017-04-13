/*
 * Created by Thomas Eichhorn
 *  (2014 DESY)
 *
 *  email:thomas.eichhorn@cern.ch
 *
 */

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTelFixDUTAlignment_h
#define EUTelFixDUTAlignment_h 1

#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include "lcio.h"
#include "EUTelEventImpl.h"
#include "EUTelAlignmentConstant.h"

// system includes <>
#include <string>
#include <vector>
#include <map>

using namespace std;

namespace eutelescope
{

	class EUTelFixDUTAlignment : public marlin::Processor
	{

		public:

		//! Returns a new instance of this processor. It is called by Marlin execution framework and it shouldn't be called/used by the final user.  @return a new EUTelFixDUTAlignment
		virtual Processor*  newProcessor() { return new EUTelFixDUTAlignment ; }

		//! Default constructor
		EUTelFixDUTAlignment() ;

		//! Called at the job beginning. This is executed only once in the whole run.
		virtual void init() ;

		//! Called for every run with @param run as the LCRunHeader of the current run
		virtual void processRunHeader( LCRunHeader* run ) ;

		//! Called for every event in the file with  @param evt as the current LCEvent event
		virtual void processEvent( LCEvent * evt ) ;

		//! Check event method. This method is called by the Marlin execution framework as soon as the processEvent is over. It can be used to fill check plots.
		virtual void check( LCEvent * evt ) ;

		//! Book histograms. Histogram pointers are stored into _aidaHistoMap so that they can be recalled and filled from anywhere in the code.
		void bookHistos();

		//! Called after data processing for clean up. Used to release memory allocated in init() step
		virtual void end() ;

		//! Silicon plane parameters as described in GEAR. This object is provided by GEAR during the init() phase and stored here for local use.
		gear::SiPlanesParameters * _siPlanesParameters;

		//! This is the real geometry description for each layer. This object is taken from _siPlanesParameters during the init() phase and stored for local use
		gear::SiPlanesLayerLayout * _siPlanesLayerLayout;

		protected:

		int _mode;

		std::string _outputfilename;

		double _shiftx;

		double _shifty;

		double _shiftz;

		double _rota;

		double _rotb;

		double _rotc;

		int _manualDUTid;

		int _manualDUTposition;

		int _eventstowrite;

		std::string _dummyhitcollectionname;

		void _EulerRotation(double* _telPos, double* _gRotation) ;

	} ;

	//! A global instance of the processor.
	EUTelFixDUTAlignment aEUTelFixDUTAlignment ;

}

#endif



