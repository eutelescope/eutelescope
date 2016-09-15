// based on EUTelDafFitter.h

#ifndef EUTELDAFFITTERANALYSIS_H
#define EUTELDAFFITTERANALYSIS_H

// marlin includes ".h"
#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>

// AIDA includes <.h>
#include <AIDA/IHistogram1D.h>

// system includes <>
#include <string>
#include <vector>
#include <map>

#include "EUTelDafBase.h"

namespace eutelescope {

	class EUTelDafFitterAnalysis : EUTelDafBase {

	public:

		// Marlin processor interface funtions
		//! Returns a new instance of EUTelDafFitter
		virtual Processor * newProcessor() { return new EUTelDafFitterAnalysis;}

		//! Default constructor
		EUTelDafFitterAnalysis ();

		//! Called at the job beginning.
		virtual void dafInit ();
		//! Called every event
		virtual void dafEvent (LCEvent* event);
		//! Called after data processing.
		virtual void dafEnd();
		//! Set fitter specific params
		virtual void dafParams();
		virtual double getZfromRefHit(int plane,int sensorID, double *pos);

		int findPairIndex(float x, float y, std::vector<double const*> vec);
		void fillDiffHistos(daffitter::TrackCandidate<float,4>& track, LCEvent* event);
		void bookDiffHistos();
 
	protected:

		//! Output track collection name
		std::string _trackCollectionName;
		std::string _trueHitCollectionName;
		//! Output track collection
		LCCollectionVec* _fittrackvec;
		LCCollectionVec* _fitpointvec;
		std::map<int, std::vector<double const*>> _trueHitMap;
		//LCCollectionVec* _trueHitCollectionVec;
		void addToLCIO(daffitter::TrackCandidate<float,4>& track, LCCollectionVec *lcvec);
		//! LCIO switch
		bool _addToLCIO, _fitDuts;
	};

	//! A global instance of the processor
	EUTelDafFitterAnalysis gEUTelDafFitterAnalysis;
}
#endif
