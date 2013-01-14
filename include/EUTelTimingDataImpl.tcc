// Version: $Id$
// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELTIMINGDATAIMPL_TCC
#define EUTELTIMINGDATAIMPL_TCC

namespace eutelescope {

	template<class DetectorType> 
	EUTelTimingDataImpl<DetectorType>::EUTelTimingDataImpl(IMPL::TrackerDataImpl * data) {

		std::auto_ptr<DetectorType> plane ( new DetectorType );
		_nElement     = plane->getNoOfElements();
		_type         = plane->getDetectorType();
		_trackerData  = data;

	} 

	template<class DetectorType>
	unsigned int EUTelTimingDataImpl<DetectorType>::size() const {
		return _trackerData->getChargeValues().size() / _nElement;
	}







}


#endif

