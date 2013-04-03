/* 
 * File:   EUTelExhaustiveTrackFinder.h
 * Contact: denys.lontkovskyi@desy.de
 *
 * Created on January 27, 2013, 6:25 PM
 */

#ifndef EXHAUSTIVETRACKFINDER_H
#define	EXHAUSTIVETRACKFINDER_H

// eutelescope includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelTrackFinder.h"
//#include "EUTelUtility.h"

// lcio includes <.h>
#include "LCIOSTLTypes.h"
#include <IMPL/TrackerHitImpl.h>

// system includes <>
#include <map>

namespace eutelescope {

    class EUTelExhaustiveTrackFinder : public EUTelTrackFinder {
    public:
        EUTelExhaustiveTrackFinder() : 
                EUTelTrackFinder( "EUTelExhaustiveTrackFinder" ),
                _allowedMissingHits(0),
                _maxTrackCandidates(9999),
                _mode(2),
		_residualsYMin(),
		_residualsXMin(),
		_residualsXMax(),
		_residualsYMax(),
		_distanceMaxVec() {};

        EUTelExhaustiveTrackFinder(std::string name) : 
               EUTelTrackFinder(name),
	       _allowedMissingHits(0), 
      	       _maxTrackCandidates(9999), 
   	       _mode(2),
   	       _residualsYMin(),
	       _residualsXMin(),
	       _residualsXMax(),
	       _residualsYMax(),
	       _distanceMaxVec()
	{
            _residualsYMin.clear();
            _residualsXMax.clear();
            _residualsXMin.clear();
            _residualsYMax.clear();
        };

        EUTelExhaustiveTrackFinder(std::string name, int allowedMissingHits, int maxTrackCandidates) : 
		EUTelTrackFinder(name), 
		_allowedMissingHits(allowedMissingHits),
		_maxTrackCandidates(maxTrackCandidates),
                _mode(2),
	        _residualsYMin(),
	        _residualsXMin(),
	        _residualsXMax(),
	        _residualsYMax(),
	        _distanceMaxVec()
	{
            _residualsYMin.clear();
            _residualsXMax.clear();
            _residualsXMin.clear();
            _residualsYMax.clear();
        };

        virtual ~EUTelExhaustiveTrackFinder() {
        };

        inline int GetAllowedMissingHits() const {
            return _allowedMissingHits;
        };

        inline void SetAllowedMissingHits(int allowedMissingHits) {
            this->_allowedMissingHits = allowedMissingHits;
        }
        
        inline int GetMaxTrackCandidates() const {
            return _maxTrackCandidates;
        };

        inline void SetMaxTrackCandidates(int maxTrackCandidates) {
            this->_maxTrackCandidates = maxTrackCandidates;
        }
        
        void SetResidualsYMin( const EVENT::FloatVec& residualsYMin ) {
            this->_residualsYMin.resize(residualsYMin.size());
            std::copy(residualsYMin.begin(), residualsYMin.end(), _residualsYMin.begin());
        }
        
        void SetResidualsYMax( const EVENT::FloatVec& residualsYMax ) {
            this->_residualsYMax.resize(residualsYMax.size());
            std::copy(residualsYMax.begin(), residualsYMax.end(), _residualsYMax.begin());
        }
        
        void SetResidualsXMin( const EVENT::FloatVec& residualsXMin ) {
            this->_residualsXMin.resize(residualsXMin.size());
            std::copy(residualsXMin.begin(), residualsXMin.end(), _residualsXMin.begin());
        }
        
        void SetResidualsXMax( const EVENT::FloatVec& residualsXMax ) {
            this->_residualsXMax.resize(residualsXMax.size());
            std::copy(residualsXMax.begin(), residualsXMax.end(), _residualsXMax.begin());
        }

        void SetDistanceMaxVec( const EVENT::FloatVec& distanceMaxVec) {
            this->_distanceMaxVec.resize(distanceMaxVec.size());
            std::copy(distanceMaxVec.begin(), distanceMaxVec.end(), _distanceMaxVec.begin());
        }
        
        inline int GetMode() const { return _mode; }

        inline void SetMode(int _mode) { this->_mode = _mode; }
        
    protected:
        EUTelTrackFinder::SearchResult DoTrackSearch();
        
    private:
        void FindTracks1(int&, std::map< int, EVENT::TrackerHitVec >&, EVENT::TrackerHitVec&, const std::vector< EVENT::TrackerHitVec> &, unsigned int, EVENT::TrackerHit*);
        void FindTracks2(int&, std::map< int, EVENT::TrackerHitVec >&, EVENT::TrackerHitVec&, const std::vector< EVENT::TrackerHitVec> &, unsigned int, EVENT::TrackerHit*);
        
        void PruneTrackCandidates( std::map< int, EVENT::TrackerHitVec >& );

    private:
        int _allowedMissingHits;
        int _maxTrackCandidates;
        int _mode;
        
        EVENT::FloatVec _residualsYMin;
        EVENT::FloatVec _residualsXMin;
        EVENT::FloatVec _residualsXMax;
        EVENT::FloatVec _residualsYMax;
        
        EVENT::FloatVec _distanceMaxVec;
    };

}
#endif	/* EXHAUSTIVETRACKFINDER_H */

