/** 
 * File:   EUTelExhaustiveTrackFinder.h
 * Contact: denys.lontkovskyi@desy.de
 *
 */

#ifndef EXHAUSTIVETRACKFINDER_H
#define	EXHAUSTIVETRACKFINDER_H 1

// eutelescope includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelTrackFinder.h"
#include "EUTelUtility.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

// lcio includes <.h>
#include "LCIOSTLTypes.h"
#include <IMPL/TrackerHitImpl.h>

// system includes <>
#include <map>
#include <cmath>

namespace eutelescope {

    /** @class EUTelExhaustiveTrackFinder perfors search of trach candidates 
     *  Search algorithm is based on enumeration of all possible combinations
     * of hits that may come from a single track
     * 
     */
    class EUTelExhaustiveTrackFinder : public EUTelTrackFinder {
    private:
        DISALLOW_COPY_AND_ASSIGN(EUTelExhaustiveTrackFinder)    // prevent users from making (default) copies of processors
        
    public:
        EUTelExhaustiveTrackFinder() : 
                EUTelTrackFinder( "EUTelExhaustiveTrackFinder" ),
                _allowedMissingHits(0),
                _maxTrackCandidates(100),
                _mode(2),
                _nEmptyPlanes(0),
		_residualsYMin(),
		_residualsXMin(),
		_residualsXMax(),
		_residualsYMax(),
		_distanceMaxVec() {};

        EUTelExhaustiveTrackFinder(std::string name) : 
               EUTelTrackFinder(name),
	       _allowedMissingHits(0), 
      	       _maxTrackCandidates(100), 
   	       _mode(2),
               _nEmptyPlanes(0),
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

        EUTelExhaustiveTrackFinder(std::string name, unsigned int allowedMissingHits, unsigned int maxTrackCandidates) : 
		EUTelTrackFinder(name), 
		_allowedMissingHits(allowedMissingHits),
		_maxTrackCandidates(maxTrackCandidates),
                _mode(2),
                _nEmptyPlanes(0),
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

        // Getters and Setters
        
    public:
        inline int GetAllowedMissingHits() const {
            return _allowedMissingHits;
        };

        inline void SetAllowedMissingHits( unsigned int allowedMissingHits) {
            this->_allowedMissingHits = allowedMissingHits;
        }
        
        
    public:
        inline int GetMaxTrackCandidates() const {
            return _maxTrackCandidates;
        };

        inline void SetMaxTrackCandidates( unsigned int maxTrackCandidates) {
            this->_maxTrackCandidates = maxTrackCandidates;
        }
        
        
    public:
        void SetResidualsYMin( const EVENT::FloatVec& residualsYMin ) {
            this->_residualsYMin.resize(residualsYMin.size());
            std::copy(residualsYMin.begin(), residualsYMin.end(), _residualsYMin.begin());
        }
        
        void SetResidualsYMax( const EVENT::FloatVec& residualsYMax ) {
            this->_residualsYMax.resize(residualsYMax.size());
            std::copy(residualsYMax.begin(), residualsYMax.end(), _residualsYMax.begin());
        }
        
        
    public:
        void SetResidualsXMin( const EVENT::FloatVec& residualsXMin ) {
            this->_residualsXMin.resize(residualsXMin.size());
            std::copy(residualsXMin.begin(), residualsXMin.end(), _residualsXMin.begin());
        }
        
        void SetResidualsXMax( const EVENT::FloatVec& residualsXMax ) {
            this->_residualsXMax.resize(residualsXMax.size());
            std::copy(residualsXMax.begin(), residualsXMax.end(), _residualsXMax.begin());
        }

    public:
        void SetDistanceMaxVec( const EVENT::FloatVec& distanceMaxVec) {
            this->_distanceMaxVec.resize(distanceMaxVec.size());
            std::copy(distanceMaxVec.begin(), distanceMaxVec.end(), _distanceMaxVec.begin());
        }
        
    public:
        inline int GetMode() const { return _mode; }

        inline void SetMode(int mode = 2) { this->_mode = mode; }
        
    public:
        inline int GetNEmptyPlanes() const {
            return _nEmptyPlanes;
        };

        inline void SetNEmptyPlanes( int nEmptyPlanes ) {
            this->_nEmptyPlanes = nEmptyPlanes;
        }
        
    protected:
        EUTelTrackFinder::SearchResult DoTrackSearch();
        
    private:
        void FindTracks( int, std::vector< EVENT::TrackerHitVec >&, std::vector< EVENT::TrackerHitVec>& );
        void FindTracks1(int&, std::vector< EVENT::TrackerHitVec >&, EVENT::TrackerHitVec&, std::map< int, EVENT::TrackerHitVec> &, int, EVENT::TrackerHit*);
        void FindTracks2(int&, std::vector< EVENT::TrackerHitVec >&, EVENT::TrackerHitVec&, std::map< int, EVENT::TrackerHitVec> &, int, EVENT::TrackerHit*);
        
        void PruneTrackCandidates( std::vector< EVENT::TrackerHitVec >& );
        
        bool IsGoodCandidate( const EVENT::TrackerHitVec& );

    private:
        int _allowedMissingHits;
        int _maxTrackCandidates;
        int _mode;
        
        int _nEmptyPlanes;
        
        EVENT::FloatVec _residualsYMin;
        EVENT::FloatVec _residualsXMin;
        EVENT::FloatVec _residualsXMax;
        EVENT::FloatVec _residualsYMax;
        
        EVENT::FloatVec _distanceMaxVec;
        
        /** @struct CombinationGenerator 
         * This structure allows to generate all possible combinations
         * of given objects.
         * 
         * Generator is defined by a set of positions that are filled 
         * with values from a given alphabet.
         * 
         * Example: suppose one has a tuple (.,.,.) with 3 units. Possible 
         * values of each unit are the following
         * 1 11 123
         * 2 22 456
         * 3 33 789
         * 
         * This structure will successively generate all possible conbinations of
         * those, like (1,11,123), (2,33,456) etc.
         * 
         * The stucture has to be initialized with a map of vectors of possible
         * objects. Afterwards one can iterate over possible combinations using
         * incrementCurrentCombination()
         * 
         * \b Usage:
         * \code{.cpp}
         * vector< int > a; a.push_back(1); a.push_back(2); a.push_back(4);
         * vector< int > b; b.push_back(11); b.push_back(22); b.push_back(33);
         * vector< int > c; c.push_back(123); c.push_back(456); c.push_back(789);
         * map< int, vector< int > > mapping;
         * mapping.insert( make_pair( 1, a ) );
         * mapping.insert( make_pair( 2, b ) );
         * mapping.insert( make_pair( 3, c ) );
         * 
         * vector< int > comb;
         * CombinationGenerator< int > st(mapping);
         * do {
         *       st.printCombination();
         *       comb = st.getCurrentCombination();
         *    } while( st.incrementCurrentCombination( ) );
         * \endcode
         */
        struct CombinationGenerator {
            
            DISALLOW_COPY_AND_ASSIGN( CombinationGenerator )     // prevent users from making (default) copies
            
            // Constructors
            
            CombinationGenerator() : _trackFinder(0), _wordAlphabet(), _currentCombination(), _endsState(), _beginsState() {
                initCurrentCombination();
            }

            CombinationGenerator( const std::vector< std::vector< EVENT::TrackerHit* > >& vec ) : _trackFinder(0), _wordAlphabet(vec), _currentCombination(),
                                                                                              _endsState(), _beginsState() {
                initCurrentCombination();
            }
            
            CombinationGenerator( EUTelExhaustiveTrackFinder* finder, const std::vector< std::vector< EVENT::TrackerHit* > >& vec ) : _trackFinder(finder),
            _wordAlphabet(vec), _currentCombination(), _endsState(), _beginsState() {
                initCurrentCombination();
            }

            /** Print current combination */
            void printCombination() {
                std::vector< std::vector< EVENT::TrackerHit* >::iterator >::iterator itr;
                for ( itr = _currentCombination.begin(); itr != _currentCombination.end(); ++itr ) std::cout << *(*itr) << " ";
                std::cout << std::endl;
            }

            /** Obtain current combination */
            std::vector < EVENT::TrackerHit* > getCurrentCombination() {
                std::vector < EVENT::TrackerHit* > result;
                std::vector< std::vector< EVENT::TrackerHit* >::iterator >::iterator itr;
                for (itr = _currentCombination.begin(); itr != _currentCombination.end(); ++itr) {
                    result.push_back( *(*itr) );
                }
                return result;
            }
            
            /** Iterate current combination */
            bool incrementCurrentCombination() {

                bool isMoreCombiantions = true;

                // check if this is the last possible combination. If yes, restart.
                {
                    bool isEnd = true;
                    std::vector< std::vector< EVENT::TrackerHit* > >::iterator itrMap = _wordAlphabet.begin();
                    std::vector< std::vector<EVENT::TrackerHit*>::iterator >::iterator itr;
                    for (itr = _currentCombination.begin(); itr != _currentCombination.end(); ++itr) {
                        if (*itr != (itrMap->end() - 1)) isEnd = false;
                        ++itrMap;
                    }

                    if (isEnd) {
                        initCurrentCombination(); // Restart if this is the last combination
                        isMoreCombiantions = false;
                        return isMoreCombiantions;
                    }
                }

                // apply constraints
		bool def_begin = false;         // flag: increment pointer of the hit in first plane or not
                // run in reverese direction from last plane to the second from the beginnig
		for ( int itrRev = _currentCombination.size()-1; itrRev >= 2; --itrRev ) {
                        // check if two hits in adjacent planes are in small window
                        if ( IsBadCandidate(*_currentCombination[itrRev], *_currentCombination[itrRev-1]) )
			{
				// if _currentCombination[itrRev-1] points to the end of it's possible hits
				if ( distance( _currentCombination[itrRev-1], _endsState[itrRev-1] ) <= 1 )  {
                                        int itrTemp = -1;
                                        // find first plane to the right with not the last hit
					for ( size_t itrX = itrRev; itrX < _currentCombination.size(); ++itrX )
                                                if ( _currentCombination[itrX] != _endsState[itrX] -1 ) { 
                                                    ++_currentCombination[itrX];
                                                    itrTemp = itrX;
                                                    break; 
                                                }
                                        if ( itrTemp == -1 ) {  // if such a hit was not found, stop. All combinations were enumerated
                                            isMoreCombiantions = false;
                                            return isMoreCombiantions;
                                        }
                                        // reset everything to  the left
					for ( int itrZero = itrRev-1; itrZero >= 1; --itrZero ) {
						_currentCombination[itrZero] = _beginsState[itrZero];
					}
                                        // set first hit in the first plane
					if ( itrRev-1 != 0 ) { _currentCombination[0] = _beginsState[0]; def_begin = true; }
                                        // move to the correct plane for the next round of checks
                                        itrRev = ( itrTemp == static_cast<int>(_currentCombination.size()-1)) ? itrTemp + 1: itrTemp + 2;
				} else {
					++_currentCombination[itrRev-1];        // increment pointer of the hit on the left
                                        // reset everything to  the left
					for ( int itrZero = itrRev-2; itrZero >= 1; --itrZero ) {
						_currentCombination[itrZero] = _beginsState[itrZero];
					}
                                        // set first hit in the first plane
					if ( itrRev-1 != 0 ) { _currentCombination[0] = _beginsState[0]; def_begin = true; }
                                        // move to the correct plane for the next round of checks
                                        itrRev = itrRev + 1;
				}
			}
		}
                
                // increment current combination
                {
			bool incrementNext = false;
			for ( size_t itr = 0; itr < _currentCombination.size(); ++itr ) {
				if ( !def_begin && itr == 0 ) _currentCombination[itr] = ++(_currentCombination[itr]);
				if ( incrementNext ) { _currentCombination[itr] = ++(_currentCombination[itr]); incrementNext = false; }
				if ( _currentCombination[itr] == _endsState[itr] ) { incrementNext = true; _currentCombination[itr] = _beginsState[itr]; }
			}
		}
                return isMoreCombiantions;
            }
            
        private:

            /** Initialise starting combination */
            void initCurrentCombination() {
                _currentCombination.resize(0);
                std::vector < std::vector< EVENT::TrackerHit* > >::iterator itr;
                for ( itr = _wordAlphabet.begin(); itr != _wordAlphabet.end(); ++itr ) {
                    _currentCombination.push_back(itr->begin());
                    _beginsState.push_back( itr->begin() );
                    _endsState.push_back( itr->end() );
                }
            }

           /** Check if track candidate satisfies selection requirements.
            *  Track candidate passes if all of it hits are in a window
            *  defined by _residualsXMin, _residualsXMax and _distanceMaxVec etc.
            *  Window size changes with distance between planes. By default 150mm
            *  is assumed.
            * 
            * @param hit1 pointer to the hit to the right
            * @param hit2 pointer to the hit to the left
            * @return true if two hits satisfy the requirements
            */
            bool IsBadCandidate( EVENT::TrackerHit* hit1, EVENT::TrackerHit* hit2 ) {
                bool isBad = true;

                const double zSpacing = 150.;   // [mm]
                const double* posHit1     = hit1->getPosition();
                const double* posHit2     = hit2->getPosition();
                const double resX = posHit1[ 0 ] - posHit2[ 0 ];
                const double resY = posHit1[ 1 ] - posHit2[ 1 ];
                const double resZ = fabs(posHit1[ 2 ] - posHit2[ 2 ]);
                const double resR = resX*resX + resY*resY;
                const int sensorID = Utility::GuessSensorID( static_cast< IMPL::TrackerHitImpl* >(hit2) );
                const int numberAlongZ = geo::gGeometry().sensorIDtoZOrder( sensorID );
                if( _trackFinder->_mode == 1 ) {
                    if( resX > _trackFinder->_residualsXMin[ numberAlongZ ] * resZ / zSpacing ) {
                        isBad = false;
                    }
                    if( resX < _trackFinder->_residualsXMax[ numberAlongZ ] * resZ / zSpacing ) {
                        isBad = false;
                    }
                    if( resY > _trackFinder->_residualsYMin[ numberAlongZ ] * resZ / zSpacing ) {
                        isBad = false;
                    }
                    if( resY < _trackFinder->_residualsYMax[ numberAlongZ ] * resZ / zSpacing ) {
                        isBad = false;
                    }
                } else {
                    if ( sqrt( resR ) < _trackFinder->_distanceMaxVec [ numberAlongZ ] * resZ / zSpacing  ) {
                        isBad = false;
                    }
                }
            
            return isBad;
        }
            
        private:
            /** Instance of enclosing class */
            EUTelExhaustiveTrackFinder* _trackFinder;
            
            /** Vector of possible values on each position */
            std::vector< std::vector< EVENT::TrackerHit* > > _wordAlphabet;
            
            /** Current combination */
            std::vector< std::vector< EVENT::TrackerHit* >::iterator > _currentCombination;
            std::vector< std::vector< EVENT::TrackerHit* >::iterator > _endsState;
            std::vector< std::vector< EVENT::TrackerHit* >::iterator > _beginsState;
        };
        
        friend struct CombinationGenerator;
        
    };

}
#endif	/* EXHAUSTIVETRACKFINDER_H */

