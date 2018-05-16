/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELMULTILINEFIT_H
#define EUTELMULTILINEFIT_H

// built only if GEAR is available
#ifdef USE_GEAR
// eutelescope includes ".h"

// marlin includes ".h"
#include "marlin/Processor.h"

// marlin util includes
#include "Mille.h"
#include "include/MilleBinary.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>
#include <IMPL/TrackerHitImpl.h>

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#endif

// system includes <>
#include <string>
#include <vector>
#include <map>

#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
//#include "TMinuit.h"
//#include <TSystem.h>
//#include <TMath.h>
#include "TVector3.h"
class TMinuit;
#endif


namespace eutelescope
{

    class EUTelMilleGBL : public marlin::Processor
    {

	public:

	    #if defined ( USE_ROOT ) || defined ( MARLIN_USE_ROOT )

	    std::map < std::string , TObject * > _rootObjectMap;

	    #endif

	    class HitsInPlane
	    {

		public:

		    HitsInPlane ( )
		    {
			measuredX = 0.0;
			measuredY = 0.0;
			measuredZ = 0.0;
		    }
		    HitsInPlane ( double x, double y, double z )
		    {
			measuredX = x;
			measuredY = y;
			measuredZ = z;
		    }
		    bool operator < ( const HitsInPlane& b ) const
		    {
			return ( measuredZ < b.measuredZ );
		    }
		    double measuredX;
		    double measuredY;
		    double measuredZ;
	    };

	    virtual Processor * newProcessor ( )
	    {
		return new EUTelMilleGBL;
	    }

	    EUTelMilleGBL ();

	    virtual void init ();

	    virtual void processRunHeader ( LCRunHeader * run );

	    virtual void processEvent ( LCEvent * evt );

	    virtual void end ( );

	    virtual bool isPointInPlane ( TVector3 point );

	    void bookHistos ( );

	    bool _BField;

	protected:

	    bool _alignmentloaded;
	    bool _dthistos;
	    bool _dutplanevectors;
	    bool _pre_alignmentloaded;
	    bool _useTrackFit;
	    bool _x0histos;

	    double _chi2ndfCut;
	    double _coordinator_x;
	    double _coordinator_y;
	    double _coordinator_z;
	    double _driCut;
	    double _driCutREFx;
	    double _driCutREFy;
	    double _eBeam;
	    double _isolationCut;
	    double _probCut;
	    double _resdutx;
	    double _resduty;
	    double _resrefx;
	    double _resrefy;
	    double _resx;
	    double _resy;
	    double _sixCut;
	    double _slopecutDUTx;
	    double _slopecutDUTy;
	    double _triCut;
	    double _triCutDUTx;
	    double _triCutDUTy;

	    double * _dut_align_x;
	    double * _dut_align_y;
	    double * _dut_align_z;
	    double * _dut_align_a;
	    double * _dut_align_b;
	    double * _dut_align_c;
	    double * _dut_align_x_error;
	    double * _dut_align_y_error;
	    double * _dut_align_z_error;
	    double * _dut_align_a_error;
	    double * _dut_align_b_error;
	    double * _dut_align_c_error;
	    double * _dut_pre_align_x;
	    double * _dut_pre_align_y;
	    double * _dut_pre_align_z;
	    double * _dut_pre_align_a;
	    double * _dut_pre_align_b;
	    double * _dut_pre_align_c;

	    int _alignMode;
	    int _chi2ndfCutCount;
	    int _cutfailsx;
	    int _cutfailsy;
	    int _doDUTAlignment;
	    int _doPreAlignment;
	    int _dutFitMethod;
	    int _generatePedeSteerfile;
	    int _manualDUTid;
	    int _maxTrackCandidatesTotal;
	    int _probCutCount;
	    int _requireDUTHit;
	    int _runPede;
	    int _useREF;

	    IntVec _FixParameter;

	    std::string _alignmentConstantCollectionName;
	    std::string _alignmentConstantLCIOFile;
	    std::string _binaryFilename;
	    std::string _coordinatorPreAlignmentCollectionName;
	    std::string _outputHitCollectionName;
	    std::string _outputTrackCollectionName;
	    std::string _pedeSteerfileName;
	    std::string _trackFitCollectionName;

	    TVector3 _Xaxis;
	    TVector3 _Yaxis;
	    TVector3 _Zaxis;

	    std::vector < int > _excludePlanes;
	    std::vector < int > _excludePlanes_sensorIDs;
	    std::vector < int > _FixedPlanes;
	    std::vector < int > _FixedPlanes_sensorIDs;
	    std::vector < int > _TelescopePlanes_sensorIDs;
	    std::vector < int > _orderedSensorID;
	    std::vector < int > _activeSensorID;
	    std::vector < int > _activeSensorZ;

	    std::vector < std::string > _alignmentCollectionName;
	    std::vector < std::string > _hitCollectionName;
	    std::vector < std::string > _pre_alignmentCollectionName;

	private:

	    int _iEvt;
	    int _nExcludePlanes;
	    int _nGoodMilleTracks;
	    int _nMilleDataPoints;
	    int _nMilleTracks;
	    int _nPlanes;

	    gear::SiPlanesParameters * _siPlanesParameters;
	    gear::SiPlanesLayerLayout * _siPlanesLayerLayout;

	    #if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
	    std::map < std::string, AIDA::IBaseHistogram * > _aidaHistoMap;
	    #endif

	    TVector3 _dutplane;
	    TVector3 _normalvector;

	    std::vector < double > _siPlaneZPosition;
	    std::vector < std::vector < double > > _xPos;
	    std::vector < std::vector < double > > _yPos;
	    std::vector < std::vector < double > > _zPos;

	    void getAlignment ( LCEvent * event );
	    void getPreAlignment ( LCEvent * event );
	    void getCoordinatorAlignment ( LCEvent * event );
	    void getDUTNormal ( );

    };

	//! A global instance of the processor
	EUTelMilleGBL gEUTelMilleGBL;

}
#endif
#endif
