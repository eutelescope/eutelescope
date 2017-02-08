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
#include "mille/Mille.h"
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
#include "TMinuit.h"
#include <TSystem.h>
#include <TMath.h>
#include "TVector3.h"
#include <TProfile.h>
class TMinuit;
#endif


namespace eutelescope
{

	class EUTelMilleGBL : public marlin::Processor
	{

		public:

		#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)

		std::map<std::string , TObject * > _rootObjectMap;

		class hit
		{

			public:

			hit() {}
			hit(double tx, double ty, double tz, double rx, double ry, double rz, int i)
			{
				x = tx;
				y = ty;
				z = tz;
				resolution_x = rx;
				resolution_y = ry;
				resolution_z = rz;
				planenumber = i;
			}

			double x;
			double y;
			double z;
			double resolution_x;
			double resolution_y;
			double resolution_z;
			int planenumber;
		};

		#endif

		class HitsInPlane
		{

			public:

			HitsInPlane()
			{
				measuredX = 0.0;
				measuredY = 0.0;
				measuredZ = 0.0;
			}
			HitsInPlane(double x, double y, double z)
			{
				measuredX = x;
				measuredY = y;
				measuredZ = z;
			}
			bool operator<(const HitsInPlane& b) const
			{
				return (measuredZ < b.measuredZ);
			}
			double measuredX;
			double measuredY;
			double measuredZ;
		};

		virtual Processor * newProcessor()
		{
			return new EUTelMilleGBL;
		}

		EUTelMilleGBL ();

		virtual void init ();

		virtual void processRunHeader (LCRunHeader * run);

		virtual void processEvent (LCEvent * evt);

		virtual void end();

		virtual bool isPointInPlane(TVector3 point);

		void bookHistos();

		bool _BField;

		protected:

		//! Ordered sensor ID
		/*! Within the processor all the loops are done up to _nPlanes and
		 *  according to their position along the Z axis (beam axis).
		 *
		 *  This vector is containing the sensorID sorted according to the
		 *  same rule.
		*/

		std::vector< int > _orderedSensorID;
		std::vector< int > _orderedSensorID_wo_excluded;

		std::vector<std::string > _hitCollectionName;

		//only for internal usage
		std::vector<int > _excludePlanes;

		//this is going to be set by the user.
		std::vector<int > _excludePlanes_sensorIDs;

		//only for internal usage
		std::vector<int > _FixedPlanes;
		std::vector<int > _FixedPlanes_sensorIDs;

		//this is going to be
		//set by the user.
		double _eBeam;
		double _triCut;
		double _triCutDUTx;
		double _triCutDUTy;
		double _driCut;
		double _driCutREFx;
		double _driCutREFy;
		double _sixCut;
		double _resx;
		double _resy;
		double _resdutx;
		double _resduty;
		double _resrefx;
		double _resrefy;

		//! Id of telescope layer which should be used as a DUT
		int _manualDUTid;

		int _dutFitMethod;

		int _doPreAlignment;

		int _doDUTAlignment;

		IntVec _FixParameter;

		double _probCut;
		double _chi2ndfCut;
		double _slopecutDUTx;
		double _slopecutDUTy;
		int _probCutCount;
		int _chi2ndfCutCount;
		int _cutfailsx;
		int _cutfailsy;

		int _maxTrackCandidatesTotal;

		std::string _binaryFilename;

		int _alignMode;

		int _generatePedeSteerfile;
		std::string _pedeSteerfileName;
		int _runPede;

		int _inputMode;

		int _requireDUTHit;

		bool _useREF;

		bool _dthistos;

		float _testModeSensorResolution;
		float _testModeXTrackSlope;
		float _testModeYTrackSlope;
		std::vector<float > _testModeSensorZPositions;
		std::vector<float > _testModeSensorXShifts;
		std::vector<float > _testModeSensorYShifts;
		std::vector<float > _testModeSensorGamma;
		std::vector<float > _testModeSensorAlpha;
		std::vector<float > _testModeSensorBeta;

		std::string _alignmentConstantLCIOFile;
		std::string _alignmentConstantCollectionName;
		std::string _coordinatorPreAlignmentCollectionName;

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
		double _coordinator_x;
		double _coordinator_y;
		double _coordinator_z;

		bool _alignmentloaded;
		bool _pre_alignmentloaded;
		bool _dutplanevectors;

		// define axis for rotations:
		TVector3 _Xaxis;
		TVector3 _Yaxis;
		TVector3 _Zaxis;

		//! Alignment collection name
		std::vector<std::string> _alignmentCollectionName;

		//! Prealignment collection name
		std::vector<std::string> _pre_alignmentCollectionName;

		std::string _outputHitCollectionName;

		std::string _outputTrackCollectionName;

		private:

		TVector3 _dutplane;
		TVector3 _normalvector;

		//! Run number
		int _iRun;

		//! Event number
		int _iEvt;

		//! Excluded planes
		int _nExcludePlanes;

		//! Statistics
		int _nMilleDataPoints;
		int _nMilleTracks;
		int _nGoodMilleTracks;

		//! The function to get the alignment from file
		void getAlignment (LCEvent * event);

		//! The function to get the prealignment from file
		void getPreAlignment (LCEvent * event);

		void getCoordinatorAlignment (LCEvent * event);

		void getDUTNormal();

		//! Silicon planes parameters as described in GEAR
		/*! This structure actually contains the following:
		 *  @li A reference to the telescope geoemtry and layout
		 *  @li An integer number saying if the telescope is w/ or w/o DUT
		 *  @li An integer number saying the number of planes in the
		 *  telescope.
		 *
		 *  This object is provided by GEAR during the init() phase and
		 *  stored here for local use.
		*/
		gear::SiPlanesParameters * _siPlanesParameters;

		//! Silicon plane layer layout
		/*! This is the real geoemetry description. For each layer
		 *  composing the telescope the relevant information are
		 *  available.
		 *
		 *  This object is taken from the _siPlanesParameters during the
		 *  init() phase and stored for local use
		 */
		gear::SiPlanesLayerLayout * _siPlanesLayerLayout;

		#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
		//! AIDA histogram map
		/*! Instead of putting several pointers to AIDA histograms as
		 *  class members, histograms are booked in the init() method and
		 *  their pointers are inserted into this map keyed by their
		 *  names.
		 *  The histogram filling can proceed recalling an object through
		 *  its name
		 */
		std::map<std::string, AIDA::IBaseHistogram * > _aidaHistoMap;
		#endif

		int _nPlanes;

		std::vector<std::vector<double> > _xPos;
		std::vector<std::vector<double> > _yPos;
		std::vector<std::vector<double> > _zPos;

		std::vector<double> _siPlaneZPosition;

	};

	//! A global instance of the processor
	EUTelMilleGBL gEUTelMilleGBL;

}
#endif
#endif
