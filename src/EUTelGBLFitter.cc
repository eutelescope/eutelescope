#ifdef USE_GBL

// its own header 
#include "EUTelGBLFitter.h"

// eutelescope includes ".h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTELESCOPE.h"
#include "EUTelNav.h"

// marlin util includes
#include "mille/Mille.h"

// ROOT
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TVector3.h"
#else
#error *** You need ROOT to compile this code.  *** 
#endif

// GBL
#include "include/GblTrajectory.h"
#include "include/GblPoint.h"
#include "include/GblData.h"
#include "include/BorderedBandMatrix.h"
#include "include/MilleBinary.h"
#include "include/VMatrix.h"


// system includes <>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <iomanip>
#include <iterator>
#include <algorithm>

namespace eutelescope {

	EUTelGBLFitter::EUTelGBLFitter() :
	_beamQ(-1),
	_eBeam(4.),
	_mEstimatorType(),
	_mille(0),
	_parameterIdXShiftsMap(),
	_parameterIdYShiftsMap(),
	_parameterIdZShiftsMap(),
	_parameterIdXRotationsMap(),
	_parameterIdYRotationsMap(),
	_parameterIdZRotationsMap()
	{}

	EUTelGBLFitter::~EUTelGBLFitter() {
	}
	//THIS IS THE SETTERS. Not simple get function but the action of all these functions it in the end to set a member variable to something.	
	void EUTelGBLFitter::setMeasurementCov(EUTelState& state){
        std::vector<double> hitcov(4);
		int izPlane = state.getLocation();
		if( _parameterIdXResolutionVec.size() > 0 && _parameterIdYResolutionVec.size() > 0 ){
			hitcov.at(0) = _parameterIdXResolutionVec[izPlane];
			hitcov.at(3) = _parameterIdYResolutionVec[izPlane];

			hitcov.at(0) *= hitcov.at(0); 
			hitcov.at(3) *= hitcov.at(3);
		}
		else{//Default in case you have not specified a variance  
			throw(lcio::Exception("There is no measurement variance specified.")); 	
		}
		streamlog_out(DEBUG0) << "SET COVARIANCE: State: " << state.getLocation() << "  Covariance X/Y  : " << hitcov[0] << " " <<hitcov[3] <<std::endl; 
		state.getHit().setCov(hitcov);
	}
	//Note that we take the planes themselfs at scatters and also add scatterers to simulate the medium inbetween. 
	void EUTelGBLFitter::setScattererGBL(gbl::GblPoint& point, EUTelState & state ) {
		streamlog_out(DEBUG1) << " setScattererGBL ------------- BEGIN --------------  " << std::endl;
		TMatrixDSym precisionMatrix =  state.getScatteringVarianceInLocalFrame();
		streamlog_out(MESSAGE1) << "The precision matrix being used for the sensor  "<<state.getLocation()<<":" << std::endl;
		streamlog_message( DEBUG0, precisionMatrix.Print();, std::endl; );
		point.addScatterer(state.getKinks(), precisionMatrix);
		streamlog_out(DEBUG1) << "  setScattererGBL  ------------- END ----------------- " << std::endl;
	}
	/// This will add measurement information to the GBL point
	void EUTelGBLFitter::setMeasurementGBL(gbl::GblPoint& point,EUTelState & state ){
		streamlog_out(DEBUG1) << " setMeasurementGBL ------------- BEGIN --------------- " << std::endl;
		TVectorD meas(2);//Remember we need to pass the same 5 since gbl expects this due to jacobian
		meas.Zero();
		meas[0] = state.getHit().getPosition()[0] - state.getPosition()[0];
		meas[1] = state.getHit().getPosition()[1] - state.getPosition()[1];
		TVectorD measPrec(2); 
        setMeasurementCov(state);
        double cov[4];
        state.getHit().getCov(cov);
		measPrec[0] = 1. / cov[0];	// cov(x,x)
		measPrec[1] = 1. / cov[3];	// cov(y,y)
		streamlog_out(DEBUG4) << "This is what we add to the measured point:" << std::endl;
		streamlog_out(DEBUG4) << "Residuals and precision matrix for the hit:" << std::endl;
		streamlog_out(DEBUG4) << "X:" << std::setw(20) << meas[0] << std::setw(20) << measPrec[0] <<"," << std::endl;
		streamlog_out(DEBUG4) << "Y:" << std::setw(20) << meas[1] << std::setw(20)  <<"," << measPrec[1] << std::endl;
		streamlog_out(DEBUG4) << "This H matrix:" << std::endl;
		streamlog_message( DEBUG0, state.getProjectionMatrix().Print();, std::endl; );
		/// The gbl library creates 5 measurement vector and 5x5 propagation matrix automatically. If  
        /// Measurement is added in local frame. Defined as a plane in the GBL fit which contains only position information. 
        /// The projection matrix which is 5x5 is defined as the transform from the global to local frame. 
        /// The last zero is the minimum precision before this is set to 0. TO DO:Remove this magic number
		point.addMeasurement(state.getProjectionMatrix(), meas, measPrec, 0);
		streamlog_out(DEBUG1) << " setMeasurementsGBL ------------- END ----------------- " << std::endl;
	}
    ///This function is the link between Millepede and the creation of global variable to the GBL points. 
	void EUTelGBLFitter::getGloPar(std::vector< gbl::GblPoint >& pointList, EUTelTrack & track, std::map< unsigned int, unsigned int> & linkMeas ){
        streamlog_out(DEBUG1) << "State measurement numbers: " << linkMeas.size()  << std::endl;
        for(size_t j=0; j < track.getStates().size(); j++){
            EUTelState state = track.getStates().at(j);
            streamlog_out(DEBUG1) << "State location" << state.getLocation() << std::endl;
            std::map<unsigned int, unsigned int>::iterator itLink = linkMeas.find(state.getLocation());
            if(itLink != linkMeas.end()){ 
                streamlog_out(DEBUG1) << "Add global derivatives" << std::endl;
                ///Internally millepede will save a jacobian and the correct global labels for this setup
                _MilleInterface->computeAlignmentGlobal(state); 
                _MilleInterface->setGlobalLabels(state);  
                ///You can then access them
                TMatrixD const& alignmentJacobian = _MilleInterface->getAlignmentJacobian();
                std::vector<int> labels =  _MilleInterface->getGlobalParameters();
                ///Now fill GBL point with global derivatives.
                pointList.at(j).addGlobals(labels, alignmentJacobian);
                streamlog_out(DEBUG0)<<"The alignment matrix after adding to point: "<<std::endl;
                streamlog_message( DEBUG0, pointList.at(j).getGlobalDerivatives().Print() ;, std::endl; );
            }
        }
        
    }

	//This set the estimate resolution for each plane in the X direction.
	void EUTelGBLFitter::setParamterIdXResolutionVec( const std::vector<float>& vector)
	{
		//We have a similar check after this to see that number of planes and elements in resolution vector are the same. We need this here since if they are different then it will just give an exception from the vector tryign to access a element that does not exist.
		if ( geo::gGeometry().sensorZOrdertoIDs().size() != vector.size() ){
			streamlog_out( ERROR5 ) << "The number of planes: "<< geo::gGeometry().sensorZOrdertoIDs().size()<<"  The size of input resolution vector: "<<vector.size()  << std::endl;
			throw(lcio::Exception("The size of the resolution vector and the total number of planes is different for x axis."));
		}

		for(size_t i=0; i < geo::gGeometry().sensorZOrdertoIDs().size(); ++i){
			_parameterIdXResolutionVec[ geo::gGeometry().sensorZOrderToID(i)] = vector.at(i);
		}
	}
	//This sets the estimated resolution for each plane in the Y direction.
	void EUTelGBLFitter::setParamterIdYResolutionVec( const std::vector<float>& vector)
	{
		if ( geo::gGeometry().sensorZOrdertoIDs().size() != vector.size() ){
			streamlog_out( ERROR5 ) << "The number of planes: "<< geo::gGeometry().sensorZOrdertoIDs().size()<<"  The size of input resolution vector: "<<vector.size()  << std::endl;
			throw(lcio::Exception("The size of the resolution vector and the total number of planes is different for y axis."));
		}
		for(size_t i=0; i < geo::gGeometry().sensorZOrdertoIDs().size(); ++i){
			_parameterIdYResolutionVec[ geo::gGeometry().sensorZOrderToID(i)] = vector.at(i);
		}
	}

	//This is used to deal with downweighting of ouliers. You must provide as input t,h or c. This specifies the function that will be used to do the downweigting.
	//TO DO: Check that this function works as expected since it has never been used. 
	void EUTelGBLFitter::setMEstimatorType( const std::string& mEstimatorType ) {
		std::string mEstimatorTypeLowerCase = mEstimatorType;
		std::transform( mEstimatorType.begin(), mEstimatorType.end(), mEstimatorTypeLowerCase.begin(), ::tolower);//Make the character lower case
		if ( mEstimatorType.size() != 1 ) {
			streamlog_out( WARNING1 ) << "More than one character supplied as M-estimator option" << std::endl;
			streamlog_out( WARNING1 ) << "No M-estimator downweighting will be used" << std::endl;
			return;
		}
		//Compare character to the ones that are accepted and if one is the same then set out member variable equal to it.
		if ( mEstimatorTypeLowerCase.compare("t") == 0 ||
			mEstimatorTypeLowerCase.compare("h") == 0 ||
			mEstimatorTypeLowerCase.compare("c") == 0   ) this->_mEstimatorType = _mEstimatorType;
		else {
			streamlog_out( WARNING1 ) << "M-estimator option " << mEstimatorType << " was not recognized" << std::endl;
			streamlog_out( WARNING1 ) << "No M-estimator downweighting will be used" << std::endl;
		}
	}
    void EUTelGBLFitter::setRadLengths(EUTelTrack & track,	std::map<const int,double>&  mapSensor, std::map<const int ,double>&  mapAir, double & rad ){
        //THE FINAL WEIGHT WE HAVE WILL BE A FRACTION PERCENTAGE OF THE TOTAL RADIATION LENGTH
        std::vector<EUTelState>& states = track.getStates();
        const double var  = pow( Utility::getThetaRMSHighland(track.getBeamEnergy(), rad) , 2);
        streamlog_out(DEBUG5)<< std::scientific << "This is the total variance of the track and energy: " << var << track.getBeamEnergy() <<std::endl; 
        for(size_t i =0; i < track.getStates().size();++i){ 
            streamlog_out(DEBUG0)<< std::scientific << " Values placed in variance using Highland formula corrected. (SENSOR) : " << (mapSensor[states.at(i).getLocation()]/rad)*var << "  (AIR)  " << (mapAir[states.at(i).getLocation()]/rad)*var <<std::endl;
            states.at(i).setRadFrac((mapSensor[states.at(i).getLocation()]/rad)*var,(mapAir[states.at(i).getLocation()]/rad)*var);//We input the fraction percentage.
        }
        track.setTotalVariance(var);
    }
    void EUTelGBLFitter::setRad(EUTelTrack & track){
		streamlog_out(DEBUG2)<<"setRadLength--BEGIN"<<std::endl;
        ///Set radiation lengths
        std::map<const int,double>  mapSensor;
        std::map<const int ,double>  mapAir;
        double rad = track.getStates().at(0).computeRadLengthsToEnd(mapSensor, mapAir);
		streamlog_out(DEBUG2)<<"Rad total: " << rad<<std::endl;
        if(rad == 0){
      //      throw(std::string("Radiation length is zero for mimosa tracks."));
              throw(std::string(""));

            _numberRadLoss++;
        }
		streamlog_out(DEBUG2)<<"Radiation lengths for this tracks geometry correcly allocated. Tracks lost:  " << _numberRadLoss <<std::endl;

        setRadLengths(track, mapSensor, mapAir, rad);
    }
	std::map<  unsigned int, std::vector<double> > EUTelGBLFitter::getScatPos(EUTelTrack& track) const {}

	void EUTelGBLFitter::getBasicList(EUTelTrack& track, std::vector< gbl::GblPoint >& pointList, std::map<  unsigned int, unsigned int>  & linkGL){
        TMatrixD jac(5, 5);
        jac.UnitMatrix();
        /// Label begins at one. 
        unsigned int label = 1;
        gbl::GblPoint point(jac);
        ///Set label, make pair, then increase. 
        point.setLabel(label); // 1
        linkGL.insert(std::make_pair(track.getStates().at(0).getLocation(), label));
        label++;
        pointList.push_back(point);
        linkGL.insert(std::make_pair(track.getStates().at(0).getLocation(), label));
        for(size_t i=0;i < (track.getStates().size()-1); i++){		
            EUTelState state = track.getStates().at(i);
            ///Need this for the location.
            EUTelState stateNext = track.getStates().at(i+1);
            streamlog_out(DEBUG5) << "The jacobian which will be used at label " << label << std::endl;
            TMatrixD jac = EUTelNav::getPropagationJacobianGlobalToGlobal(state.getArcLengthToNextState(), state.getDirGlobal());
            streamlog_message( DEBUG0, jac.Print();, std::endl; );
            gbl::GblPoint point(jac);
            point.setLabel(label);
            linkGL.insert(std::make_pair(stateNext.getLocation(), label));
            label++;
            pointList.push_back(point);
        } 
    }   
	void EUTelGBLFitter::getMeas(EUTelTrack& track, std::vector< gbl::GblPoint >& pointList, std::map< unsigned int,unsigned int > & linkGL , std::map< unsigned int, unsigned int>  & linkMeas){
        for(size_t i=0;i < track.getStates().size(); i++){		
            EUTelState state = track.getStates().at(i);
            //Use state to get correct label
            streamlog_out(DEBUG5) << "State location " << state.getLocation() << std::endl;
            if(state.getStateHasHit()){
                streamlog_out(DEBUG5) << "State has hit" << std::endl;
                unsigned int label = linkGL[state.getLocation()];
                setMeasurementGBL(pointList.at(label-1), state );
                linkMeas.insert(std::make_pair(state.getLocation(),label));
            }
        }
        streamlog_out(DEBUG5) << "Link size: " << linkMeas.size() << std::endl;

    }
	void EUTelGBLFitter::getScat(EUTelTrack& track, std::vector< gbl::GblPoint >& pointList, std::map< unsigned int,unsigned int > & linkGL){
        for(size_t i;i < track.getStates().size(); i++){		
            EUTelState state = track.getStates().at(i);
            //Use state to get correct label
            unsigned int label = linkGL[state.getLocation()];
            setScattererGBL(pointList.at(label-1),state ); 
        }
    }

    ///This is the work horse of the GBL fitter. It creates GBL points from EUTelTracks and returns the relations between the two.
	void EUTelGBLFitter::getGBLPointsFromTrack(EUTelTrack& track, std::vector< gbl::GblPoint >& pointList, std::map<  unsigned int,unsigned int > & linkGL,std::map< unsigned int, unsigned int > & linkMeas ){
		streamlog_out(DEBUG0)<<"EUTelGBLFitter::getGBLPointsFromTrack-------------------------------------BEGIN"<<std::endl;
        ///get radiation length of the track. This is added to the track internally.
        setRad(track);
        /// At the moment this does not pass anything
//        std::map<  unsigned int, std::vector<double> > scatPos = getScatPos(track);
        /// This is the minimum needed to create a GBL trajectory. It basically only relates each GBL point to each other at this stage.
		streamlog_out(DEBUG0)<<"Create basic traj. "<<std::endl;
        getBasicList(track,pointList,linkGL);
        /// Now add measurement. Connect from local to global will be made here. 
		streamlog_out(DEBUG0)<<"Add measurement. "<<std::endl;
        getMeas(track,pointList,linkGL,linkMeas);
        /// Now add scattering.  
		streamlog_out(DEBUG0)<<"Add scattering information. "<<std::endl;
        getScat(track,pointList,linkGL);
		streamlog_out(DEBUG0)<<"EUTelGBLFitter::getGBLPointsFromTrack-------------------------------------END"<<std::endl;

	}

  void EUTelGBLFitter::getResLoc(gbl::GblTrajectory* traj, std::map< unsigned int, unsigned int> & linkMeas, std::vector< gbl::GblPoint > pointList,std::map< int, std::map< float, float > > &  SensorResidual, std::map< int, std::map< float, float > >& sensorResidualError){
    for(std::map<unsigned int, unsigned int>::iterator iter = linkMeas.begin(); iter != linkMeas.end(); ++iter){
        unsigned int numData; 
        TVectorD aResiduals(2);
        TVectorD aMeasErrors(2);
        TVectorD aResErrors(2);
        TVectorD aDownWeights(2); 
        traj->getMeasResults(iter->second, numData, aResiduals, aMeasErrors, aResErrors, aDownWeights);
        streamlog_out(DEBUG0) <<"State location: "<<iter->first<<" The residual x " <<aResiduals[0]<<" The residual y " <<aResiduals[1]<<std::endl;
        std::map<float, float> res; //This is create on the stack but will pass thisa by value to the new map so it 
        res.insert(std::make_pair(aResiduals[0],aResiduals[1]));
        SensorResidual.insert(std::make_pair(iter->first, res));		
        std::map<float, float> resError; //This is create on the stack but will pass thisa by value to the new map so it 
        resError.insert(std::make_pair(aResErrors[0],aResErrors[1]));
        sensorResidualError.insert(std::make_pair(iter->first, resError));	
	}
  }
	std::string EUTelGBLFitter::getMEstimatorType( ) const {
			return _mEstimatorType;
	}
	//COMPUTE
	void EUTelGBLFitter::computeTrajectoryAndFit(gbl::GblTrajectory* traj, double* chi2, int* ndf, int & ierr){
		streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::computeTrajectoryAndFit-- BEGIN " << std::endl;
		double loss = 0.;
		streamlog_out ( DEBUG0 ) << "This is the trajectory we are just about to fit: " << std::endl;
		streamlog_message( DEBUG0, traj->printTrajectory(10);, std::endl; );
		streamlog_out ( DEBUG0 ) << "This is the points in that trajectory " << std::endl;
		streamlog_message( DEBUG0, traj->printPoints(10);, std::endl; );


		if ( !_mEstimatorType.empty( ) ) ierr = traj->fit( *chi2, *ndf, loss, _mEstimatorType );
		else ierr = traj->fit( *chi2, *ndf, loss );

		if( ierr != 0 ){
			streamlog_out(MESSAGE0) << "Fit failed!" << " Track error: "<< ierr << " and chi2: " << *chi2 << std::endl;
		}
		else{
            streamlog_out(MESSAGE0) << "Fit Successful!" << " Track error; "<< ierr << " and chi2: " << *chi2 << std::endl;
		}
		streamlog_out ( DEBUG4 ) << " EUTelGBLFitter::computeTrajectoryAndFit -- END " << std::endl;
	}
	//TEST
	void EUTelGBLFitter::testUserInput(){
		if(_parameterIdXResolutionVec.size() != _parameterIdYResolutionVec.size()){
				throw(lcio::Exception("The vector for resolutions for X and Y are different sizes."));
		}
		if(_parameterIdXResolutionVec.size() != geo::gGeometry().nPlanes() ){
				throw(lcio::Exception("The total number of planes and the resolution of the planes vector are different sizes."));
		}
	}
	void EUTelGBLFitter::testTrack(EUTelTrack& track){
		streamlog_out(DEBUG4)<<"EUTelGBLFitter::testTrack------------------------------------BEGIN"<<std::endl;
		if(track.getStates().size() == 0 ){
			throw(lcio::Exception("The number of states is zero."));
		}
		///Note we do not use excluded planes here. This should be dealt with in pattern recognition.
		if (track.getNumberOfHitsOnTrack() > geo::gGeometry().nPlanes() ){
			throw(lcio::Exception("The number of hits on the track is greater than the number of planes.")); 	
		}
		streamlog_out(DEBUG5)<< "Input track passed tests!" <<std::endl;
		streamlog_out(DEBUG4)<<"EUTelGBLFitter::testTrack------------------------------------END"<<std::endl;

	} 
	void EUTelGBLFitter::getCorr(gbl::GblTrajectory* traj,EUTelTrack &track, std::map<unsigned int ,unsigned int >& linkGL , std::map<int, std::vector<double> > &  mapSensorIDToCorrectionVec){
        ///Only state which are created before will be updated. Scattering planes within GBL are not saved.
		for(size_t i = 0;i < track.getStates().size(); i++){		
			EUTelState& state = track.getStates().at(i);
			TVectorD corrections(5);
			TMatrixDSym cov(5);
            /// Get the corrections in the global frame!!!! 
            /// This is corrected internally by EUTelTrack and EUTelState.
            traj->getResults(linkGL[state.getLocation()], corrections, cov );
            streamlog_out(DEBUG3) << std::endl << "State before we have added corrections: " << std::endl;
            state.print();
            streamlog_out(DEBUG3) << std::endl << "Correction: " << std::endl;
            streamlog_message( DEBUG3, corrections.Print();, std::endl; );			
            state.setStateUsingCorrection(corrections);
            track.setTrackUsingCorrection(corrections);
            state.setCov(cov);
            unsigned int numData;
            /// Scattering is for every plane and is added here. 
            ///Measurement - Prediction is the residual. Initial M-P is always 0
            TVectorD aResidualsKink(2);
            TVectorD aMeasErrorsKink(2);
            TVectorD aResErrorsKink(2);
            TVectorD aDownWeightsKink(2); 
            traj->getScatResults( linkGL[state.getLocation()], numData, aResidualsKink, aMeasErrorsKink, aResErrorsKink, aDownWeightsKink);
            state.setKinks(aResidualsKink);
        }
    }
}


#endif
