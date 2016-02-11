/*
 * EUTelProcessorTrackSelection.h 
 * 
 * Created on: April 19th 2015 
 *     author:Alexander Morton 
 * 
 * This processor will take in GBLTrack objects and remove tracks not meeting all cut requirements.  
 *
 */




// LCIO
#include <EVENT/LCCollection.h>
#include "lcio.h"

// MARLIN
#include "marlin/Exceptions.h"
#include "marlin/Global.h"
#include "marlin/Processor.h"
#include "marlin/VerbosityLevels.h"

//EUTelescope
#include "EUTelUtility.h"
#include "EUTelEventImpl.h"
#include "EUTelTrack.h"
#include "EUTelState.h"
#include "EUTelTrackSelection.h"
#include "EUTelGBLFitter.h"
#include "EUTelReaderGenericLCIO.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "TAxis.h"
#include <TFile.h>
#include <TH1F.h>
#include <THStack.h>
#include <TH1.h>
#include <TGaxis.h>
#include <TCanvas.h> 
#include "TStyle.h"
#include "TLegend.h"
#include "TPad.h"
#include "TPaveStats.h"

namespace eutelescope {

	class  EUTelProcessorTrackSelection : public marlin::Processor {

  	private:
  	DISALLOW_COPY_AND_ASSIGN(EUTelProcessorTrackSelection)
        
    public:

    virtual Processor* newProcessor() {
    	return new  EUTelProcessorTrackSelection;
    }

    EUTelProcessorTrackSelection();

    /** Called at the begin of the job before anything is read.
    * Use to initialize the processor, e.g. book histograms.
    */
    virtual void init();

    /** Called for every event - the working horse.*/
    virtual void processEvent(LCEvent * evt);

    /** Called after data processing for clean up. **/
		virtual void end();

    void outputLCIO(LCEvent* evt, std::vector<EUTelTrack>  track);
    bool stateSelection(EUTelState & state);
    void trackPlot(EUTelTrack & track,bool& pass);
    bool trackSelection(EUTelTrack  track);

//    EUTelTrackSelection* _selector;
    unsigned int _trackCountPlot;
    std::string _outHist; 
    std::string _trackInputCollectionName;
    std::string _tracksOutputCollectionName;
    std::vector<int> _sensors;
    std::vector<float> _events;
    std::vector<float> _chi2;
    std::vector<float> _event;
    std::vector<float> _curv;
    std::vector<float> _weightsX;
    std::vector<float> _weightsY;
    std::vector<float> _residualX;
    std::vector<float> _residualY;
    std::vector<float> _positionX;
    std::vector<float> _positionY;
    std::vector<float> _angleX;
    std::vector<float> _angleY;
    std::vector<float> _kinksX;
    std::vector<float> _kinksY;
    std::vector<float> _covSlopeX;
    std::vector<float> _covSlopeY;
    std::vector<float> _covPosX;
    std::vector<float> _covPosY;
    bool _hit;
    std::vector<float> _residualErrorX;
    std::vector<float>  _residualErrorY;
    std::vector<std::string> _labels;
    std::unique_ptr<TFile> _outFile;
    std::pair< std::unique_ptr<TH1F> , std::unique_ptr<TH1F> > _chi2Hists;
    std::map< std::string,  std::pair< std::unique_ptr<TH1F> , std::unique_ptr<TH1F> > > _hists;
    const bool compareXYToCut( double X, double Y, float cutLowX, float cutHighX, float cutLowY, float cutHighY){
        streamlog_out(DEBUG5)<<"Compare " << X << " " << Y << " "  <<std::endl;
        if( cutLowX  > X or X > cutHighX){
            return false;
        }
        if( cutLowY  > Y or Y > cutHighY){
            return false;

        }
        streamlog_out(DEBUG5)<<"true"  <<std::endl;
        return true;
    }
    const bool compareXYToCut( float & X, float & Y, float & cutLowX, float & cutHighX, float & cutLowY, float & cutHighY){
        streamlog_out(DEBUG5)<<"Compare " << X << " " << Y << " "  <<std::endl;
        if( cutLowX  > X or X > cutHighX){
            return false;
        }
        if( cutLowY  > Y or Y > cutHighY){
            return false;
        }
        streamlog_out(DEBUG5)<<"true"  <<std::endl;
        return true;
    }
    const bool compareXYToCut( int X, float & cutLowX, float & cutHighX){
        if( cutLowX  > X or X > cutHighX){
            return false;
        }
        streamlog_out(DEBUG5)<<"true"  <<std::endl;
        return true;
    }

    const bool compareXYToCut( float & X, float & cutLowX, float & cutHighX){
        if( cutLowX  > X or X > cutHighX){
            return false;
        }
        streamlog_out(DEBUG5)<<"true"  <<std::endl;
        return true;
    }
    const bool compareXYToCut( double X, float & cutLowX, float & cutHighX){
        if( cutLowX  > X or X > cutHighX){
            return false;
        }
        streamlog_out(DEBUG5)<<"true"  <<std::endl;
        return true;
    }


	};

    EUTelProcessorTrackSelection gEUTelProcessorTrackSelection;
}
