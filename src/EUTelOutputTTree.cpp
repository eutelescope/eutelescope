#include "EUTELESCOPE.h"
#include "EUTelEventImpl.h"
#include "EUTelExceptions.h"
#include "EUTelRunHeaderImpl.h"
#include "EUTelTrackerDataInterfacerImpl.h"

#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelGenericPixGeoDescr.h"

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <UTIL/CellIDDecoder.h>

#include <algorithm>
#include "marlin/Processor.h"
#include <iostream>

#include "EVENT/TrackerPulse.h"

using namespace eutelescope;
using namespace std;





#include "marlin/Processor.h"

// system includes <>
#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include "EUTelTrack.h"
#include "EUTelHit.h"
#include "EUTelReaderGenericLCIO.h"


class GBL_trackOutput;
TTree* m_tree = new TTree("gbl", "gbl");

namespace eutelescope {


    class EUTelOutputTTree : public  marlin::Processor{

        public:
        virtual Processor*  newProcessor() { return new EUTelOutputTTree; }
        std::string _trackCol;
        EUTelOutputTTree();
        ~EUTelOutputTTree(){}
        virtual void init();
        virtual void processRunHeader(LCRunHeader* run);
        virtual void processEvent(LCEvent * evt);
        virtual void check(LCEvent * /*evt*/){ ; };
        virtual void end();

        protected:
        GBL_trackOutput* m_gbl;
      //  TTree *m_tree;
     //   TFile *m_file;
        std::string path;
        private:
    };
}




class GBL_trackOutput {
    public:

    GBL_trackOutput(const std::string& name )  {
    m_tree->Branch("ID", &m_id);
    m_tree->Branch("x", &m_x);
    m_tree->Branch("y", &m_y);
    m_tree->Branch("chi2", &m_chi2);
    m_tree->Branch("ndf", &m_ndf);
    m_tree->Branch("phi", &m_phi);
    m_tree->Branch("theta", &m_theta);

    m_tree->Branch("event_nr", &m_event_nr);
}
    ~GBL_trackOutput(){
      }
    virtual void pushCollection( EVENT::LCEvent* ev, std::string colName ) {
        beginEvent();
        std::vector<EUTelTrack> tr = lc_reader.getTracks(ev, colName);

        for (size_t i = 0; i < tr.size();++i){

            processTrack(tr[i]);
        }
        endEvent();
    };
    void beginEvent(){
        m_x.clear();
        m_y.clear();
        m_id.clear();
        m_phi.clear();
        m_theta.clear();
        m_ndf.clear();
        m_chi2.clear();
        ++m_event_nr;

    }
    void endEvent(){
        m_tree->Fill();
//        m_tree->Print();

    }
    void processTrack(EUTelTrack& trc){

        std::vector<EUTelState>& planes = trc.getStates();

        chi2 = trc.getChi2();
        ndf = trc.getNdf();

        for (size_t i = 0; i <  planes.size(); ++i){

            processPlanes(planes[i]);
        }

    }
    void processPlanes(EUTelState& pln){

        x = pln.getPosition()[0];
        y = pln.getPosition()[1];

        ID = static_cast<double>(pln.getLocation());

        phi = pln.getDirLocalX() / pln.getDirLocalZ();
        theta = pln.getDirLocalY() / pln.getDirLocalZ();
        pushHit();
    }
    void pushHit(){
        m_x.push_back(x);
        m_y.push_back(y);
        m_id.push_back(ID);
        m_phi.push_back(phi);
        m_theta.push_back(theta);
        m_ndf.push_back(ndf);
        m_chi2.push_back(chi2);
    }


    double x;
    double y;
    double ID;
    double phi;
    double theta;
    double chi2;
    double ndf;
    EUTelReaderGenericLCIO lc_reader;
    int m_event_nr;
    std::string m_name, m_type;
    std::vector<double> m_x, m_y, m_id,
    m_chi2,m_ndf,m_phi,m_theta;
};
EUTelOutputTTree::EUTelOutputTTree() :Processor("EUTelOutputTTree"),m_gbl(NULL){
    registerProcessorParameter("OutputPath", "Path/File where root-file should be stored",path, std::string("NTuple.root"));
    registerProcessorParameter("TrackCol", "Input track collection name",_trackCol,std::string("TrackCol"));

}
void EUTelOutputTTree::init()
{

    std::string name("test.root");
    geo::gGeometry().initializeTGeoDescription(name, true);
    m_gbl = new GBL_trackOutput("GBL_tracks");
}

void EUTelOutputTTree::processRunHeader(LCRunHeader* run)
{

}

void EUTelOutputTTree::processEvent(LCEvent * evt)
{
    try{

        m_gbl->pushCollection(evt,_trackCol );
    }catch(...){
        throw marlin::SkipEventException(this);
    }
}

void EUTelOutputTTree::end()
{
            TFile* file =  new  TFile("fitter.root", "RECREATE");

            m_tree->Write();

            file->Write();

            file->Close();
            delete file;
            delete m_tree;

  //  if (gFile_)
  //      {
  //          std::cout << "closing  file " << std::endl;
  //          gFile_->Write();
  //          gFile_->Close();
  //          delete gFile_;
  //          gFile_ = NULL;
  //      }
        if (m_gbl)
        {

            delete m_gbl;
            m_gbl = NULL;
        }
}


EUTelOutputTTree gEUTelOutputTTree;
