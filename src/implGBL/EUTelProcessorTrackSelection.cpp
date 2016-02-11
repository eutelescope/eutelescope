#include "EUTelProcessorTrackSelection.h"

using namespace eutelescope;
/// Combined input can create a particular type of track, with hits on certain planes, and a particular chi2.
/**
 * \param [in] TrackInputCollectionName Name of collection containing GBL tracks.
 * \param [in] TrackOutputCollectionName Name of collection containing GBL tracks to output.
 * \param [in] sensorsIDsToPass The sensors must pass these selection. 
 * \param [in] sensorsIDsMustNotHaveHit If plane specified then this plane must not have a hit.
 */
EUTelProcessorTrackSelection::EUTelProcessorTrackSelection() :
Processor("EUTelProcessorTrackSelection") {
    registerInputCollection(LCIO::TRACK, "TrackInputCollectionName", "Input track collection name",_trackInputCollectionName,std::string("TrackCandidatesCollection"));
    registerOutputCollection(LCIO::TRACK,"TracksOutputCollectionName","Output tracks collection name",_tracksOutputCollectionName,std::string("TrackCollection"));
    registerOptionalParameter("outhist", "Name and location of output histogram ", _outHist, std::string("selection.root"));
	registerOptionalParameter("sensorsIDsToPass", "The must pass selection",_sensors,std::vector<int> () );
	registerOptionalParameter("chi2", "These sensors must have a hit to pass selection",_chi2,std::vector<float> () );
	registerOptionalParameter("events", "Cut ranges for variable.",_event,std::vector<float> () );
	registerOptionalParameter("curv", "Cut ranges for variable.",_curv,std::vector<float> () );
	registerOptionalParameter("weightsX ", "Cut ranges for variable.",_weightsX,std::vector<float> () );
	registerOptionalParameter("weightsY", "Cut ranges for variable.",_weightsY,std::vector<float> () );
	registerOptionalParameter("residualX", "Cut ranges for variable.",_residualX,std::vector<float> () );
	registerOptionalParameter("residualY", "Cut ranges for variable.",_residualY,std::vector<float> () );
    registerOptionalParameter("residualErrorX", "Cut ranges for variable.",_residualErrorX,std::vector<float> () );
    registerOptionalParameter("residualErrorY", "Cut ranges for variable.",_residualErrorY,std::vector<float> () );
	registerOptionalParameter("positionX", "Cut ranges for variable.",_positionX,std::vector<float> () );
    registerOptionalParameter("positionY", "Cut ranges for variable.",_positionY,std::vector<float> () );
    registerOptionalParameter("angleX", "Cut ranges for variable.",_angleX,std::vector<float> () );
    registerOptionalParameter("angleY", "Cut ranges for variable.",_angleY,std::vector<float> () );
    registerOptionalParameter("kinksX", "Cut ranges for variable.",_kinksX,std::vector<float> () );
    registerOptionalParameter("kinksY", "Cut ranges for variable.",_kinksY,std::vector<float> () );
    registerOptionalParameter("hit", "Cut ranges for variable.",_hit,bool(false) );
    registerOptionalParameter("covSlopeX", "Cut ranges for variable.",_covSlopeX,std::vector<float> () );
    registerOptionalParameter("covSlopeY", "Cut ranges for variable.",_covSlopeY,std::vector<float> () );
    registerOptionalParameter("covPosX", "Cut ranges for variable.",_covPosX,std::vector<float> () );
    registerOptionalParameter("covPosY", "Cut ranges for variable.",_covPosY,std::vector<float> () );

}
           // _hists[key].first.get()->Fill(state.getRes().at(0));
           // _hists[key].first.get()->Fill(state.getRes().at(1));
           // _hists[key].first.get()->Fill(state.getKinks()[0]);
           // _hists[key].first.get()->Fill(state.getKinks()[1]);
           // _hists[key].first.get()->Fill(state.getSlopeX());
           // _hists[key].first.get()->Fill(state.getSlopeY());
           // _hists[key].first.get()->Fill(state.getPosition()[1]);
           // _hists[key].first.get()->Fill(state.getPosition()[0]);
           // _hists[key].first.get()->Fill(state.getHit().getWeight().at(0));
           // _hists[key].first.get()->Fill(state.getHit().getWeight().at(1));


void EUTelProcessorTrackSelection::init(){
	try{
        streamlog_out(DEBUG5)<<"Init selection...." <<std::endl;
        _trackCountPlot= 0;
        std::vector<int> none = {};
        EUTelExcludedPlanes::setRelativeComplementSet(none);
       // _selector = new EUTelTrackSelection();
		std::string name("test.root");
//		geo::gGeometry().initializeTGeoDescription(name,false);
        //Is this a good idea? Could confuse ownership?
        std::unique_ptr<TFile> outfile( new TFile(_outHist.c_str(),"RECREATE"));
        _outFile =std::move(outfile);
        std::unique_ptr<TH1F> chi2Pass( new TH1F("PassChi2Ndof","Pass Chi2 Normalised", 100,0,20) );
        chi2Pass.get()->SetFillColor(kGreen);
        std::unique_ptr<TH1F> chi2Fail( new TH1F("FailChi2Ndof","Fail Chi2 Normalised", 100,0,20) );
        chi2Fail.get()->SetFillColor(kRed);
        _chi2Hists = std::move(std::make_pair<std::unique_ptr<TH1F>,std::unique_ptr<TH1F>>(std::move(chi2Pass),std::move(chi2Fail)));
        std::vector <std::string> labels{"Residual X ", "Residual Y ","KinkX ","KinkY ","SlopeX ","SlopeY ","PosX ","PosY ", "Weight X ","Weight Y ","Residual Error X ", "Residual Error Y ", "Cov Sqrt Y ", "Cov Sqrt X ","Cov Sqrt Slope Y ","Cov Sqrt Slope X" };
        //Add number to label to get sensor ID.
        for(size_t i = 0 ; i < EUTelExcludedPlanes::_senInc.size(); ++i){
            for(auto label: labels){ 
                streamlog_out(DEBUG5)<<"label: " << label << "Sen " << std::to_string(EUTelExcludedPlanes::_senInc.at(i)) <<std::endl;
                std::string titlePass ="Passed " + label + " of Sensor " +  std::to_string(EUTelExcludedPlanes::_senInc.at(i));
                std::string titleFail ="Failed " + label + " of Sensor " +  std::to_string(EUTelExcludedPlanes::_senInc.at(i));
                std::string namePass = titlePass;
                std::string nameFail = titleFail;
                namePass.erase( std::remove_if( namePass.begin(), namePass.end(), ::isspace ), namePass.end() );
                nameFail.erase( std::remove_if( nameFail.begin(), nameFail.end(), ::isspace ), nameFail.end() );
                std::unique_ptr<TH1F> histPass( new TH1F(namePass.c_str(),titlePass.c_str(), 10000,0,0) );
                histPass.get()->SetFillColor(kGreen);
                std::unique_ptr<TH1F> histFail( new TH1F(nameFail.c_str(),titleFail.c_str(), 10000,0,0) );
                histFail.get()->SetFillColor(kRed);
                if(i == 0){///Silly but we only need it for one run. We attach the number again as id in processEvent.
                    _labels.push_back(label);
                }
                label = label + std::to_string(EUTelExcludedPlanes::_senInc.at(i));
                _hists[label] =std::make_pair<std::unique_ptr<TH1F>,std::unique_ptr<TH1F>>(std::move(histPass),std::move(histFail));
//                std::cout << "Sen: " << i <<std::endl;
            }
        }
        }catch(...){	
            streamlog_out(ERROR5)<<"There is an unknown error in EUTelProcessorTrackSelection-init" <<std::endl;
            throw marlin::StopProcessingException( this ) ;

        }

}

void EUTelProcessorTrackSelection::processEvent(LCEvent * evt){
//    std::cout << "Event" <<std::endl;
	try{
		EUTelEventImpl * event = static_cast<EUTelEventImpl*> (evt); ///We change the class so we can use EUTelescope functions

		if (event->getEventType() == kEORE) {
			streamlog_out(DEBUG4) << "EORE found: nothing else to do." << std::endl;
			return;
		}else if (event->getEventType() == kUNKNOWN) {
			streamlog_out(WARNING2) << "Event number " << event->getEventNumber() << " in run " << event->getRunNumber() << " is of unknown type. Continue considering it as a normal Data Event." << std::endl;
		}
        std::vector<EUTelTrack> tracksOut;
        EUTelReaderGenericLCIO reader = EUTelReaderGenericLCIO();
        std::vector<EUTelTrack> tracks = reader.getTracks(evt, _trackInputCollectionName);
        for (size_t i = 0; i < tracks.size(); ++i){
//            std::cout << "tracks in "<< tracks.size() <<std::endl;
//            std::cout << "Tracks: " << i <<std::endl;
            streamlog_out(DEBUG5)<<"Track " << i << " Size: "<<tracks.size() <<" Planes: " << geo::gGeometry().sensorIDsVec().size() << std::endl;
            bool pass = false;
            EUTelTrack track = tracks.at(i); 
            streamlog_out(DEBUG5)<<"Selection... track chi2  " << track.getChi2() << " " << std::endl;
            //////////////////////////////
            // Is in event range??
            /////////////////////////////////
            streamlog_out(DEBUG5)<<"Found "<< event->getEventNumber() << " events range " << _event.at(0) << " " << _event.at(1)  <<std::endl;
            const bool inEventRange = compareXYToCut(event->getEventNumber() , _event.at(0), _event.at(1));
            if(inEventRange){
                pass = trackSelection(track);
//                std::cout << "Pass value: " << pass <<std::endl;
            }else{
                pass = false;
            }
//            std::cout << "Passed " << i+1 <<std::endl;
            _trackCountPlot++;
            trackPlot(track,pass);

            if(pass){
                    tracksOut.push_back(track);
            }
        }
  //      std::cout << "tracks out "<< tracksOut.size() <<std::endl;

        outputLCIO(evt,tracksOut);
    }	
    catch (DataNotAvailableException e) {
		streamlog_out(DEBUG5) << "Data not avaliable skip event. " << std::endl;
		throw marlin::SkipEventException(this);
	}
	catch(std::string &e){
		streamlog_out(DEBUG5) << e << std::endl;
		throw marlin::SkipEventException( this ) ;
	}
	catch(lcio::Exception& e){
		streamlog_out(ERROR5) << e.what() <<std::endl;
		throw marlin::StopProcessingException( this ) ;
	}
}

void EUTelProcessorTrackSelection::end(){
    {
        std::cout << "The number of tracks passed to plotter: " <<_trackCountPlot <<std::endl;
        _outFile->cd();
        TDirectory* subDir = (TDirectory*) _outFile->mkdir("Global");
        subDir->cd();
        //Make the bins the same
//        chi2HistFailSameBinning =  _chi2Hists.first.get()->Rebin(24,"hnew",xbins);
//        ###########
//        ##
//        ## Do not need to rebin. Assume same axis.
//        ##
//        ###############
    //    _chi2Hists.second.get()->Write();
    //    _chi2Hists.first.get()->Write();
///        TH1F chi2Total("Chi2 Total","Chi2 Total",100,0,20);
        TCanvas ctot("TotChi2","Total Chi2", 1260, 500);
        ctot.cd();
        gStyle->SetOptStat(1100);
//        gStyle->SetOptFit(1111);
        TH1F totalChi("Total Chi2","Total Chi2", 100,0,20);
        totalChi.SetFillColor(kBlue);
        totalChi.Add(_chi2Hists.first.get());
        totalChi.Add(_chi2Hists.second.get());
//        totalChi.Fit("Chi2","L");
        totalChi.SetFillColor(kBlue);
        totalChi.Write();
        totalChi.Draw("F");

     //   gStyle->SetOptStat(00011210);
        double leg_pos[4] = {0.72, 0.68, 0.84, 0.78};///Should use array.
     //   TLegend leg(leg_pos[0],leg_pos[1],leg_pos[2],leg_pos[3]);
     //   leg.SetFillColor(0);
     //   leg.SetFillStyle(0);
    //    leg.SetBorderSize(0);
    //    leg.AddEntry(_chi2Hists.first.get(), "Pass", "f");
    //    leg.AddEntry(_chi2Hists.second.get(), "Fail", "f");
    //    leg.SetTextSize(0.05);
        TCanvas c("PFChi2","Pass/Fail Chi2", 1260, 500);
        c.cd();
        gStyle->SetOptStat(1100);
 //       gStyle->SetOptFit(1111);
        THStack stack("Chi2","Pass/Fail Chi2");
        stack.Add(_chi2Hists.first.get());
        stack.Add(_chi2Hists.second.get());
//   stack.GetXaxis()->SetTitle("Tracks");
        c.SetLeftMargin(0.13);
        c.SetRightMargin(0.075);
        c.SetGrid();
    //    stack.Write();
        stack.Draw();
//        leg.Draw();
        ///
     //   c.Update();
      //  gPad->Update();

      //  TPaveStats *tps1 = (TPaveStats*) _chi2Hists.first.get()->FindObject("stats");
      //  tps1->SetName("Pass");
      //  tps1->SetTextColor(kGreen);
      //  tps1->SetLineColor(kGreen);
      //  double X1 = tps1->GetX1NDC();
      //  double Y1 = tps1->GetY1NDC();
      //  double X2 = tps1->GetX2NDC();
      //  double Y2 = tps1->GetY2NDC();

      //  TPaveStats *tps2 = (TPaveStats*) _chi2Hists.second.get()->FindObject("stats");
      //  tps2->SetTextColor(kRed);
      //  tps2->SetLineColor(kRed);
      //  tps2->SetLineColor(kRed);
      //  tps2->SetX1NDC(X1);
      //  tps2->SetX2NDC(X2);
      //  tps2->SetY1NDC(Y1-(Y2-Y1));
      //  tps2->SetY2NDC(Y1);
      //  tps1->Draw("same");
      //  tps2->Draw("same");

        c.Write();
    }

    for(size_t i = 0 ; i < EUTelExcludedPlanes::_senInc.size(); ++i){
        _outFile->cd();
        std::string dir =  "Detector " + std::to_string(EUTelExcludedPlanes::_senInc.at(i));
        TDirectory* subDir = (TDirectory*) _outFile->mkdir(dir.c_str());
        subDir->cd();
        for(auto& label: _labels){
            std::string name = label + std::to_string(EUTelExcludedPlanes::_senInc.at(i));
            auto& hist = _hists[name];
            ///If you do not set a range the histogram entries will be kept in a buffer until 1000 entries. Then the entries are used to determine the X axis range.
            ///TH1::SetDefaultBufferSize will change the entries which the histogram will be filled. Before this the histogram has not entries.
            ///You can also empty the buffer into the histogram by calling BufferEmpty
            hist.first.get()->BufferEmpty();
            hist.second.get()->BufferEmpty();
            //Check we have stuff in the histograms
            if(hist.first.get()->GetEntries() > 1 and hist.second.get()->GetEntries() > 1){
                 std::cout<< "Hists do have enough hits to be plotted together. " << label <<" Entries " <<hist.first.get()->GetEntries() << " " <<hist.second.get()->GetEntries()  << std::endl;

            }else{
                 std::cout<< "Hists do not have enough hits to be plotted together. " << label <<" Entries " <<hist.first.get()->GetEntries() << " " <<hist.second.get()->GetEntries()  << std::endl;
          //      hist.first.get()->Write();
          //      hist.second.get()->Write();
                continue;
            }
            //Pass bin
            Int_t binsPass =   hist.first.get()->GetXaxis()->GetNbins(); 
            std::vector<Double_t> vbinsPass(binsPass);
            Double_t* xbinsPass = &vbinsPass[0];
            hist.first.get()->GetLowEdge(xbinsPass);
            Double_t upEdgePass  = hist.first.get()->GetXaxis()->GetBinUpEdge(hist.first.get()->GetXaxis()->GetNbins());
    //        std::cout<< "Num pass: " << binsPass << " low/high: " << vbinsPass.at(0) << " " << upEdgePass  <<std::endl;
            ///
            //Failed bin
            Int_t binsFail =   hist.second.get()->GetXaxis()->GetNbins(); 
            std::vector<Double_t> vbinsFail(binsFail);
            Double_t* xbinsFail = &vbinsFail[0];
            hist.second.get()->GetLowEdge(xbinsFail);
            Double_t upEdgeFail  = hist.second.get()->GetXaxis()->GetBinUpEdge(hist.second.get()->GetXaxis()->GetNbins());
  //          std::cout<< "Num fail: " << binsFail << " low/high: " << vbinsFail.at(0) << " " << upEdgeFail  <<std::endl;
            double min = std::min(vbinsPass.at(0),vbinsFail.at(0));
            double max = std::max(upEdgePass,upEdgeFail);
//            std::cout<< "Min/Max Edge"  << min << " " <<max << std::endl;
            if(min == max or min>max){
                std::cout<< "The histograms can not be create. Bin information is wrong" << std::endl;
                std::cout<< "Hist not produced " << label  << std::endl;
                continue;
            }
            ///Create new histograms
            TH1F newBinPass(hist.first.get()->GetTitle(),hist.first.get()->GetTitle(), 300,min,max);
            newBinPass.SetFillColor(kGreen);
            TH1F newBinFail(hist.second.get()->GetTitle(),hist.second.get()->GetTitle(), 300,min,max);
            newBinFail.SetFillColor(kRed);
            ///Assume fine binning so if you have an entry then add this to the histogram at the right location
            for(unsigned int i = 0 ;i<  hist.second.get()->GetXaxis()->GetNbins();i++){ 
                Double_t binContent  = hist.second.get()->GetBinContent(i);
                Double_t binCenter = hist.second.get()->GetXaxis()->GetBinCenter(i);
    //            std::cout<< "Bins Fail"  << i<< " " <<binContent << " " << binCenter << std::endl;
                for(unsigned int j = 0 ; j < binContent; j++){
//                    std::cout<<"found"<<std::endl;
                    newBinFail.Fill(binCenter);
                }
            }
            for(unsigned int i = 0 ;i<  hist.first.get()->GetXaxis()->GetNbins();i++){ 
                Double_t binContent  = hist.first.get()->GetBinContent(i);
                Double_t binCenter = hist.first.get()->GetXaxis()->GetBinCenter(i);
  //              std::cout<< "Bins Pass"  << i<< " " <<binContent << " " << binCenter << std::endl;
                for(unsigned int j = 0 ; j < binContent; j++){
                    newBinPass.Fill(binCenter);
                }
            }
            ///Set total for further analysis. 
            name =" Total " +  label + " Sensor " + std::to_string(EUTelExcludedPlanes::_senInc.at(i));
            TH1F totalVar(name.c_str(),name.c_str(), 300,min,max);
            totalVar.SetFillColor(kBlue);
            totalVar.Add(&newBinPass);
            totalVar.Add(&newBinFail);
            totalVar.Write();

            ///Now we can see the background after selection.
    //        newBinPass.Write();
    //        newBinFail.Write();
    //        hist.first.get()->Write();
    //        hist.second.get()->Write();
            streamlog_out(DEBUG4) << "Name of hist:" << name << std::endl;
            //Combine with number of get histograms again.
       //     gStyle->SetOptStat(00011210);
            double leg_pos[4] = {0.62, 0.48, 0.84, 0.88};///Should use array.
         //   TLegend leg(leg_pos[0],leg_pos[1],leg_pos[2],leg_pos[3]);
         //   leg.SetFillColor(0);
           // leg.SetFillStyle(4000);
          //  leg.SetBorderSize(0);
         //   leg.AddEntry(hist.first.get(), "Pass", "f");
        //    leg.AddEntry(rebinHist, "Fail", "f");
        //    leg.SetTextSize(0.025);
            name =" Pass/Fail " +  label + " Sensor " + std::to_string(EUTelExcludedPlanes::_senInc.at(i));
            THStack stack(name.c_str(),name.c_str());
            stack.Add(&newBinPass);
            stack.Add(&newBinFail);
            TCanvas c(name.c_str(),name.c_str(), 1260, 500);
            c.cd();
            gStyle->SetOptStat(1100);
            c.SetLeftMargin(0.13);
            c.SetRightMargin(0.075);
            c.SetGrid();
//            stack.Write();
            stack.Draw();
          //  leg.Draw();
            c.Write();
        }
        _outFile->cd();
    }
}

void EUTelProcessorTrackSelection::outputLCIO(LCEvent* evt, std::vector<EUTelTrack>  tracks)
{
    if(!tracks.empty()){
        for(unsigned int i=0 ; i< tracks.size(); i++){
            streamlog_out(DEBUG1)<<"Found "<<tracks.size()<<" track for event " << evt->getEventNumber() <<".  Track number  " << i <<std::endl;
            tracks.at(i).print();
        }
        EUTelReaderGenericLCIO reader = EUTelReaderGenericLCIO();
        reader.getColVec(tracks, evt, _tracksOutputCollectionName);
    }

}
void EUTelProcessorTrackSelection::trackPlot(EUTelTrack & track,bool& pass)
{
    if(pass){
        _chi2Hists.first.get()->Fill(track.getChi2()/track.getNdf());
    }else{
        _chi2Hists.second.get()->Fill(track.getChi2()/track.getNdf());
    }

    for(auto state : track.getStates()){
        TMatrixD cov = state.getHit().getCov();
        TMatrixD covState = state.getCov();
        if(pass){
            ////STUPID WAY!!!
            std::string key =  _labels.at(0)+std::to_string(state.getLocation());
            _hists[key].first.get()->Fill(state.getRes().at(0));
            key =  _labels.at(1)+std::to_string(state.getLocation());
            _hists[key].first.get()->Fill(state.getRes().at(1));
            key =  _labels.at(2)+std::to_string(state.getLocation());
            _hists[key].first.get()->Fill(state.getKinks()[0]);
            key =  _labels.at(3)+std::to_string(state.getLocation());
            _hists[key].first.get()->Fill(state.getKinks()[1]);
            key =  _labels.at(4)+std::to_string(state.getLocation());
            _hists[key].first.get()->Fill(state.getSlopeX());
            key =  _labels.at(5)+std::to_string(state.getLocation());
            _hists[key].first.get()->Fill(state.getSlopeY());
            key =  _labels.at(6)+std::to_string(state.getLocation());
            _hists[key].first.get()->Fill(state.getPosition()[1]);
            key =  _labels.at(7)+std::to_string(state.getLocation());
            _hists[key].first.get()->Fill(state.getPosition()[0]);
            key =  _labels.at(8)+std::to_string(state.getLocation());
            _hists[key].first.get()->Fill(state.getHit().getWeight().at(0));
            key =  _labels.at(9)+std::to_string(state.getLocation());
            _hists[key].first.get()->Fill(state.getHit().getWeight().at(1));
            ////
            key =  _labels.at(10)+std::to_string(state.getLocation());
            _hists[key].first.get()->Fill(sqrt(cov[0][0]));
            key =  _labels.at(11)+std::to_string(state.getLocation());
            _hists[key].first.get()->Fill(sqrt(cov[1][1]));
            key =  _labels.at(12)+std::to_string(state.getLocation());
            _hists[key].first.get()->Fill(sqrt(covState[4][4]));
            key =  _labels.at(13)+std::to_string(state.getLocation());
            _hists[key].first.get()->Fill(sqrt(covState[3][3]));
            key =  _labels.at(14)+std::to_string(state.getLocation());
            _hists[key].first.get()->Fill(sqrt(covState[2][2]));
            key =  _labels.at(15)+std::to_string(state.getLocation());
            _hists[key].first.get()->Fill(sqrt(covState[1][1]));



  //          std::cout << "The weights passed: "<< state.getHit().getWeight().at(0) << " " << state.getHit().getWeight().at(1)<<std::endl;
        }else{
            std::string key =  _labels.at(0)+std::to_string(state.getLocation());
            _hists[key].second.get()->Fill(state.getRes().at(0));
            key =  _labels.at(1)+std::to_string(state.getLocation());
            _hists[key].second.get()->Fill(state.getRes().at(1));
            key =  _labels.at(2)+std::to_string(state.getLocation());
            _hists[key].second.get()->Fill(state.getKinks()[0]);
            key =  _labels.at(3)+std::to_string(state.getLocation());
            _hists[key].second.get()->Fill(state.getKinks()[1]);
            key =  _labels.at(4)+std::to_string(state.getLocation());
            _hists[key].second.get()->Fill(state.getSlopeX());
            key =  _labels.at(5)+std::to_string(state.getLocation());
            _hists[key].second.get()->Fill(state.getSlopeY());
            key =  _labels.at(6)+std::to_string(state.getLocation());
            _hists[key].second.get()->Fill(state.getPosition()[1]);
            key =  _labels.at(7)+std::to_string(state.getLocation());
            _hists[key].second.get()->Fill(state.getPosition()[0]);
            key =  _labels.at(8)+std::to_string(state.getLocation());
            _hists[key].second.get()->Fill(state.getHit().getWeight().at(0));
            key =  _labels.at(9)+std::to_string(state.getLocation());
            _hists[key].second.get()->Fill(state.getHit().getWeight().at(1));
            ////
            key =  _labels.at(10)+std::to_string(state.getLocation());
            _hists[key].second.get()->Fill(sqrt(cov[0][0]));
            key =  _labels.at(11)+std::to_string(state.getLocation());
            _hists[key].second.get()->Fill(sqrt(cov[1][1]));
            key =  _labels.at(12)+std::to_string(state.getLocation());
            _hists[key].second.get()->Fill(sqrt(covState[4][4]));
            key =  _labels.at(13)+std::to_string(state.getLocation());
            _hists[key].second.get()->Fill(sqrt(covState[3][3]));
            key =  _labels.at(14)+std::to_string(state.getLocation());
            _hists[key].second.get()->Fill(sqrt(covState[2][2]));
            key =  _labels.at(15)+std::to_string(state.getLocation());
            _hists[key].second.get()->Fill(sqrt(covState[1][1]));
//            std::cout << "The weights failed: "<< state.getHit().getWeight().at(0) << " " << state.getHit().getWeight().at(1)<<std::endl;
        }
    }


}


bool EUTelProcessorTrackSelection::trackSelection(EUTelTrack  track)
{

    //////////////////////////////
    //In curvature/chi2 range?
    /////////////////////////////////
   streamlog_out(DEBUG1)<<"Track selection... "  <<std::endl;

    const bool inCurvRange = compareXYToCut( track.getChi2()/track.getNdf(), _curv.at(0), _curv.at(1));
    const bool inChiRange = compareXYToCut( track.getChi2()/track.getNdf(), _chi2.at(0), _chi2.at(1));

    if(inCurvRange and inChiRange){
        bool passAllStateChecks = true;
        for(auto state : track.getStates()){        
            auto inList = std::find(std::begin(_sensors), std::end(_sensors), state.getLocation());
            //Only apply cuts to planes specified.    
            if(inList != std::end(_sensors)) {
                bool statePass = stateSelection(state);
//                bool statePass=true;
                if(!statePass){
                    passAllStateChecks = false;
                }
            }
        }
        if(passAllStateChecks){
            return true;
        }
    }
    return false;
}
bool EUTelProcessorTrackSelection::stateSelection(EUTelState & state){
        //////////////////////////////
        //In state variable range?
        /////////////////////////////////
       streamlog_out(DEBUG5)<<"State... "<< std::endl;

        const bool inPosRange =  compareXYToCut( state.getPosition()[0],state.getPosition()[1], _positionX.at(0), _positionX.at(1),_positionY.at(0), _positionY.at(1));
        const bool inAngleRange =   compareXYToCut( state.getSlopeX(),state.getSlopeY(), _angleX.at(0), _angleX.at(1),_angleY.at(0), _angleY.at(1));
        const bool inKinksRange =   compareXYToCut( state.getKinks()[0],state.getKinks()[1], _kinksX.at(0), _kinksX.at(1),_kinksY.at(0), _kinksY.at(1));
        const bool hasHit =  state.getStateHasHit();
        //Must check we have a hit before we cut on residual. Will not try and look if no hit required.
        if(hasHit and _hit){
            streamlog_out(DEBUG5)<<"State has hit... "<< std::endl;
            streamlog_out(DEBUG5)<<"State has hit... " << " Res: " << state.getRes().at(0) <<"  "<<state.getRes().at(1)  << std::endl;
            TMatrixD cov = state.getHit().getCov();
            TMatrixD covState = state.getCov();
            const bool inResRange = compareXYToCut( state.getRes().at(0),state.getRes().at(1), _residualX.at(0), _residualX.at(1),_residualY.at(0), _residualY.at(1));
            streamlog_out(DEBUG5)<<"Residual Error... "<< std::endl;
            const bool inResErrorRange = compareXYToCut( sqrt(cov[0][0]),sqrt(cov[1][1]), _residualErrorX.at(0),_residualErrorX.at(1),_residualErrorY.at(0), _residualErrorY.at(1));
            streamlog_out(DEBUG5)<<"Cov pos Error... "<< std::endl;
            const bool inCovPosErrorRange = compareXYToCut(sqrt( covState[4][4]),sqrt(covState[3][3]), _covPosX.at(0),_covPosX.at(1),_covPosY.at(0), _covPosY.at(1));
//            std::cout<<"Cov pos Error, pass? "<< inCovPosErrorRange << " "<< sqrt(covState[3][3])<<" "<< sqrt(covState[4][4]) <<std::endl;
            streamlog_out(DEBUG5)<<"Cov slope Error... "<< std::endl;
            const bool inCovSlopeErrorRange = compareXYToCut(sqrt( covState[4][4]),sqrt(covState[3][3]), _covSlopeX.at(0),_covSlopeX.at(1),_covSlopeY.at(0), _covSlopeY.at(1));
            //Can not pass weightsX ?????
//            streamlog_out(DEBUG5)<<"Weights size: " << _weightsX.size() << " " << _weightsY.size()  << std::endl;
//            streamlog_out(DEBUG5)<<"Weights cuts: " << _weightsX.at(0) << " " << _weightsX.at(1) << " " << _weightsY.at(0) << "  " <<_weightsY.at(1) << std::endl;
//            streamlog_out(DEBUG5)<<"Weights: " << state.getHit().getWeight().at(0)<< " " << state.getHit().getWeight().at(0)<<  std::endl;

//            const bool inWeightRange =compareXYToCut( state.getHit().getWeight().at(0),state.getHit().getWeight().at(1), _weightsX.at(0), _weightsX.at(1), _weightsY.at(0), _weightsY.at(1) );
            if(inPosRange and inAngleRange and inKinksRange and hasHit and inResRange and inResErrorRange and inCovPosErrorRange and inCovSlopeErrorRange){
                return true;
            }else{
                return false;
            }
        }else{
            if(inPosRange and inAngleRange and inKinksRange){
                return true;
            }else{
                return false;
            }
        }
}
