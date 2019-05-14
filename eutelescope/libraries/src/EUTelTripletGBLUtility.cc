// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#if defined USE_GEAR

// EUTelescope includes:
#include "EUTelTripletGBLUtility.h"
#include "EUTelGeometryTelescopeGeoDescription.h"

#include "EUTELESCOPE.h"
#include <cmath>
#include <memory>
#include <type_traits>

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IAxis.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/IProfile2D.h>
#include <AIDA/ITree.h>
#endif

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/ProcessorMgr.h"
#include "marlin/Global.h"

using namespace std;
using namespace eutelescope;
using namespace marlin;


EUTelTripletGBLUtility::EUTelTripletGBLUtility(){}

Eigen::Matrix<double, 5,5> EUTelTripletGBLUtility::JacobianPointToPoint( double ds ) {
  /* for GBL:
     Jacobian for straight line track
     track = q/p, x', y', x, y
     0,   1,  2,  3, 4
     */
  Eigen::Matrix<double, 5,5> jac = Eigen::Matrix<double, 5,5>::Identity();
  jac(3,1) = ds; // x = x0 + xp * ds
  jac(4,2) = ds; // y = y0 + yp * ds
  return jac;
}

void EUTelTripletGBLUtility::bookHistos(){
     
  //cut plots
  marlin::AIDAProcessor::tree(parent)->mkdir("Cuts");
  
  upstreamTripletSlopeX = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/upstreamTripletSlopeCutX", 1000, -3, 3 ); //binning to be reviewed
  upstreamTripletSlopeX->setTitle( "Upstream Triplet Slope X;Upstream Triplet Slope X [mrad];Counts" );
  
  upstreamTripletSlopeY = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/upstreamTripletSlopeCutY", 1000, -3, 3 ); //binning to be reviewed
  upstreamTripletSlopeY->setTitle( "Upstream Triplet Slope Y;Upstream Triplet Slope Y [mrad];Counts" );
  
  downstreamTripletSlopeX = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/downstreamTripletSlopeCutX", 1000, -3, 3 ); //binning to be reviewed
  downstreamTripletSlopeX->setTitle( "Downstream Triplet Slope X;Downstream Triplet Slope X [mrad];Counts" );
  
  downstreamTripletSlopeY = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/downstreamTripletSlopeCutY", 1000, -3, 3 ); //binning to be reviewed
  downstreamTripletSlopeY->setTitle( "Downstream Triplet Slope Y;Downstream Triplet Slope Y [mrad];Counts" );
  
  upstreamTripletResidualX = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/upstreamTripletResidualCutX", 1000, -3, 3 ); //binning to be reviewed
  upstreamTripletResidualX->setTitle( "Upstream Triplet Residual X;Upstream Triplet Residual X [mm];Counts" );
  
  upstreamTripletResidualY = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/upstreamTripletResidualCutY", 1000, -3, 3 ); //binning to be reviewed
  upstreamTripletResidualY->setTitle( "Upstream Triplet Residual Y;Upstream Triplet Residual Y [mm];Counts" );
  
  downstreamTripletResidualX = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/downstreamTripletResidualCutX", 1000, -3, 3 ); //binning to be reviewed
  downstreamTripletResidualX->setTitle( "Downstream Triplet Residual X;Downstream Triplet Residual X [mm];Counts" );
  
  downstreamTripletResidualY = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/downstreamTripletResidualCutY", 1000, -3, 3 ); //binning to be reviewed
  downstreamTripletResidualY->setTitle( "Downstream Triplet Residual Y;Downstream Triplet Residual Y [mm];Counts" );
  
  tripletMatchingResidualX = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/tripletMatchingResidualCutX", 1000, -3, 3 ); //binning to be reviewed
  tripletMatchingResidualX->setTitle( "Triplet Matching Residual X;Triplet Matching Residual X [mm];Counts" );
  
  tripletMatchingResidualY = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/tripletMatchingResidualCutY", 1000, -3, 3 ); //binning to be reviewed
  tripletMatchingResidualY->setTitle( "Triplet Matching Residual Y;Triplet Matching Residual Y [mm];Counts" );
  
  DUTMatchingResidualX = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/DUTMatchingResidualCutX", 1000, -3, 3 ); //binning to be reviewed
  DUTMatchingResidualX->setTitle( "DUT Matching Residual local X;DUT Matching Residual local X [mm];Counts" );
  
  DUTMatchingResidualY = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/DUTMatchingResidualCutY", 1000, -3, 3 ); //binning to be reviewed
  DUTMatchingResidualY->setTitle( "DUT Matching Residual local Y;DUT Matching Residual local Y [mm];Counts" );
  
  DUTHitNumber = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/DUTHitNumber", 21, -0.5, 20.5 ); //binning to be reviewed
  DUTHitNumber->setTitle( "Number of Hits matched to a track per DUT ID;DUT ID;Number of Hits matched to a track" );

}

void EUTelTripletGBLUtility::MatchTriplets(std::vector<triplet> const & up, std::vector<EUTelTripletGBLUtility::triplet> const & down, double z_match, double trip_matching_cut, std::vector<EUTelTripletGBLUtility::track> &tracks) {

  // Cut on the matching of two triplets [mm]

  for( auto& trip: up ){

    // Track impact position at Matching Point from Upstream:
    double xA = trip.getx_at(z_match); // triplet impact point at matching position
    double yA = trip.gety_at(z_match);

    // check if trip is isolated. use at least double the trip_machting_cut for isolation in order to avoid double matching
    bool IsolatedTrip = IsTripletIsolated(trip, up, z_match, trip_matching_cut*2.0001);
    streamlog_out(DEBUG4) << "  Is triplet isolated? " << IsolatedTrip << std::endl;

    for( auto& drip: down ){

      // Track impact position at Matching Point from Downstream:
      double xB = drip.getx_at(z_match); // triplet impact point at matching position
      double yB = drip.gety_at(z_match);

      // check if drip is isolated
      bool IsolatedDrip = IsTripletIsolated(drip, down, z_match, trip_matching_cut*2.0001);
      streamlog_out(DEBUG4) << "  Is driplet isolated? " << IsolatedDrip << std::endl;


      // Build a track candidate from one upstream and one downstream triplet:
      EUTelTripletGBLUtility::track newtrack(trip,drip);

      // driplet - triplet
      double dx = xB - xA; 
      double dy = yB - yA;

      
      //cut plots
      tripletMatchingResidualX->fill(dx);
      tripletMatchingResidualY->fill(dy);
      // match driplet and triplet:
      streamlog_out(DEBUG4) << "  Distance for matching x: " << fabs(dx)<< std::endl;
      streamlog_out(DEBUG4) << "  Distance for matching y: " << fabs(dy)<< std::endl;
      if( fabs(dx) > trip_matching_cut) continue;
      if( fabs(dy) > trip_matching_cut) continue;
      streamlog_out(DEBUG4) << "  Survived matching " << std::endl;

      // check isolation
      if( !IsolatedTrip || !IsolatedDrip ) {
	continue;
      }
      streamlog_out(DEBUG4) << "  Trip and Drip isolated " << std::endl;      

      // Add the track to the vector if trip/drip are isolated, the triplets are matched, and all other cuts are passed
      tracks.push_back(newtrack);

    } // Downstream
  } // Upstream

  streamlog_out(DEBUG2) << "Found " << tracks.size() << " tracks from matched triplets." << std::endl;
  //return tracks;
}

bool EUTelTripletGBLUtility::IsTripletIsolated(EUTelTripletGBLUtility::triplet const & it, std::vector<EUTelTripletGBLUtility::triplet> const & trip, double z_match, double isolation_cut) { // isolation_cut is defaulted to 0.3 mm
  bool IsolatedTrip = true;

  // check first if trip is isolated
  double xA = it.getx_at(z_match); // triplet impact point at matching position
  double yA = it.gety_at(z_match);

  double ddAMin = -1.0;
  for( auto& tripIsoCheck: trip ) {
    if( &it != &tripIsoCheck ){
      double xAIsoCheck = tripIsoCheck.getx_at(z_match);
      double yAIsoCheck = tripIsoCheck.gety_at(z_match);
      double ddA = sqrt( fabs(xAIsoCheck - xA)*fabs(xAIsoCheck - xA) 
	  + fabs(yAIsoCheck - yA)*fabs(yAIsoCheck - yA) );
      if(ddAMin < 0 || ddA < ddAMin) ddAMin = ddA;
	}
  }

  if(ddAMin < isolation_cut && ddAMin > -0.5) IsolatedTrip = false; // if there is only one triplet, ddAmin is still -1.

  return IsolatedTrip;
}

bool EUTelTripletGBLUtility::AttachDUT(EUTelTripletGBLUtility::triplet & triplet, std::vector<EUTelTripletGBLUtility::hit> const & hits, unsigned int dutID,  std::vector<float> dist_cuts){

	auto zPos = geo::gGeometry().getPlaneZPosition(dutID);
	int minHitIx = -1;
	double minDist = std::numeric_limits<float>::max();

	auto trX = triplet.getx_at(zPos);
	auto trY = triplet.gety_at(zPos);
	size_t ix = 0;	
	for(auto& hit: hits) {
		if(hit.plane == dutID) {
			auto hitX = hit.x;
			auto hitY = hit.y;
			auto distX = fabs(trX-hitX);
			auto distY = fabs(trY-hitY);
            double dist = distX*distX + distY*distY;
			DUTMatchingResidualX->fill(trX-hitX);
			DUTMatchingResidualY->fill(trY-hitY);
			if(distX <= dist_cuts.at(0) && distY <= dist_cuts.at(1) && dist < minDist ){
				minHitIx = static_cast<int>(ix);
				DUTHitNumber->fill(dutID);
				minDist = dist;
			}
		}
		++ix;
	}

	if(minHitIx != -1) {
		triplet.push_back_DUT(hits[minHitIx].plane, hits[minHitIx]);
		return true;
	}
return false;
}

EUTelTripletGBLUtility::track::track(triplet up, triplet down) : upstream(up), downstream(down) {}

double EUTelTripletGBLUtility::track::kink_x() {
  return (downstream.slope().x - upstream.slope().x);
}

double EUTelTripletGBLUtility::track::kink_y() {
  return (downstream.slope().y - upstream.slope().y);
}

EUTelTripletGBLUtility::hit EUTelTripletGBLUtility::track::intersect() {
  hit inter;
  // Re-check what this actually is...
  // and simplify using triplet class members...
  inter.x = ( upstream.base().x - upstream.slope().x * upstream.base().z - downstream.base().x + downstream.slope().x * downstream.base().z ) / kink_x();
  inter.y = ( upstream.base().y - upstream.slope().y * upstream.base().z - downstream.base().y + downstream.slope().y * downstream.base().z ) / kink_y();
  return inter;
}

EUTelTripletGBLUtility::triplet& EUTelTripletGBLUtility::track::get_upstream() {
  return upstream;
}

EUTelTripletGBLUtility::triplet& EUTelTripletGBLUtility::track::get_downstream() {
  return downstream;
}

EUTelTripletGBLUtility::triplet::triplet() : linked_dut(false), hits() {
  // Empty default constructor
}

EUTelTripletGBLUtility::triplet::triplet(hit hit0, hit hit1, hit hit2) : linked_dut(false), hits() {
  triplet();
  filltriplet(hit0, hit1, hit2);
}

EUTelTripletGBLUtility::hit EUTelTripletGBLUtility::triplet::getpoint_at(double z) const{
  hit impact;
  impact.z = z - base().z;
  impact.x = base().x + slope().x * impact.z;
  impact.y = base().y + slope().y * impact.z;
  return impact;
}

double EUTelTripletGBLUtility::triplet::getx_at(double z) const {
  return base().x + slope().x * (z - base().z);
}

double EUTelTripletGBLUtility::triplet::getdx() const {
  return hits.rbegin()->second.x - hits.begin()->second.x;
}

double EUTelTripletGBLUtility::triplet::getdx(int ipl) const {
  return hits.at(ipl).x - base().x - slope().x * (hits.at(ipl).z - base().z);
}

double EUTelTripletGBLUtility::triplet::getdx(hit point) const {
  return point.x - base().x - slope().x * (point.z - base().z);
}

double EUTelTripletGBLUtility::triplet::gety_at(double z) const {
  return base().y + slope().y * (z - base().z);
}

double EUTelTripletGBLUtility::triplet::getdy() const {
  return hits.rbegin()->second.y - hits.begin()->second.y;
}

double EUTelTripletGBLUtility::triplet::getdy(int ipl) const {
  return hits.at(ipl).y - base().y - slope().y * (hits.at(ipl).z - base().z);
}

double EUTelTripletGBLUtility::triplet::getdy(hit point) const {
  return point.y - base().y - slope().y * (point.z - base().z);
}

double EUTelTripletGBLUtility::triplet::getdz() const {
  return hits.rbegin()->second.z - hits.begin()->second.z;
}

EUTelTripletGBLUtility::hit const & EUTelTripletGBLUtility::triplet::gethit(int plane) const {
  return hits.at(plane);
}

EUTelTripletGBLUtility::hit EUTelTripletGBLUtility::triplet::base() const {
  hit center;
  center.x = 0.5*( hits.begin()->second.x + hits.rbegin()->second.x );
  center.y = 0.5*( hits.begin()->second.y + hits.rbegin()->second.y );
  center.z = 0.5*( hits.begin()->second.z + hits.rbegin()->second.z );
  return center;
}

EUTelTripletGBLUtility::hit EUTelTripletGBLUtility::triplet::slope() const {
  hit sl;
  double dz = (hits.rbegin()->second.z - hits.begin()->second.z);
  sl.x = (hits.rbegin()->second.x - hits.begin()->second.x) / dz;
  sl.y = (hits.rbegin()->second.y - hits.begin()->second.y) / dz;
  return sl;
}

std::pair<double,double> EUTelTripletGBLUtility::doIterativeGaussianFit(AIDA::IHistogram1D* in_hist, int rebinFactor) const {
    
    //--- First convert from IHistogram to TH1 so that we could use ROOT's fitters
    int nb_bins = in_hist->axis().bins();
    streamlog_out (DEBUG5) << in_hist->title() << " has " << in_hist->allEntries() << " entries" << std::endl;
    
    TH1D current(in_hist->title().data(), in_hist->title().data(), nb_bins, in_hist->axis().lowerEdge(), in_hist->axis().upperEdge());
    for(int ibin = 0 ; ibin < nb_bins ; ibin++) current.SetBinContent(ibin+1, in_hist->binEntries(ibin));
    current.Rebin(rebinFactor);
    
    double startForMaxFraction = 0.2;
    double max = current.GetMaximum();
    int min_bin = current.FindFirstBinAbove(startForMaxFraction*max);
    int max_bin = current.FindLastBinAbove(startForMaxFraction*max);
    if(min_bin < 0.2*current.GetNbinsX()) {
        startForMaxFraction = 0.35;
        min_bin = current.FindFirstBinAbove(startForMaxFraction*max);
        max_bin = current.FindLastBinAbove(startForMaxFraction*max);
    }
    if(min_bin < 0.2*current.GetNbinsX()) {
        startForMaxFraction = 0.5;
        min_bin = current.FindFirstBinAbove(startForMaxFraction*max);
        max_bin = current.FindLastBinAbove(startForMaxFraction*max);
    }
    for(int ibin = min_bin ; ibin < max_bin-1 ; ibin++) {
        if( (current.GetBinContent(ibin) > 0.7*current.GetBinContent(ibin+1)) && (current.GetBinContent(ibin+2) > 0.7*current.GetBinContent(ibin+1)) ) {
            current.Smooth(2);
            min_bin = current.FindFirstBinAbove(startForMaxFraction*max);
            max_bin = current.FindLastBinAbove(startForMaxFraction*max);
        }
    }
    
    double bound_low = current.GetBinCenter(min_bin);
    double bound_hi = current.GetBinCenter(max_bin);
    double mean = (bound_hi+bound_low)/2.;
    double sigma = bound_hi-bound_low;
    
    TF1 gausfit("gausfit", "([3]+[0]*exp(-0.5*( ((x-[1])/[2])*((x-[1])/[2]) ) ))", bound_low, bound_hi);
    
    gausfit.SetParameter(0, max); //-current.GetBinContent(1));
    gausfit.SetParameter(1, mean);
    gausfit.SetParameter(2, sigma);
    gausfit.SetParameter(3, current.GetBinContent(1));
    
    gausfit.SetParLimits(0, 0.1*max, 2*max);
    gausfit.SetParLimits(1, bound_low, bound_hi);
    gausfit.SetParLimits(2, 0.1*(bound_hi-bound_low), 1.5*(bound_hi-bound_low));
    gausfit.SetParLimits(3, 0, 0.2*max);
             
    TFitResultPtr fitresult = current.Fit(&gausfit, "SQ", "", bound_low, bound_hi);
    //FIXME: why refit here with same settings? brute force?
    if(!fitresult->IsValid()) fitresult = current.Fit(&gausfit, "SQ", "", bound_low, bound_hi);
        
    mean = fitresult->GetParams()[1];
    sigma = fitresult->GetParams()[2];

    /*--- This is for DEBUGGING ONLY !!!   
    TCanvas can; 
    can.cd();
    current.Draw();
    can.SaveAs(TString(in_hist->title()).ReplaceAll(" ","")+".pdf");*/    
    
    return std::pair<double,double> (mean, sigma);
}

void EUTelTripletGBLUtility::determineBestCuts() const {

    std::pair<double, double> res_upstreamTripletSlopeX      = doIterativeGaussianFit(upstreamTripletSlopeX, 10);
    std::pair<double, double> res_upstreamTripletSlopeY      = doIterativeGaussianFit(upstreamTripletSlopeY, 10);    
    std::pair<double, double> res_downstreamTripletSlopeX    = doIterativeGaussianFit(downstreamTripletSlopeX, 10);
    std::pair<double, double> res_downstreamTripletSlopeY    = doIterativeGaussianFit(downstreamTripletSlopeY, 10);
    std::pair<double, double> res_upstreamTripletResidualX   = doIterativeGaussianFit(upstreamTripletResidualX, 1);
    std::pair<double, double> res_upstreamTripletResidualY   = doIterativeGaussianFit(upstreamTripletResidualY, 1);
    std::pair<double, double> res_downstreamTripletResidualX = doIterativeGaussianFit(downstreamTripletResidualX, 1);
    std::pair<double, double> res_downstreamTripletResidualY = doIterativeGaussianFit(downstreamTripletResidualY, 1);
    std::pair<double, double> res_tripletMatchingResidualX   = doIterativeGaussianFit(tripletMatchingResidualX, 1);
    std::pair<double, double> res_tripletMatchingResidualY   = doIterativeGaussianFit(tripletMatchingResidualY, 1);

    streamlog_out (DEBUG5) << "__________ *** An iterative gaussian fit determined the following parameters for the distributions we will cut on : __________" << std::endl;
    streamlog_out (DEBUG5) << "upstreamTripletResidualX --- Mean = " << res_upstreamTripletResidualX.first << " ; Sigma = " << res_upstreamTripletResidualX.second << std::endl;
    streamlog_out (DEBUG5) << "upstreamTripletSlopeX --- Mean = " << res_upstreamTripletSlopeX.first << " ; Sigma = " << res_upstreamTripletSlopeX.second << std::endl;
    streamlog_out (DEBUG5) << "upstreamTripletSlopeY --- Mean = " << res_upstreamTripletSlopeY.first << " ; Sigma = " << res_upstreamTripletSlopeY.second << std::endl;
    streamlog_out (DEBUG5) << "downstreamTripletSlopeX --- Mean = " << res_downstreamTripletSlopeX.first << " ; Sigma = " << res_downstreamTripletSlopeX.second << std::endl;
    streamlog_out (DEBUG5) << "downstreamTripletSlopeY --- Mean = " << res_downstreamTripletSlopeY.first << " ; Sigma = " << res_downstreamTripletSlopeY.second << std::endl;
    streamlog_out (DEBUG5) << "upstreamTripletResidualY --- Mean = " << res_upstreamTripletResidualY.first << " ; Sigma = " << res_upstreamTripletResidualY.second << std::endl;
    streamlog_out (DEBUG5) << "downstreamTripletResidualX --- Mean = " << res_downstreamTripletResidualX.first << " ; Sigma = " << res_downstreamTripletResidualX.second << std::endl;
    streamlog_out (DEBUG5) << "downstreamTripletResidualY --- Mean = " << res_downstreamTripletResidualY.first << " ; Sigma = " << res_downstreamTripletResidualY.second << std::endl;
    streamlog_out (DEBUG5) << "tripletMatchingResidualX --- Mean = " << res_tripletMatchingResidualX.first << " ; Sigma = " << res_tripletMatchingResidualX.second << std::endl;
    streamlog_out (DEBUG5) << "tripletMatchingResidualY --- Mean = " << res_tripletMatchingResidualY.first << " ; Sigma = " << res_tripletMatchingResidualY.second << std::endl;

    //calculate recommended cut values
    double nb_sigma = 4.;
    double upstreamTripletCut = std::max( fabs(res_upstreamTripletResidualX.first)+nb_sigma*res_upstreamTripletResidualX.second, 
					  fabs(res_upstreamTripletResidualY.first)+nb_sigma*res_upstreamTripletResidualY.second );
    double downstreamTripletCut = std::max( fabs(res_downstreamTripletResidualX.first)+nb_sigma*res_downstreamTripletResidualX.second,
					    fabs(res_downstreamTripletResidualY.first)+nb_sigma*res_downstreamTripletResidualY.second );
    double upstreamSlopeCut = std::max( fabs(res_upstreamTripletSlopeX.first)+nb_sigma*res_upstreamTripletSlopeX.second, 
					fabs(res_upstreamTripletSlopeY.first)+nb_sigma*res_upstreamTripletSlopeY.second );
    double downstreamSlopeCut = std::max( fabs(res_downstreamTripletSlopeX.first)+nb_sigma*res_downstreamTripletSlopeX.second,
					  fabs(res_downstreamTripletSlopeY.first)+nb_sigma*res_downstreamTripletSlopeY.second );
    double tripletMatchingCut = std::max( fabs(res_tripletMatchingResidualX.first)+nb_sigma*res_tripletMatchingResidualX.second, 
					  fabs(res_tripletMatchingResidualY.first)+nb_sigma*res_tripletMatchingResidualY.second );

    streamlog_out (MESSAGE8) << "We recommend the following cuts, although they should be checked ! :" << std::endl;
    streamlog_out (MESSAGE8) << "UpstreamTripletCut 	= " << upstreamTripletCut  << std::endl;
    streamlog_out (MESSAGE8) << "DownstreamTripletCut 	= " << downstreamTripletCut << std::endl;
    streamlog_out (MESSAGE8) << "UpstreamSlopeCut   	= " << upstreamSlopeCut << std::endl;
    streamlog_out (MESSAGE8) << "DownstreamSlopeCut 	= " << downstreamSlopeCut << std::endl;
    streamlog_out (MESSAGE8) << "TripletMatchingCut 	= " << tripletMatchingCut << std::endl;
    //we may have no DUT hits
    if(DUTMatchingResidualX->allEntries() > 0  && DUTMatchingResidualY->allEntries() > 0){
      std::pair<double, double> res_DUTMatchingResidualX = doIterativeGaussianFit(DUTMatchingResidualX, 1);
      std::pair<double, double> res_DUTMatchingResidualY = doIterativeGaussianFit(DUTMatchingResidualY, 1);
      streamlog_out (DEBUG5) << "DUTMatchingResidualX --- Mean = " << res_DUTMatchingResidualX.first << " ; Sigma = " << res_DUTMatchingResidualX.second << std::endl;
      streamlog_out (DEBUG5) << "DUTMatchingResidualY --- Mean = " << res_DUTMatchingResidualY.first << " ; Sigma = " << res_DUTMatchingResidualY.second << std::endl;
      double dutXCut = fabs(res_DUTMatchingResidualX.first)+nb_sigma*res_DUTMatchingResidualX.second;
      double dutYCut = fabs(res_DUTMatchingResidualY.first)+nb_sigma*res_DUTMatchingResidualY.second;
      streamlog_out (MESSAGE8) << "DUTCuts 		= " << dutXCut << " " << dutYCut << std::endl; 
    } else {
      streamlog_out (MESSAGE8) << "The DUT is propably not an active device - no hits detected" << std::endl;	
    }
}

#endif
