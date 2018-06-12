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
//#include "EUTelTripletGBLDUTscatInstance.h"

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
    
  marlin::AIDAProcessor::tree(parent)->mkdir("GBLUtility");

  sixkxHisto = AIDAProcessor::histogramFactory(parent)->
    createHistogram1D( "GBLUtility/sixkx", 100, -10, 10 );
  sixkxHisto->setTitle( "kink x;kink x [mrad];triplet-driplet pairs" );

  sixkyHisto = AIDAProcessor::histogramFactory(parent)->
    createHistogram1D( "GBLUtility/sixky", 100, -10, 10 );
  sixkyHisto->setTitle( "kink y;kink y [mrad];triplet-driplet pairs" );

  sixdxHisto = AIDAProcessor::histogramFactory(parent)->
    createHistogram1D( "GBLUtility/sixdx", 100, -1, 1 );
  sixdxHisto->setTitle( "six match x;match x [mm];triplet-driplet pairs" );

  sixdyHisto = AIDAProcessor::histogramFactory(parent)->
    createHistogram1D( "GBLUtility/sixdy", 100, -1, 1 );
  sixdyHisto->setTitle( "six match y;match y [mm];triplet-driplet pairs" );

  sixdxcHisto = AIDAProcessor::histogramFactory(parent)->
    createHistogram1D( "GBLUtility/sixdxc", 100, -250, 250 );
  sixdxcHisto->setTitle( "six match x;track #Deltax[#mum];triplet-driplet pairs" );

  sixdycHisto = AIDAProcessor::histogramFactory(parent)->
    createHistogram1D( "GBLUtility/sixdyc", 100, -250, 250 );
  sixdycHisto->setTitle( "six match y;track #Deltay[#mum];triplet-driplet pairs" );

  sixkxcHisto = AIDAProcessor::histogramFactory(parent)->
    createHistogram1D( "GBLUtility/sixkxc", 100, -10, 10 );
  sixkxcHisto->setTitle( "kink x, x-y matched;kink x [mrad];triplet-driplet pairs" );

  sixkycHisto = AIDAProcessor::histogramFactory(parent)->
    createHistogram1D( "GBLUtility/sixkyc", 100, -10, 10 );
  sixkycHisto->setTitle( "kink y, x-y matched;kink y [mrad];triplet-driplet pairs" );

  sixxHisto = AIDAProcessor::histogramFactory(parent)->
    createHistogram1D( "GBLUtility/sixx", 240, -12, 12 );
  sixxHisto->setTitle( "six x at DUT;six x_{out} at DUT [mm];six-plane tracks" );

  sixyHisto = AIDAProcessor::histogramFactory(parent)->
    createHistogram1D( "GBLUtility/sixy", 120, -6, 6 );
  sixyHisto->setTitle( "six y at DUT;six y_{up} at DUT [mm];six-plane tracks" );

  sixxyHisto = AIDAProcessor::histogramFactory(parent)->
    createHistogram2D( "GBLUtility/sixxy", 240, -12, 12, 120, -6, 6 );
  sixxyHisto->setTitle( "six at DUT;six x_{out} at DUT [mm];six y_{up} at DUT [mm];six-plane tracks" );

  sixxycHisto = AIDAProcessor::histogramFactory(parent)->
    createHistogram2D( "GBLUtility/sixxyc", 240, -12, 12, 120, -6, 6 );
  sixxycHisto->setTitle( "six large kink;six x_{out} at DUT [mm];six y_{up} at DUT [mm];large kink tracks" );

  /*sixxylkHisto = AIDAProcessor::histogramFactory(parent)->
    createHistogram2D( "GBLUtilityTracks/sixxylk", 240, -12, 12, 120, -6, 6 );
  sixxylkHisto->setTitle( "six with REF link at DUT;six x_{out} at DUT [mm];six y_{up} at DUT [mm];six-plane tracks with REF link" );
  */

  kinkxvsxy = AIDAProcessor::histogramFactory(parent)->
    createProfile2D( "GBLUtility/kinkxvsxy", 120, -12, 12, 60, -6, 6, 0, 100 );
  kinkxvsxy->setTitle( "kink x;six x_{out} at DUT [mm];six y_{up} at DUT [mm];sqrt(<kink^{2}>) [mrad]" );

  kinkyvsxy = AIDAProcessor::histogramFactory(parent)->
    createProfile2D( "GBLUtility/kinkyvsxy", 120, -12, 12, 60, -6, 6, 0, 100 );
  kinkyvsxy->setTitle( "kink y;six x_{out} at DUT [mm];six y_{up} at DUT [mm];sqrt(<kink^{2}>) [mrad]" );

  kinkxyvsxy = AIDAProcessor::histogramFactory(parent)->
    createProfile2D( "GBLUtility/kinkxyvsxy", 120, -12, 12, 60, -6, 6, 0, 100 );
  kinkxyvsxy->setTitle( "kink <|x|+|y|>;six x_{out} at DUT [mm];six y_{up} at DUT [mm]; (|x| + |y|)/2 [mrad]" );

  kinkx = AIDAProcessor::histogramFactory(parent)->
    createHistogram1D( "GBLUtility/kinkx", 500, -5, 5 );
  kinkx->setTitle( "triplet kink x angle at DUT;x angle at DUT [mrad];tracks" );

  kinky = AIDAProcessor::histogramFactory(parent)->
    createHistogram1D( "GBLUtility/kinky", 500, -5, 5 );
  kinky->setTitle( "triplet kink y angle at DUT;y angle at DUT [mrad];tracks" );

  kinkxy = AIDAProcessor::histogramFactory(parent)->
    createHistogram1D( "GBLUtility/kinkxy", 500, -5, 5 );
  kinkxy->setTitle( "triplet kink xy angle at DUT;xy angle at DUT [mrad];tracks" );

  // for both triplet and driplet
  triddaMindutHisto = AIDAProcessor::histogramFactory(parent)->
    createHistogram1D( "triddaMindut", 1400, -2, 5 );
  triddaMindutHisto->setTitle( "minimal triplet distance at DUT;triplet distance at DUT [mm];telescope triplets" );

}

void EUTelTripletGBLUtility::MatchTriplets(std::vector<triplet> const & up, std::vector<EUTelTripletGBLUtility::triplet> const & down, double z_match, double trip_matching_cut, std::vector<EUTelTripletGBLUtility::track> &tracks) {

  // Cut on the matching of two triplets [mm]
  //double intersect_residual_cut = 0.1;

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

      // Track kinks as difference in triplet slopes:
      double kx = newtrack.kink_x();
      double ky = newtrack.kink_y();

      // driplet - triplet
      double dx = xB - xA; 
      double dy = yB - yA;

      sixkxHisto->fill( kx*1E3 );
      sixkyHisto->fill( ky*1E3 );
      sixdxHisto->fill( dx );
      sixdyHisto->fill( dy );

      if( abs(dy) < 0.5 ) sixdxcHisto->fill( dx*1E3 );
      if( abs(dx) < 0.5 ) sixdycHisto->fill( dy*1E3 );
      

      // match driplet and triplet:
      streamlog_out(DEBUG4) << "  Distance for matching x: " << fabs(dx)<< std::endl;
      streamlog_out(DEBUG4) << "  Distance for matching y: " << fabs(dy)<< std::endl;
      if( fabs(dx) > trip_matching_cut) continue;
      if( fabs(dy) > trip_matching_cut) continue;
      streamlog_out(DEBUG4) << "  Survived matching " << std::endl;

      // check isolation
      if( !IsolatedTrip || !IsolatedDrip ) {
	//hIso->fill(0);
	//std::cout << " me so non-isolated" << std::endl;
	continue;
      }
      streamlog_out(DEBUG4) << "  Trip and Drip isolated " << std::endl;
      //else hIso->fill(1);
      

      sixkxcHisto->fill( kx*1E3 );
      sixkycHisto->fill( ky*1E3 );
      sixxHisto->fill( -xA ); // -xA = x_DP = out
      sixyHisto->fill( -yA ); // -yA = y_DP = up
      sixxyHisto->fill( -xA, -yA ); // DP: x_out, y_up
      // Fill kink map histogram:
      if( abs( kx ) > 0.002 || abs( ky ) > 0.002 ) sixxycHisto->fill( -xA, -yA );      

      // apply fiducial cut
      if ( fabs(xA) >  9.0) continue;
      if (     -yA  < -4.0) continue;
      kinkx->fill( kx*1E3 ); //sqrt(<kink^2>) [mrad]
      kinky->fill( ky*1E3 ); //sqrt(<kink^2>) [mrad]
      kinkxy->fill( (fabs(kx)+fabs(ky))/2*1E3 ); // [mrad]
      kinkxvsxy->fill( -xA, -yA, fabs(kx)*1E3 ); //sqrt(<kink^2>) [mrad]
      kinkyvsxy->fill( -xA, -yA, fabs(ky)*1E3 ); //sqrt(<kink^2>) [mrad]
      kinkxyvsxy->fill( -xA, -yA, (fabs(kx) + fabs(ky))/2*1E3 ); // [mrad]

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

  triddaMindutHisto->fill(ddAMin);
  if(ddAMin < isolation_cut && ddAMin > -0.5) IsolatedTrip = false; // if there is only one triplet, ddAmin is still -1.

  return IsolatedTrip;
}

bool EUTelTripletGBLUtility::AttachDUT(EUTelTripletGBLUtility::triplet & triplet, std::vector<EUTelTripletGBLUtility::hit> const & hits, unsigned int dutID,  double dist_cut){

	auto zPos = geo::gGeometry().getPlaneZPosition(dutID);
	int minHitIx = -1;
	double minDist = std::numeric_limits<double>::max();
	double cut_squared = dist_cut*dist_cut;	

	auto trX = triplet.getx_at(zPos);
	auto trY = triplet.gety_at(zPos);

//	std::cout << "Triplet x/y: " << trX << "|" << trY << '\n'; 

	size_t ix = 0;	
	for(auto& hit: hits) {
		if(hit.plane == dutID) {
			auto hitX = hit.x;
			auto hitY = hit.y;
			auto dist = (trX-hitX)*(trX-hitX)+(trY-hitY)*(trY-hitY);
//			std::cout << "Hit x/y: " << hitX << "|" << hitY << " dist: " << dist << '\n'; 
			if(dist <= cut_squared && dist < minDist ){
				minHitIx = static_cast<int>(ix);
				//std::cout << "Dist: " << dist << std::endl;
			}
		}
		++ix;
	}

	if(minHitIx != -1) {
		/*
    auto hitX = hits[minHitIx].x;
    auto hitY = hits[minHitIx].y;
		auto dist = (trX-hitX)*(trX-hitX)+(trY-hitY)*(trY-hitY);
		std::cout << "Added hit with: " << dist << " on plane: " << hits[minHitIx].plane << '\n';
    */
		triplet.push_back_DUT(hits[minHitIx].plane, hits[minHitIx]);
		return true;
	}
return false;
}
/*
bool EUTelTripletGBLUtility::AttachDUT(std::vector<EUTelTripletGBLUtility::triplet> & triplets, std::vector<EUTelTripletGBLUtility::hit> const & hits, unsigned int dutID, double trip_res_cut, double trip_slope_cut){

	auto zPos = geo::gGeometry().getPlaneZPosition(dutID);

	struct distance {
		EUTelTripletGBLUtility::triplet * triplet;
		EUTelTripletGBLUtility::hit const * hit;		
		double dist;
		distance(EUTelTripletGBLUtility::triplet & triplet, EUTelTripletGBLUtility::hit const & hit, double const & dist): triplet(&triplet), hit(&hit), dist(dist){};
		distance(distance const & other): triplet(other.triplet), hit(other.hit), dist(other.dist){};
		distance & operator= (distance const & other){ triplet = other.triplet; hit=other.hit; dist=other.dist;};
	};

	std::vector<distance> distances;
	for(auto& triplet: triplets){
		auto trX = triplet.getx_at(zPos);
		auto trY = triplet.gety_at(zPos);
		for(auto& hit: hits){
			if(hit.plane == dutID) {
				auto hitX = hit.x;
				auto hitY = hit.x;
				auto dist = (trX-hitX)*(trX-hitX)+(trY-hitY)*(trY-hitY);
				distances.emplace_back(triplet, hit, dist);
				std::cout << "Distance: " << dist << '\n';
			}
		}
	}

	std::sort(distances.begin(), distances.end(), [](distance const & a, distance const & b) -> bool {return a.dist > b.dist;});

	while(!distances.empty()) {
		auto curTriplet = distances.front().triplet;
		auto curHit = distances.front().hit;
		//curTriplet->addDUTHit(curHit);

		for(auto it = distances.begin(); it != distances.end(); ) { 
			auto otherTriplet = it->triplet; 
			auto otherHit  = it->hit;
			if( ( curTriplet == otherTriplet ) || ( curHit == otherHit) ) {
				it = distances.erase(it);
			}  else {
				++it;
			}
		}
	}
}
*/
double EUTelTripletGBLUtility::PlaneEfficiency(std::vector<EUTelTripletGBLUtility::triplet> &eff_triplets_UP, std::vector<EUTelTripletGBLUtility::triplet> &eff_triplets_DOWN, std::vector<EUTelTripletGBLUtility::hit> &hits, unsigned int PUT, double track_match_z, double DUTz, double track_match_cut, double eff_radius, std::vector<AIDA::IProfile1D*> &profile) {

  std::vector<EUTelTripletGBLUtility::triplet> eff_triplets;

  // Iterate over all found eff-triplets to match them to the DUT (plane3):
  //int n_matched_trips = 0;
  //int n_unmatched_trips = 0;

  //std::cout << " n eff triplets UP   = " << eff_triplets_UP->size() << std::endl;
  //std::cout << " n eff triplets DOWN = " << eff_triplets_DOWN->size() << std::endl;

  for( auto& trip: eff_triplets_UP ) {

    // Track impact position at Matching Point from Upstream:
    double xA = trip.getx_at(track_match_z); // triplet impact point at matching position
    double yA = trip.gety_at(track_match_z);

    // check if trip is isolated
    bool IsolatedTrip = IsTripletIsolated(trip, eff_triplets_UP, track_match_z, track_match_cut*2.0001);

    for( auto& drip: eff_triplets_DOWN ){

      // Track impact position at Matching Point from Downstream:
      double xB = drip.getx_at(track_match_z); // triplet impact point at matching position
      double yB = drip.gety_at(track_match_z);

      // check if drip is isolated
      bool IsolatedDrip = IsTripletIsolated(drip, eff_triplets_DOWN, track_match_z, track_match_cut*2.0001);

      // driplet - triplet
      double dx = xB - xA; 
      double dy = yB - yA;

      // match driplet and triplet:
      if( fabs(dx) > track_match_cut) continue;
      if( fabs(dy) > track_match_cut) continue;

      //std::cout << " intersec ";

      // check isolation
      if( !IsolatedTrip || !IsolatedDrip ) continue;
      //std::cout << " , isolated ";

      // apply fiducial cut
      if ( fabs(xA) >  9.0) continue;
      if (     -yA  < -4.0) continue;

      //std::cout << " , fiducial " << std::endl;

      eff_triplets.emplace_back(trip);


    } // Downstream
  } // Upstream

  //std::cout << " n eff triplets = " << eff_triplets.size() << std::endl;

  int n_sum = 0;
  int n_matched = 0;
  for( auto& trip: eff_triplets ) {

    double ddAMin = -1.0;
    // extrapolate triplet to plane  under test
    double xTrip = trip.getx_at(DUTz);
    double yTrip = trip.gety_at(DUTz);

    for( auto& lhit: hits ){

      // Fill residuals of triplet and hit in the selected plane:
      if( lhit.plane == PUT ) {
	double xHit = lhit.x;
	double yHit = lhit.y;

	double ddA = sqrt( fabs(xHit - xTrip)*fabs(xHit - xTrip) 
	    + fabs(yHit - yTrip)*fabs(yHit - yTrip) );
	if(ddAMin < 0 || ddA < ddAMin) ddAMin = ddA;
      }
    } // end loop over hits

    // if distance is smaller then limit, accept this as matched Hit
    if(fabs(ddAMin) < eff_radius) {
      //n_matched_trips++;
      profile.at(0)->fill(-xTrip, 1.);
      profile.at(1)->fill(-yTrip, 1.);
      n_matched++;
      n_sum++;
    } else {
      profile.at(0)->fill(-xTrip, 0.);
      profile.at(1)->fill(-yTrip, 0.);
      n_sum++;
    }

  }
  
  eff_triplets_UP.clear();
  eff_triplets_DOWN.clear();

  double eff = static_cast<double>(n_matched)/static_cast<double>(n_sum);
  return eff;

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

EUTelTripletGBLUtility::hit const & EUTelTripletGBLUtility::track::gethit(int plane) {
  if(plane < 3) return upstream.gethit(plane);
  else return downstream.gethit(plane);
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
#endif
