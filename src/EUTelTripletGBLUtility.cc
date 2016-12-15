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
//#include "EUTelTripletGBLDUTscatInstance.h"

#include "EUTELESCOPE.h"
#include <cmath>

using namespace std;
using namespace eutelescope;


EUTelTripletGBLUtility::EUTelTripletGBLUtility(){}

TMatrixD EUTelTripletGBLUtility::JacobianPointToPoint( double ds ) {
  /* for GBL:
     Jacobian for straight line track
     track = q/p, x', y', x, y
     0,   1,  2,  3, 4
     */
  TMatrixD jac( 5, 5 );
  jac.UnitMatrix();
  jac[3][1] = ds; // x = x0 + xp * ds
  jac[4][2] = ds; // y = y0 + yp * ds
  return jac;
}

void EUTelTripletGBLUtility::MatchTriplets(std::vector<triplet> &up, std::vector<EUTelTripletGBLUtility::triplet> &down, double z_match, double trip_matching_cut, std::vector<EUTelTripletGBLUtility::track> &tracks) {

  // Cut on the matching of two triplets [mm]
  //double intersect_residual_cut = 0.1;

  for( std::vector<EUTelTripletGBLUtility::triplet>::iterator trip = up.begin(); trip != up.end(); trip++ ){

    // Track impact position at Matching Point from Upstream:
    double xA = (*trip).getx_at(z_match); // triplet impact point at matching position
    double yA = (*trip).gety_at(z_match);

    // check if trip is isolated. use at least double the trip_machting_cut for isolation in order to avoid double matching
    bool IsolatedTrip = IsTripletIsolated(trip, up, z_match, trip_matching_cut*2.0001);

    for( std::vector<EUTelTripletGBLUtility::triplet>::iterator drip = down.begin(); drip != down.end(); drip++ ){

      // Track impact position at Matching Point from Downstream:
      double xB = (*drip).getx_at(z_match); // triplet impact point at matching position
      double yB = (*drip).gety_at(z_match);

      // check if drip is isolated
      bool IsolatedDrip = IsTripletIsolated(drip, down, z_match, trip_matching_cut*2.0001);


      // Build a track candidate from one upstream and one downstream triplet:
      EUTelTripletGBLUtility::track newtrack((*trip),(*drip));

      // Track kinks as difference in triplet slopes:
      //double kx = newtrack.kink_x();
      //double ky = newtrack.kink_y();

      // driplet - triplet
      double dx = xB - xA; 
      double dy = yB - yA;

      /*sixkxHisto->fill( kx*1E3 );
      sixkyHisto->fill( ky*1E3 );
      sixdxHisto->fill( dx );
      sixdyHisto->fill( dy );


      if( abs(dy) < 0.4 ) sixdxcHisto->fill( dx*1E3 ); // sig = 17 um at 5 GeV
      if( abs(dx) < 0.4 ) sixdycHisto->fill( dy*1E3 );
      */

      // match driplet and triplet:
      if( fabs(dx) > trip_matching_cut) continue;
      if( fabs(dy) > trip_matching_cut) continue;

      // check isolation
      //if( !IsolatedTrip || !IsolatedDrip ) {
	//hIso->fill(0);
	//continue;
      //}
      //else hIso->fill(1);
      

      /*sixkxcHisto->fill( kx*1E3 );
      sixkycHisto->fill( ky*1E3 );
      sixxHisto->fill( -xA ); // -xA = x_DP = out
      sixyHisto->fill( -yA ); // -yA = y_DP = up
      sixxyHisto->fill( -xA, -yA ); // DP: x_out, y_up
      */

      // Fill kink map histogram:
      /*if( abs( kx ) > 0.002 || abs( ky ) > 0.002 ) sixxycHisto->fill( -xA, -yA );
       */

      //kinkvsxy->fill( -xA, -yA, (kx*kx + ky*ky)*1E6 ); //<kink^2> [mrad^2]

      // apply fiducial cut
      if ( fabs(xA) >  9.0) continue;
      if (     -yA  < -4.0) continue;

      // Add the track to the vector if trip/drip are isolated, the triplets are matched, and all other cuts are passed
      tracks.push_back(newtrack);

    } // Downstream
  } // Upstream

  streamlog_out(DEBUG2) << "Found " << tracks.size() << " tracks from matched triplets." << std::endl;
  //return tracks;
}

bool EUTelTripletGBLUtility::IsTripletIsolated(std::vector<EUTelTripletGBLUtility::triplet>::iterator it, std::vector<EUTelTripletGBLUtility::triplet> &trip, double z_match, double isolation_cut) { // isolation_cut is defaulted to 0.3 mm
  bool IsolatedTrip = true;

  // check first if trip is isolated
  double xA = (*it).getx_at(z_match); // triplet impact point at matching position
  double yA = (*it).gety_at(z_match);


  double ddAMin = -1.0;
  for( std::vector<EUTelTripletGBLUtility::triplet>::iterator tripIsoCheck = trip.begin(); tripIsoCheck != trip.end(); tripIsoCheck++ ) {
    if(it != tripIsoCheck){
      double xAIsoCheck = (*tripIsoCheck).getx_at(z_match);
      double yAIsoCheck = (*tripIsoCheck).gety_at(z_match);
      double ddA = sqrt( fabs(xAIsoCheck - xA)*fabs(xAIsoCheck - xA) 
	  + fabs(yAIsoCheck - yA)*fabs(yAIsoCheck - yA) );
      if(ddAMin < 0 || ddA < ddAMin) ddAMin = ddA;
    }
  }

  //triddaMindutHisto->fill(ddAMin);
  if(ddAMin < isolation_cut) IsolatedTrip = false;

  return IsolatedTrip;
}

void EUTelTripletGBLUtility::FindTriplets(std::vector<EUTelTripletGBLUtility::hit> &hits, unsigned int plane0, unsigned int plane1, unsigned int plane2, double trip_res_cut, double slope_cut, std::vector<EUTelTripletGBLUtility::triplet> &triplets) {

  // get all hit is plane = plane0
  for( std::vector<EUTelTripletGBLUtility::hit>::iterator ihit = hits.begin(); ihit != hits.end(); ihit++ ){
    if( (*ihit).plane != plane0 ) continue; // First plane

    // get all hit is plane = plane2
    for( std::vector<EUTelTripletGBLUtility::hit>::iterator jhit = hits.begin(); jhit != hits.end(); jhit++ ){
      if( (*jhit).plane != plane2 ) continue; // Last plane

      // get all hit is plane = plane1
      for( std::vector<EUTelTripletGBLUtility::hit>::iterator khit = hits.begin(); khit != hits.end(); khit++ ){
	if( (*khit).plane != plane1 ) continue; // Middle plane

	// Create new preliminary triplet from the three hits:
	EUTelTripletGBLUtility::triplet new_triplet((*ihit),(*khit),(*jhit));

	// Setting cuts on the triplet track angle:
	if( fabs(new_triplet.getdx()) > slope_cut * new_triplet.getdz()) continue;
	if( fabs(new_triplet.getdy()) > slope_cut * new_triplet.getdz()) continue;

	// Setting cuts on the triplet residual on the middle plane
	if( fabs(new_triplet.getdx(plane1)) > trip_res_cut) continue;
	if( fabs(new_triplet.getdy(plane1)) > trip_res_cut) continue;

	/*
	// For low threshold (high noise) and/or high occupancy, use only the triplet with the smallest sum of residuals on plane1
	sum_res = abs(new_triplet.getdx(plane1)) + abs(new_triplet.getdy(plane1));
	if(sum_res < sum_res_old || IsFirst){

	// Remove the last one since it fits worse, not if its the first
	if(!IsFirst) triplets->pop_back();
	// The triplet is accepted, push it back:
	triplets->push_back(new_triplet);
	IsFirst = false;
	streamlog_out(DEBUG2) << new_triplet;
	sum_res_old = sum_res;
	}
	*/

	// The triplet is accepted, push it back:
	triplets.push_back(new_triplet);
	streamlog_out(DEBUG2) << new_triplet;

      }//loop over hits
    }//loop over hits
  }// loop over hits

  //return triplets;
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

EUTelTripletGBLUtility::triplet EUTelTripletGBLUtility::track::get_upstream() {
  return upstream;
}

EUTelTripletGBLUtility::triplet EUTelTripletGBLUtility::track::get_downstream() {
  return downstream;
}

EUTelTripletGBLUtility::hit EUTelTripletGBLUtility::track::gethit(int plane) {
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

EUTelTripletGBLUtility::hit EUTelTripletGBLUtility::triplet::getpoint_at(double z) {
  hit impact;
  impact.z = z - base().z;
  impact.x = base().x + slope().x * impact.z;
  impact.y = base().y + slope().y * impact.z;
  return impact;
}

double EUTelTripletGBLUtility::triplet::getx_at(double z) {
  return this->base().x + this->slope().x * (z - this->base().z);
}

double EUTelTripletGBLUtility::triplet::getdx() {
  return this->hits.rbegin()->second.x - this->hits.begin()->second.x;
}

double EUTelTripletGBLUtility::triplet::getdx(int ipl) {
  return this->hits[ipl].x - this->base().x - this->slope().x * (this->hits[ipl].z - this->base().z);
}

double EUTelTripletGBLUtility::triplet::getdx(hit point) {
  return point.x - this->base().x - this->slope().x * (point.z - this->base().z);
}

double EUTelTripletGBLUtility::triplet::gety_at(double z) {
  return this->base().y + this->slope().y * (z - this->base().z);
}

double EUTelTripletGBLUtility::triplet::getdy() {
  return this->hits.rbegin()->second.y - this->hits.begin()->second.y;
}

double EUTelTripletGBLUtility::triplet::getdy(int ipl) {
  return this->hits[ipl].y - this->base().y - this->slope().y * (this->hits[ipl].z - this->base().z);
}

double EUTelTripletGBLUtility::triplet::getdy(hit point) {
  return point.y - this->base().y - this->slope().y * (point.z - this->base().z);
}

double EUTelTripletGBLUtility::triplet::getdz() {
  return this->hits.rbegin()->second.z - this->hits.begin()->second.z;
}

EUTelTripletGBLUtility::hit EUTelTripletGBLUtility::triplet::gethit(int plane) {
  return this->hits[plane];
}

EUTelTripletGBLUtility::hit EUTelTripletGBLUtility::triplet::base() {
  hit center;
  center.x = 0.5*( this->hits.begin()->second.x + this->hits.rbegin()->second.x );
  center.y = 0.5*( this->hits.begin()->second.y + this->hits.rbegin()->second.y );
  center.z = 0.5*( this->hits.begin()->second.z + this->hits.rbegin()->second.z );
  return center;
}

EUTelTripletGBLUtility::hit EUTelTripletGBLUtility::triplet::slope() {
  hit sl;
  double dz = (this->hits.rbegin()->second.z - this->hits.begin()->second.z);
  sl.x = (this->hits.rbegin()->second.x - this->hits.begin()->second.x) / dz;
  sl.y = (this->hits.rbegin()->second.y - this->hits.begin()->second.y) / dz;
  return sl;
}

#endif
