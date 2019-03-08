#ifndef EUTelTripletGBLUtility_tcc
#define EUTelTripletGBLUtility_tcc
namespace eutelescope {

template<typename T>
void EUTelTripletGBLUtility::FindTriplets(std::vector<EUTelTripletGBLUtility::hit> const & hits, T const & triplet_sensor_ids, double trip_res_cut, double slope_cut, std::vector<EUTelTripletGBLUtility::triplet> & found_triplets, bool only_best_triplet, bool upstream) {

  if(triplet_sensor_ids.size() != 3){
    throw std::runtime_error("EUTelTripletGBLUtility::FindTriplets called with an invalid set of triplet_sensor_ids (size should be three entries)");
  }

  auto plane0 = static_cast<unsigned>(triplet_sensor_ids[0]);
  auto plane1 = static_cast<unsigned>(triplet_sensor_ids[1]);
  auto plane2 = static_cast<unsigned>(triplet_sensor_ids[2]);

  // get all hit is plane = plane0
  for( auto& ihit: hits ){
    if( ihit.plane != plane0 ) continue; // First plane

    // get all hit is plane = plane2
    for( auto& jhit: hits ){
      if( jhit.plane != plane2 ) continue; // Last plane

      double sum_res_old = -1.;
      // get all hit is plane = plane1
      for( auto& khit: hits ){
	if( khit.plane != plane1 ) continue; // Middle plane

	// Create new preliminary triplet from the three hits:
	EUTelTripletGBLUtility::triplet new_triplet(ihit,khit,jhit);

    //Create triplet slope plots
    if(upstream == 1){
		upstreamTripletSlopeX->fill(new_triplet.getdx()*1E3/new_triplet.getdz()); //factor 1E3 to convert from rad to mrad. To be checked
		upstreamTripletSlopeY->fill(new_triplet.getdy()*1E3/new_triplet.getdz());
	} else {
		downstreamTripletSlopeX->fill(new_triplet.getdx()*1E3/new_triplet.getdz());
		downstreamTripletSlopeY->fill(new_triplet.getdy()*1E3/new_triplet.getdz());
	}
	// Setting cuts on the triplet track angle:
	if( fabs(new_triplet.getdx()) > slope_cut * new_triplet.getdz()) continue;
	if( fabs(new_triplet.getdy()) > slope_cut * new_triplet.getdz()) continue;
    
    //Create triplet residual plots
    if(upstream == 1){
		upstreamTripletResidualX->fill(new_triplet.getdx(plane1));
		upstreamTripletResidualY->fill(new_triplet.getdy(plane1));
	} else {
		downstreamTripletResidualX->fill(new_triplet.getdx(plane1));
		downstreamTripletResidualY->fill(new_triplet.getdy(plane1));
	}
	// Setting cuts on the triplet residual on the middle plane
	if( fabs(new_triplet.getdx(plane1)) > trip_res_cut) continue;
	if( fabs(new_triplet.getdy(plane1)) > trip_res_cut) continue;

    // Edo: This (best triplets) is hardcoded as false in EUTelAlignGBL. Is this really useful?
    if(only_best_triplet) {
		// For low threshold (high noise) and/or high occupancy, use only the triplet with the smallest sum of residuals on plane1
		double sum_res = sqrt(new_triplet.getdx(plane1)*new_triplet.getdx(plane1) + new_triplet.getdy(plane1)*new_triplet.getdy(plane1));
		if(sum_res < sum_res_old){
	
		  // Remove the last one since it fits worse, not if its the first
		  found_triplets.pop_back();
		  // The triplet is accepted, push it back:
		  found_triplets.emplace_back(new_triplet);
		  streamlog_out(DEBUG2) << new_triplet;
		  sum_res_old = sum_res;
		}

		// update sum_res_old on first iteration
		if(sum_res_old < 0.) {
		  // The triplet is accepted, push it back:
		  found_triplets.emplace_back(new_triplet);
		  streamlog_out(DEBUG2) << new_triplet;
		  sum_res_old = sum_res;
		}
	} else {	
		found_triplets.emplace_back(new_triplet);
	}
      }//loop over hits
    }//loop over hits
  }// loop over hits

  //return triplets;
}

}//namespace
#endif
