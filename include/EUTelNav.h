#ifndef EUTELNAV_H
#define EUTELNAV_H

#include "EUTelGeometryTelescopeGeoDescription.h"
#include "TMatrix.h"
#include "TVector3.h"
#include "gear/BField.h"

namespace eutelescope 
{

class EUTelNav
{
	public: 
		static TMatrix getPropagationJacobianF( float x0, float y0, float z0, float px, float py, float pz, float beamQ, float dz);
		static TMatrixD getLocalToCurvilinearTransformMatrix(TVector3 globalMomentum, int  planeID, float charge);
		static TMatrixD getPropagationJacobianCurvilinear(float ds, float qbyp, TVector3 t1w, TVector3 t2w);
		static TVector3 getXYZfromArcLength(TVector3 pos, TVector3 pVec, float beamQ, double s);
		static TVector3 getXYZMomentumfromArcLength(TVector3 momentum, TVector3 globalPositionStart, float charge, float arcLength);
  	
	private:
		EUTelNav();
};

}
#endif
