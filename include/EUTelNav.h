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
		static TMatrixD getLocalToCurvilinearTransformMatrixLimit(TVector3 globalMomentum, int  planeID, float charge);
		static TMatrixD getMeasToGlobal(TVector3 t1w, int  planeID);

		static TMatrixD getPropagationJacobianCurvilinear(float ds, float qbyp, TVector3 t1w, TVector3 t2w);
		static TMatrixD getPropagationJacobianGlobalToGlobal(float ds, TVector3 t1w);
		static TVector3 getPositionfromArcLength(TVector3 pos, TVector3 pVec, float beamQ, double s);
		static TVector3 getMomentumfromArcLength(TVector3 momentum, float charge, float arcLength);
		static TVector3 getMomentumfromArcLengthLocal(TVector3 pVec, TVector3 pos, float beamQ, float s, int  planeID);
        static bool findIntersectionWithCertainID(	float x0, float y0, float z0, float px, float py, float pz, float beamQ, int nextPlaneID, float outputPosition[],
TVector3& outputMomentum, float& arcLength, int& newNextPlaneID);

	
	private:
		EUTelNav();
};

}
#endif
