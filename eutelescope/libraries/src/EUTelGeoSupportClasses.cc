#include "EUTelGeoSupportClasses.h"
#include <Eigen/LU>


namespace eutelescope {
namespace geo{


void EUTelLayer::reseatPosition(Eigen::Vector3d const & posVec) {
	_posVector = posVec;
	//To reseat we need to keep the plane (i.e. its position) at the same location,
	//i.e. the new offset is the position minus the new layer position
	for(auto& plane: _activeVec) {
		auto position = plane->getPosition();
		plane->setOffset(position-_posVector);
	}
	for(auto& plane: _passiveVec) {
		auto position = plane->getPosition();
		plane->setOffset(position-_posVector);
	}
}



void EUTelLayer::reseatRotation(Eigen::Matrix3d const & rotMat) {
	
	_rotMatrix = rotMat;
	_angleVector = Utility::getRotationAnglesFromMatrix( rotMat );

	Eigen::Matrix3d invNewRotMat = rotMat.inverse();
	
	for(auto& plane: _activeVec) {
		auto oldGlobalRotMat = plane->getGlobalRotationMatrix();
		plane->setDeltaRotation( oldGlobalRotMat*invNewRotMat );
	}
	for(auto& plane: _passiveVec) {
		auto oldGlobalRotMat = plane->getGlobalRotationMatrix();
		plane->setDeltaRotation( oldGlobalRotMat*invNewRotMat );
	}
}

}}//namespaces
