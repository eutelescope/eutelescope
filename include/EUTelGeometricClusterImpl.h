/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELGEOMETRICCLUSTERIMPL_H
#define EUTELGEOMETRICCLUSTERIMPL_H

#include "EUTelGenericSparseClusterImpl.h"
#include "EUTelGeometricPixel.h"

namespace eutelescope {

//! Implementation of a cluster made of sparsified pixels
class EUTelGeometricClusterImpl : public EUTelGenericSparseClusterImpl<EUTelGeometricPixel>
{
  public:
	//! Default constructor
	EUTelGeometricClusterImpl(IMPL::TrackerDataImpl* data);

	//! Destructor
	virtual ~EUTelGeometricClusterImpl();

	void getClusterGeomInfo(float& xPos, float& yPos, float& xSize, float& ySize) const;

	void getGeometricCenterOfGravity(float& xCoG, float& yCoG) const;

  private:
	DISALLOW_COPY_AND_ASSIGN(EUTelGeometricClusterImpl)
};

}
#endif
