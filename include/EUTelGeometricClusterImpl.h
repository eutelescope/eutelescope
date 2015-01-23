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

  private:
	#ifndef DISALLOW_COPY_AND_ASSIGN
	//Following #define stops the accidental creation of a copy or assignment operator by causing a link error. 
	//Copy and Assignment operators not allowed because they are unnecessary and the cause of many bugs
	#define DISALLOW_COPY_AND_ASSIGN(EUTelGeometricClusterImpl) \
	EUTelGeometricClusterImpl(const EUTelGeometricClusterImpl&); \
	void operator=(const EUTelGeometricClusterImpl&);
	//Private Functions
	DISALLOW_COPY_AND_ASSIGN(EUTelGeometricClusterImpl)//See #define just above
	#endif
};

}
#endif
