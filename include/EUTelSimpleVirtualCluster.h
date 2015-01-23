//TODO: DOCUMENTATION
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 */

#ifndef EUTELSIMPLEVIRTUALCLUSTER_H
#define EUTELSIMPLEVIRTUALCLUSTER_H

// personal includes ".h"
#include "EUTELESCOPE.h"

// lcio includes <.h>
#include <IMPL/TrackerDataImpl.h>

// system includes <>
#include <iostream>
#include <vector>

namespace eutelescope {

//! Virtual class to describe cluster from ZS data in EUTelescope
/* Currently only implemented in @class EUTelGenericSparseClusterImpl 
 * it should be used as a base class for any cluster classes.
 */
class EUTelSimpleVirtualCluster {

public:
    //! Default constructor
    EUTelSimpleVirtualCluster(IMPL::TrackerDataImpl* data) : _trackerData(data) { } 

    //! Default destructor
    virtual ~EUTelSimpleVirtualCluster() {;}

    //! Get the cluster dimensions
    /*! For each cluster type is always possible to define the
     *  external sizes. 
     *
     *  @param xSize The size along x
     *  @param ySize The size along y
     */ 
    virtual void getClusterSize(int& xSize, int& ySize) const = 0;

    virtual void getCenterOfGravity(float& xCoG, float& yCoG) const = 0;

    //virtual unsigned int size() const = 0;

    //! Return the total charge
    /*!
     *  @return The total integrated charge
     */
    virtual float getTotalCharge() const = 0;

    //! Return a pointer to the TrackerDataImpl
    /*! This method is used to expose to the public the
     *  TrackerDataImpl member.
     *
     *  @return The pointer of _trackerData
     */
    virtual IMPL::TrackerDataImpl* trackerData() = 0;

    //! Print
    /*! This method is used to print out the content of the clusters
     * 
     *  @param os The input output stream
     */
    virtual void print(std::ostream& os) const = 0;
      
    //! Overload of operator<<
    /*! This friend function is the overload of the operator << for
     *  the virtual cluster class. It uses the print method that is
     *  virtually defined for all cluster subclasses.
     *
     *  @param os The input output stream as modified by the print
     *  method
     *  @param clu The cluster to be stream out
     *  @return The output stream
     */ 
    friend  std::ostream& operator<< (std::ostream& os , const EUTelSimpleVirtualCluster & clu )  { clu.print(os); return os; }

private:
	DISALLOW_COPY_AND_ASSIGN(EUTelSimpleVirtualCluster)

protected:
    //! The tracker data member
    /*! This is the core of the decorator pattern design. Whenever an
     *  object deriving from this virtual class is created from an
     *  already existing TrackerData object, this is assigned to
     *  _trackerData by the constructor. 
     */ 
    IMPL::TrackerDataImpl* _trackerData;
  };
}
#endif
