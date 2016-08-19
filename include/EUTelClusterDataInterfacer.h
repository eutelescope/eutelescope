/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTelClusterDataInterfacer_H
#define EUTelClusterDataInterfacer_H

#include "EUTelTrackerDataInterfacerImpl.h"
#include "EUTELESCOPE.h"

namespace eutelescope {

template<class PixelType> 
class EUTelClusterDataInterfacer
{
  public:
    EUTelClusterDataInterfacer(IMPL::TrackerDataImpl* data): _rawDataInterfacer(data) {
    } 
  protected:
    //! The interfacer to the raw data
    EUTelTrackerDataInterfacerImpl<PixelType> _rawDataInterfacer;

  public:
    //! Get one of the sparse pixel
    /*! This method is used to get one of the sparse pixel contained
     *  into the TrackerData. Not mutant version.
     */ 
    template <typename ...Params> 
    void getSparsePixelAt(Params&&... params) const {
	_rawDataInterfacer.getSparsePixelAt(std::forward<Params>(params)...);
    }

    //! Add a sparse pixel
    /*! This method is used to add to the current TrackerDataImpl a
     *  new sparse pixel with all the pieces of information.
     */
    template <typename ...Params>
    void addSparsePixel(Params&&... params) {
	_rawDataInterfacer.addSparsePixel(std::forward<Params>(params)...);
    }

    template <typename ...Params>
    void push_back(Params&&... params) {
	_rawDataInterfacer.push_back(std::forward<Params>(params)...);
    }

    template <typename ...Params>
    void emplace_back(Params&&... params) {
	_rawDataInterfacer.emplace_back(std::forward<Params>(params)...);
    }

    template <typename ...Params>
    auto at(Params&&... params) const -> decltype(_rawDataInterfacer.at(std::forward<Params>(params)...)){
	return _rawDataInterfacer.at(std::forward<Params>(params)...);
    }

    template <typename ...Params>
    auto operator[](Params&&... params) const -> decltype(_rawDataInterfacer.operator[](std::forward<Params>(params)...)){
	return _rawDataInterfacer.operator[](std::forward<Params>(params)...);
    }

    auto begin() const -> decltype(_rawDataInterfacer.begin()) {
	return _rawDataInterfacer.begin();
    }

    auto end() const -> decltype(_rawDataInterfacer.end()) {
	return _rawDataInterfacer.end();
    }

    std::vector<PixelType> const & getPixels() const {
	return _rawDataInterfacer.getPixels();
    }
    auto size() const -> decltype(_rawDataInterfacer.size()) {
		return _rawDataInterfacer.size();
    }

    auto empty() const -> decltype(_rawDataInterfacer.empty()) {
		return _rawDataInterfacer.empty();
    }

   };

}
#endif
