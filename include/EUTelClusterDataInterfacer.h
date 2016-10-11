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

#include "EUTelTrackerDataInterfacer.h"
#include "EUTelTrackerDataInterfacerImpl.h"
#include "EUTELESCOPE.h"

namespace eutelescope {

class EUTelClusterDataInterfacerBase
{
  public:
	EUTelClusterDataInterfacerBase() = default;
	virtual ~EUTelClusterDataInterfacerBase() = default;

  protected:
	EUTelTrackerDataInterfacer* _rawDataInterfacerRef;
	void setEUTelTrackerDataInterfacerPtr(EUTelTrackerDataInterfacer* ptr) {
		_rawDataInterfacerRef = ptr;
	}

  public:
	void push_back(EUTelBaseSparsePixel const & pixel){
		_rawDataInterfacerRef->push_back(pixel);
	}
	virtual auto at(size_t i) const -> decltype(_rawDataInterfacerRef->at(i)) { //throws std::out_of_range
		return _rawDataInterfacerRef->at(i);
	}

	virtual auto operator[](size_t i) const -> decltype(_rawDataInterfacerRef->operator[](i)){
		return _rawDataInterfacerRef->operator[](i);
	}
	virtual auto size() const -> decltype(_rawDataInterfacerRef->size()) {
		return _rawDataInterfacerRef->size();
	}

	virtual auto empty() const -> decltype(_rawDataInterfacerRef->empty()) {
		return _rawDataInterfacerRef->empty();
	}
};

template<class PixelType> 
class EUTelClusterDataInterfacer: public EUTelClusterDataInterfacerBase 
{
  public:
	EUTelClusterDataInterfacer(IMPL::TrackerDataImpl* data): _rawDataInterfacer(data) {
		setEUTelTrackerDataInterfacerPtr(&_rawDataInterfacer);
	}
	virtual ~EUTelClusterDataInterfacer() = default; 
  protected:
	//! The interfacer to the raw data
	EUTelTrackerDataInterfacerImpl<PixelType> _rawDataInterfacer;

  public:
	template <typename ...Params>
	void push_back(Params&&... params) {
		_rawDataInterfacer.push_back(std::forward<Params>(params)...);
	}

	template <typename ...Params>
	void emplace_back(Params&&... params) {
		_rawDataInterfacer.emplace_back(std::forward<Params>(params)...);
	}

	auto at(size_t i) const -> decltype(_rawDataInterfacer.at(i)) override {
		return _rawDataInterfacer.at(i);
	}

	auto operator[](size_t i) const -> decltype(_rawDataInterfacer.operator[](i)) override {
		return _rawDataInterfacer.operator[](i);
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
	auto size() const -> decltype(_rawDataInterfacer.size()) override {
		return _rawDataInterfacer.size();
	}

	auto empty() const -> decltype(_rawDataInterfacer.empty()) override {
		return _rawDataInterfacer.empty();
	}
};
}
#endif
