/*
* This source code is part of the Eutelescope package of Marlin.
* You are free to use this source files for your own development as
* long as it stays in a public research context. You are not
* allowed to use it for commercial purpose. You must put this
* header with author names in all development based on this file.
*
*/
#ifndef EUTELTRACKERDATAINTERFACERIMPL_TCC
#define EUTELTRACKERDATAINTERFACERIMPL_TCC

namespace eutelescope {

	//default constructor
	template<class PixelType>
	EUTelTrackerDataInterfacerImpl<PixelType>::EUTelTrackerDataInterfacerImpl(IMPL::TrackerDataImpl* data): _trackerData(data), _nElement(), _type(), _pixelVec()
	{
		auto pixel = std::make_unique<PixelType>();
		_nElement = pixel->getNoOfElements();
		_type = pixel->getSparsePixelType();
		_pixelVec.clear();
		fillPixelVec();
	}

	//the amount of pixels is the same as the ones stored in the local copy, the _pixelVec
	template<class PixelType>
	size_t EUTelTrackerDataInterfacerImpl<PixelType>::size() const
	{
		return _pixelVec.size();
	}

	//the amount of pixels is the same as the ones stored in the local copy, the _pixelVec
	template<class PixelType>
	void EUTelTrackerDataInterfacerImpl<PixelType>::addSparsePixel(EUTelBaseSparsePixel* pixel)
	{
		this->addSparsePixel(dynamic_cast<PixelType*>(pixel));
	}



	
} //namespace
#endif
