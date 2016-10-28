#ifndef EUTELGENERICPIXGEODESCR_H
#define	EUTELGENERICPIXGEODESCR_H

  /** @class EUTelGenericPixGeoDescr
	* This class is the base class for any description of the pixel geometry
	* used in the pixel geometry framework. It contains three essential virtual
	* methods which have to be implemented by the user information. Additionally
	* it stores information on the pixel count as well as the dimensions of the
	* sensitive area (in mm) and the radiation length of the sensor.
    */

//STL
#include <string>
#include <utility>

//ROOT
#include "TGeoManager.h"

namespace eutelescope {
namespace geo {

class EUTelGenericPixGeoDescr {
	
	public:	

	  /** The only constructor which is public or protected
        * @param are the dimensions of the sensor (size) as well as the minimum and maximum pixel count (min&max) as well as the radiation length
		* Since every sensor has this information, this constructor has to be called whenever creating an instance of a derived class
		*/
		EUTelGenericPixGeoDescr(double sizeX, double sizeY, double sizeZ, int minX, int maxX, int minY, int maxY, double radLen);

	  /** Default deconstructor declared to be virtual */
		virtual ~EUTelGenericPixGeoDescr() {};

	  /** Takes three references to floats and stores the sensitive size of X,Y,Z in them */
		void getSensitiveSize(float& x, float& y, float& z)
		{
			x = _sizeSensitiveAreaX;
			y = _sizeSensitiveAreaY;
			z = _sizeSensitiveAreaZ;
		}

	  /** Signature overloaded version of @see getSensitiveSize() omitting Z information */
		void getSensitiveSize(float& x, float& y)
		{
			x = _sizeSensitiveAreaX;
			y = _sizeSensitiveAreaY;
		}
	
	  /** Takes references to four ints and returns the max and min pixel index in X and Y */
		void getPixelIndexRange(int& minX, int& maxX, int& minY, int& maxY)
		{
			minX= _minIndexX;
			maxX = _maxIndexX;
			minY = _minIndexY;
			maxY = _maxIndexY;
		}

	  /** Takes the char* as a path name for the plane
		* and creates the nodes for the pixel representation
		* in it */
		virtual void createRootDescr(char const *) = 0;


	  /** Signature overloaded version to also take
		* a std::string as an argument */
		void createRootDescr(std::string plane)
		{
			this->createRootDescr( plane.c_str() );	
		};

	  /** Returns the path of pixel @param pixel index as it is
		* represented in the TGeo description	*/
		virtual std::string getPixName(int, int) = 0;


	  /** From a given path (char*) the pixel index is returned */
		virtual std::pair<int, int> getPixIndex(char const *) = 0;

	  /** Signature overloaded version to also take a std::string */
		std::pair<int, int> getPixIndex(std::string path)
		{
			return this->getPixIndex( path.c_str() );
		};

	protected:
		TGeoManager* _tGeoManager;

		double _sizeSensitiveAreaX, _sizeSensitiveAreaY, _sizeSensitiveAreaZ;
		int _minIndexX, _minIndexY;
		int _maxIndexX, _maxIndexY;
		double _radLength;

	private:
	  /** Empty constructor is private, no need to ever call it */
		EUTelGenericPixGeoDescr();
};

} //namespace geo
} //namespace eutelescope

#endif	//EUTELGENERICPIXGEODESCR_H
