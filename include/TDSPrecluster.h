/* 
   Description: Precluster for Tracker Detailed Simulation.

   Author: Piotr Niezurawski

   Date: 2008-10-30
*/

#ifndef TDSPRECLUSTER_H
#define TDSPRECLUSTER_H

#include <vector>

namespace TDS {
  class TDSPrecluster {

    friend class TDSPixelsChargeMap;

    public:
    // Constructor
    inline TDSPrecluster(bool val_empty=true): empty(val_empty) {};

    // Destructor
    inline ~TDSPrecluster() {};
  
    inline bool isEmpty()
      {
        return empty;
      };

    inline unsigned long int getPixelL()
      {
        return pixelL;
      };

    inline unsigned long int getPixelW()
      {
        return pixelW;
      };

    inline double getCoordL()
      {
        return coordL;
      };

    inline double getCoordW()
      {
        return coordW;
      };

    inline unsigned long int getRectLmin()
      {
        return rectLmin;
      };

    inline unsigned long int getRectLmax()
      {
        return rectLmax;
      };

    inline unsigned long int getRectWmin()
      {
        return rectWmin;
      };

    inline unsigned long int getRectWmax()
      {
        return rectWmax;
      };

    inline std::vector<double > getPixelsCharges()
      {
        return pixelsCharges;
      }

    void print()
      {
        std::cout << "pixelL=" << pixelL << " " << "pixelW=" << pixelW << std::endl;
        std::cout << "coordL=" << coordL << " " << "coordW=" << coordW << std::endl;
        std::cout << "pixelsCharges.size()=" << pixelsCharges.size() << std::endl;
      };

    protected:
  
    // Flag
    bool empty;
    // Core pixel (greatest |charge|) - integer coordinates. This is the center of the rectangle (rectLength/2 pixels on the left, rectLength/2 on the right etc.; sometimes - at the layer border some pixels can be absent).
    unsigned long int pixelL, pixelW;
    // Coordinates (at the beginning just center of the core pixel)
    double coordL, coordW;
    // Sides of the rectangle (as you wish)
    unsigned int rectLength, rectWidth;
    // Actual (take into account borders)
    unsigned long int rectLmin, rectLmax, rectWmin, rectWmax;
    // Vector of pixels' charges. Integer coords: (L,W). First in the vector is (1,1) then (1,2), (1,3), (1,rectWidth), (2,1) ...
    std::vector<double > pixelsCharges;
  };

}

#endif
