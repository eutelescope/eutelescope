/* 
   Description: Pixel for Tracker Detailed Simulation - mainly for output

   Author: Piotr Niezurawski

   Date: 2008-11-02

   Comments:
   Coordinates used:
   L - length coordinate
   W - width coordinate

   ^
   |---------------
   |              |
   W |   sensitive  |
   |              |
   0 ---------------->  Corner of pixel (1,1) has coordinates (0.,0.)
   0      L

*/

#ifndef TDSPIXEL_H
#define TDSPIXEL_H

namespace TDS {


  class TDSPixel {

    friend class TDSPixelsChargeMap;

    public:
    // Constructor
    inline TDSPixel() {};

    // Destructor
    inline ~TDSPixel() {};
  
  
    inline unsigned long int getIndexAlongL()
      {
        return indexAlongL;
      };

    inline unsigned long int getIndexAlongW()
      {
        return indexAlongW;
      };

    inline double getCoordL()
      {
        return coordL;
      };

    inline double getCoordW()
      {
        return coordW;
      };

    inline double getCharge()
      {
        return charge;
      }

    void print()
      {
        std::cout << "indexAlongL=" << indexAlongL << " " << "indexAlongW=" << indexAlongW << std::endl;
        std::cout << "coordL=" << coordL << " " << "coordW=" << coordW << std::endl;
        std::cout << "charge=" << charge << std::endl;
      };

    protected:
  
    // Pixel index (L,W)
    unsigned long int indexAlongL, indexAlongW;
    // Coordinates of the geometrical center
    double coordL, coordW;
    // Charge
    double charge;
  };

}

#endif 
