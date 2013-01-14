// Version: $Id$
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
 0 +--------------->  Corner of pixel (0,0) has coordinates (0.,0.)
   0      L

*/

#ifndef TDSPIXEL_H
#define TDSPIXEL_H

namespace TDS {


//! Pixel for Tracker Detailed Simulation
/*! 
   Pixel for Tracker Detailed Simulation - mainly for output

   <pre>
   Coordinates used:
   L - length coordinate
   W - width coordinate

   ^
   |---------------
   |              |
 W |   sensitive  |
   |              |
 0 +--------------->  By default corner of the pixel (0,0) has coordinates (0.,0.)
   0      L
   </pre>

   @author Piotr Niezurawski

   Date: 2008-12-10

*/
  class TDSPixel {

    friend class TDSPixelsChargeMap;

    public:

    //! Constructor
    inline TDSPixel() {};


    //! Constructor
    inline TDSPixel(unsigned long int v_indexAlongL, unsigned long int v_indexAlongW, double v_coordL, double v_coordW, double v_charge) 
      :  indexAlongL(v_indexAlongL), indexAlongW(v_indexAlongW), coordL(v_coordL), coordW(v_coordW), charge(v_charge)
      {};


    //! Destructor
    inline ~TDSPixel() {};
  
  
    //! Returns pixel's index along L
    inline unsigned long int getIndexAlongL()
      {
        return indexAlongL;
      };


    //! Returns pixel's index along W
    inline unsigned long int getIndexAlongW()
      {
        return indexAlongW;
      };


    //! Returns pixel's coordinate L (middle of the pixel)
    inline double getCoordL()
      {
        return coordL;
      };


    //! Returns pixel's coordinate W (middle of the pixel)
    inline double getCoordW()
      {
        return coordW;
      };


    //! Returns charge/deposit in the pixel
    inline double getCharge()
      {
        return charge;
      }


    //! Print to stdout info about pixel
    void print()
    {
      std::cout << "indexAlongL=" << indexAlongL << " " << "indexAlongW=" << indexAlongW << std::endl;
      std::cout << "coordL=" << coordL << " " << "coordW=" << coordW << std::endl;
      std::cout << "charge=" << charge << std::endl;
    }
    
    //! Comparison function - used for sorting
    static bool greaterCharge( TDSPixel a, TDSPixel b)
    {
      return a.charge > b.charge;
    }


    protected:
  
    //! Pixel index along L
    unsigned long int indexAlongL;

    //! Pixel index along W
    unsigned long int indexAlongW;

    //! Coordinate along L of the geometrical center
    double coordL;

    //! Coordinate along W of the geometrical center
    double coordW;

    //! Charge
    double charge;
  };

}

#endif 
