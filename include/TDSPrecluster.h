// Version: $Id$
/* 
   Description: Precluster for Tracker Detailed Simulation.

   Author: Piotr Niezurawski

   Date: 2008-10-30
*/

#ifndef TDSPRECLUSTER_H
#define TDSPRECLUSTER_H

#include <vector>
#include <iostream>
#include <map>
#include <iterator>
#include <cmath>

#include <TDSPixel.h>

namespace TDS {


//! Precluster for Tracker Detailed Simulation
/*! 
   Simple rectangular precluster from the charge deposits
   stored in map.
   A few useful methods for charge distribution analysis are present.

   @author Piotr Niezurawski

   Date: 2009-01-13
*/
  class TDSPrecluster {

    friend class TDSPixelsChargeMap;
    
  public:


    //! Constructor
    inline TDSPrecluster(bool val_empty=true): empty(val_empty) {};


    //! Destructor
    inline ~TDSPrecluster() {};
  

    //! Is empty?
    inline bool isEmpty()
      {
        return empty;
      };


    //! Returns seed-pixel's index along L
    inline unsigned long int getSeedIndexAlongL()
      {
        return pixelL;
      };


    //! Returns seed-pixel's index along W
    inline unsigned long int getSeedIndexAlongW()
      {
        return pixelW;
      };


    //! Returns seed-pixel's coordinate L (middle of the seed pixel)
    inline double getSeedCoordL()
      {
        return coordL;
      };


    //! Returns seed-pixel's coordinate W (middle of the seed pixel)
    inline double getSeedCoordW()
      {
        return coordW;
      };


    //! Returns index along L of the leftmost pixels
    inline unsigned long int getRectLmin()
      {
        return rectLmin;
      };


    //! Returns index along L of the rightmost pixels
    inline unsigned long int getRectLmax()
      {
        return rectLmax;
      };


    //! Returns index along W of the bottom pixels
    inline unsigned long int getRectWmin()
      {
        return rectWmin;
      };


    //! Returns index along W of the top pixels
    inline unsigned long int getRectWmax()
      {
        return rectWmax;
      };


    //! Returns L coordinate of the charge center
    inline double getCoordL_chargeCenter()
    {
      return coordL_chargeCenter;
    }


    //! Returns W coordinate of the charge center
    inline double getCoordW_chargeCenter()
    {
      return coordW_chargeCenter;
    }


    //! Returns charge of the precluster
    inline double getCharge()
    {
      return charge;
    }


    //! Vector of pixels belonging to the precluster
    inline std::vector<TDSPixel> getVectorOfPixels()
      {
        return vectorOfPixels;
      }


    //! Vector of pixels' charges sorted in charge in descending order
    std::vector<double> getVecCharges_DescendingInCharge();


    //! Vector of pixels' charges sorted in |charge| in descending order
    std::vector<double> getVecCharges_DescendingInAbsCharge();


    //! Vector of pixels' charges sorted in charge/distance_from_seed in descending order
    std::vector<double> getVecCharges_DescendingInChargeByDistance();


    //! Vector of pixels' charges sorted in |charge|/distance_from_seed in descending order
    std::vector<double> getVecCharges_DescendingInAbsChargeByDistance();


    //! Print to stdout some info about precluster
    void print();


    //! Comparison function - used for sorting
    static bool greaterCharge( TDSPrecluster a, TDSPrecluster b)
    {
      return a.charge > b.charge;
    }


    protected:
  
    //! Flag
    bool empty;


    //! Index of the seed pixel along L. 
    /*! This is the center of the rectangle (rectLength/2 pixels on the left, rectLength/2 on the right etc.; 
        sometimes - at the layer border some pixels can be absent). 
     */
    unsigned long int pixelL;


    //! Index of the seed pixel along W. 
    /*! This is the center of the rectangle (rectLength/2 pixels on the left, rectLength/2 on the right etc.; 
        sometimes - at the layer border some pixels can be absent). 
     */
    unsigned long int pixelW;


    //! Coordinate L of the center of the seed pixel
    double coordL;


    //! Coordinate W of the center of the seed pixel
    double coordW;


    //! Side L of the rectangle (expected)
    unsigned int rectLength;


    //! Side W of the rectangle (expected)
    unsigned int rectWidth;


    //! Actual leftmost-pixels' index of pixels' rectangle (borders are taken into account)
    unsigned long int rectLmin;


    //! Actual rightmost-pixels' index of pixels' rectangle (borders are taken into account)
    unsigned long int rectLmax;


    //! Actual bottom-pixels' index of pixels' rectangle (borders are taken into account)
    unsigned long int rectWmin;


    //! Actual top-pixels' index of pixels' rectangle (borders are taken into account)
    unsigned long int rectWmax;


    //! L coordinate of the charge center 
    /*! Mean of the charge distribution over pixels.
        In calculation it is assumed that charge in the pixel is in the geometrical center of that pixel.
    */
    double coordL_chargeCenter;


    //! W coordinate of the charge center
    /*! Mean of the charge distribution over pixels.
        In calculation it is assumed that charge in the pixel is in the geometrical center of that pixel.
    */
    double coordW_chargeCenter;

    
    //! Total charge of the precluster
    double charge;


    //! Vector of pixels belonging to the precluster
    /*! Provide pixels sorted in charge in descending order (e.g.: +10., +5., -1., -20.).
     */

    std::vector<TDSPixel> vectorOfPixels;


    //! Vector of pixels' charges. Integer coords: (L,W). First in the vector is (0,0) then (0,1), (0,2), (1,rectWidth-1), (1,0) ...
    // std::vector<double> pixelsCharges;
    
    
  };

}

#endif
