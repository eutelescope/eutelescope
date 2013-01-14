// Version: $Id$
/* 
   Description: Store of integration results - may be used for many layers with the same pixel sizes, lambda, height etc.
   Divide each pixel into sectors -- integration results are stored and reused.
   Only even numbers should be considered for L and W (therefore (val/2)*2).

   Author: Piotr Niezurawski

   Date: 2008-10-29

*/

#ifndef TDSIntegrationStorage_H
#define TDSIntegrationStorage_H 1

#include <map>

namespace TDS {

//! Storage of integration results for Tracker Detailed Simulation
/*! 
   A storage for integration results - may be used for many layers with the same pixel sizes, lambda, height etc.
   <br>
   Divide each pixel into sectors (segments) -- integration results are stored and reused.
   <br>
   Only even numbers should be considered for L and W (therefore (val/2)*2).

   @author Piotr Niezurawski

   Date: 2008-10-29

*/
  class TDSIntegrationStorage {

    friend class TDSPixelsChargeMap;
  
    public:

    //! Constructor
    TDSIntegrationStorage(const unsigned int val_integPixelSegmentsAlongL=0, const unsigned int val_integPixelSegmentsAlongW=0, const unsigned int val_integPixelSegmentsAlongH=0)  
      { 
        // Number of one pixel segments - for integration-results storage (DEFAULT: No storage)
        integPixelSegmentsAlongL = (val_integPixelSegmentsAlongL/2)*2;
        if (integPixelSegmentsAlongL > 2000)
          {
            std::cout << "Too many pixel segments along length (for integration storage)!" << std::endl;
            exit(1);
          }
        integPixelSegmentsAlongW = (val_integPixelSegmentsAlongW/2)*2;
        if (integPixelSegmentsAlongW > 2000)
          {
            std::cout << "Too many pixel segments along width (for integration storage)!" << std::endl;
            exit(1);
          }
        integPixelSegmentsAlongH = val_integPixelSegmentsAlongH;
        if (integPixelSegmentsAlongH > 1000)
          {
            std::cout << "Too many pixel segments along height (for integration storage)!" << std::endl;
            exit(1);
          }

        // If only one is != 0 then other dimensions have just 1 segment each!!!
      };

    //! Destructor
    ~TDSIntegrationStorage() 
      { 
      };

    //! Number of segments for ONE pixel along L (integration results are stored for each segment)
    inline void setIntegPixelSegmentsAlongL(const unsigned int val = 0 ) 
      { 
        integPixelSegmentsAlongL = (val/2)*2; 
        if (integPixelSegmentsAlongL > 2000)
          {
            std::cout << "Too many pixel segments along length (for integration storage)!" << std::endl;
            exit(1);
          }
      };


    //! Number of segments for ONE pixel along W (integration results are stored for each segment)
    inline void setIntegPixelSegmentsAlongW(const unsigned int val = 0 ) 
      { 
        integPixelSegmentsAlongW = (val/2)*2; 
        if (integPixelSegmentsAlongW > 2000)
          {
            std::cout << "Too many pixel segments along width (for integration storage)!" << std::endl;
            exit(1);
          }
      };


    //! Number of segments for ONE pixel along H (integration results are stored for each segment)
    inline void setIntegPixelSegmentsAlongH(const unsigned int val = 0 ) 
      { 
        integPixelSegmentsAlongH = val; 
        if (integPixelSegmentsAlongH > 1000)
          {
            std::cout << "Too many pixel segments along height (for integration storage)!" << std::endl;
            exit(1);
          }
      };


    //! Is the result already stored?
    inline bool isResultStored(const unsigned long long int integSegmentID)
      {
        return integResultsForSegments.find(integSegmentID) != integResultsForSegments.end();
      };

    //! Store a charge deposit from the segment in the pixel
    /*! Caller should determine pixel's segment for integration results storage
     */
    inline void rememberResult(const unsigned long long int integSegmentID, double integrationResult)
      {
        integResultsForSegments[ integSegmentID ] = integrationResult;
      };


    //! Return a stored result
    inline double getResult(const unsigned long long int integSegmentID)
      {
        return integResultsForSegments[ integSegmentID ];
      };


    private:


    //! For integration-results storage - number of segments/divisions of ONE pixel 
    unsigned int integPixelSegmentsAlongL, integPixelSegmentsAlongW, integPixelSegmentsAlongH;

    //! Map for storing the results of integration
    /*! Function integSegmentID() is used as a key 
      map< segmentL*1000000000000 + segmentW*1000000000 + segmentH*1000000 + pixelL*1000 + pixelW , charge> 
     */
    std::map<unsigned long long int, double> integResultsForSegments;

  };


//! Key/hash for storage of integration results. Global TDS function
  inline unsigned long long int integSegmentID(const unsigned int segmentL, const unsigned int segmentW, const unsigned int segmentH, const unsigned int pixelL, const unsigned int pixelW)
    {
      if(0)
      {
      std::cout << " W " <<  pixelW << std::endl;
      std::cout << " LW " <<  pixelL*1000 + pixelW << std::endl;
      std::cout << " HLW " <<  segmentH << " = " << segmentH*1000000ULL << " " << pixelL*1000  << " " << pixelW << std::endl;
      std::cout << " WHLW " <<  segmentW*1000000000ULL + segmentH*1000000ULL + pixelL*1000 + pixelW << std::endl;
      std::cout << " LWHLW " << segmentL*1000000000000ULL + segmentW*1000000000ULL + segmentH*1000000ULL + pixelL*1000 + pixelW << std::endl;
      }

      return segmentL*1000000000000ULL + segmentW*1000000000ULL + segmentH*1000000ULL + pixelL*1000 + pixelW;
    }
  

}

#endif
