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

  // Main class 
  class TDSIntegrationStorage {

    friend class TDSPixelsChargeMap;
  
    public:

    // Constructor
    TDSIntegrationStorage(const unsigned int val_integPixelSegmentsAlongL=0, const unsigned int val_integPixelSegmentsAlongW=0, const unsigned int val_integPixelSegmentsAlongH=0)  
      { 
        // Number of one pixel segments - for integration-results storage (DEFAULT: No storage)
        integPixelSegmentsAlongL = (val_integPixelSegmentsAlongL/2)*2;
        if (integPixelSegmentsAlongL >= 2000)
          {
            std::cout << "Too many pixel segments along length (for integration storage)!" << std::endl;
            exit(1);
          }
        integPixelSegmentsAlongW = (val_integPixelSegmentsAlongW/2)*2;
        if (integPixelSegmentsAlongW >= 2000)
          {
            std::cout << "Too many pixel segments along width (for integration storage)!" << std::endl;
            exit(1);
          }
        integPixelSegmentsAlongH = val_integPixelSegmentsAlongH;
        if (integPixelSegmentsAlongH >= 1000)
          {
            std::cout << "Too many pixel segments along height (for integration storage)!" << std::endl;
            exit(1);
          }

        // If only one is != 0 then other dimensions have just 1 segment each!!!
        // if ( integPixelSegmentsAlongL > 0 || integPixelSegmentsAlongW > 0 || integPixelSegmentsAlongH > 0 ) ...

      };

    // Destructor
    ~TDSIntegrationStorage() 
      { 
      };


    inline void setIntegPixelSegmentsAlongL(const unsigned int val = 0 ) 
      { 
        integPixelSegmentsAlongL = (val/2)*2; 
        if (integPixelSegmentsAlongL >= 2000)
          {
            std::cout << "Too many pixel segments along length (for integration storage)!" << std::endl;
            exit(1);
          }
      };
    inline void setIntegPixelSegmentsAlongW(const unsigned int val = 0 ) 
      { 
        integPixelSegmentsAlongW = (val/2)*2; 
        if (integPixelSegmentsAlongW >= 2000)
          {
            std::cout << "Too many pixel segments along width (for integration storage)!" << std::endl;
            exit(1);
          }
      };
    inline void setIntegPixelSegmentsAlongH(const unsigned int val = 0 ) 
      { 
        integPixelSegmentsAlongH = val; 
        if (integPixelSegmentsAlongH >= 1000)
          {
            std::cout << "Too many pixel segments along height (for integration storage)!" << std::endl;
            exit(1);
          }
      };

    // Is result already stored?
    inline bool isResultStored(const unsigned long long int integSegmentID)
      {
        return integResultsForSegments.find(integSegmentID) != integResultsForSegments.end();
      };

    // Function which adds charge contribution to pixels
    // Caller should determine pixel segment for integration results storage
    inline void rememberResult(const unsigned long long int integSegmentID, double integrationResult)
      {
        integResultsForSegments[ integSegmentID ] = integrationResult;
      };

    // Return stored result
    inline double getResult(const unsigned long long int integSegmentID)
      {
        return integResultsForSegments[ integSegmentID ];
      };


    private:


    // For integration-results storage (division of ONE pixel)
    unsigned int integPixelSegmentsAlongL, integPixelSegmentsAlongW, integPixelSegmentsAlongH;

    // Map for storing the results of integration
    // < segmentL*1000000000000 + segmentW*1000000000 + segmentH*1000000 + pixelL*1000 + pixelW , charge> 
    std::map<unsigned long long int, double> integResultsForSegments;

  };

// Outside class definition
// Hash for storage
  inline unsigned long long int integSegmentID(const unsigned int segmentL, const unsigned int segmentW, const unsigned int segmentH, const unsigned int pixelL, const unsigned int pixelW)
    {
      return segmentL*1000000000000ULL + segmentW*1000000000 + segmentH*1000000 + pixelL*1000 + pixelW;
    }
  

}

#endif
