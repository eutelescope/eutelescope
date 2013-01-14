// Version: $Id$
#ifndef EUTelHotPixel_h
#define EUTelHotPixel_h 1

//C++
#include <iostream>
#include <string>
#include <utility>

//LCIO
#include "lcio.h"
#include "UTIL/LCFixedObject.h"

#define EUTelHotPixelNINTVals 2       // HotPixel coordinates: Row, Column
#define EUTelHotPixelNFLOATVals 1     // HotPixel firing frequency        
#define EUTelHotPixelNDOUBLEVals 0

namespace eutelescope
{

    class EUTelHotPixel;

    typedef EUTelHotPixel LCHotPixel;

    /** Class that combines information on hot pixels for storage.
     *  The list of pixel Row,Column,Frequency-Firing is stored as   
     *  HotPixel into a single object of type LCGenericObject.
     *  Based on the LCFixedObject template.
     *  This can be used to be stored in a LCCollectionVec.
     * @author I.Rubinskiy (DESY) 
     *  -- based on the Pedestal class for TPC written by M.E.Jansenn and R.Diener                    
     */
    
    class EUTelHotPixel : public UTIL::LCFixedObject<EUTelHotPixelNINTVals,EUTelHotPixelNFLOATVals,EUTelHotPixelNDOUBLEVals>
    {

        public:

            /** Convenient constructor using to integers per pixel and one float to mark its firing frequency
             */
            EUTelHotPixel(int pixelRow, int pixelColumn, float pixelFiringFrequency);                

            /** Convenient constructor using to integers per pixel and one float to mark its firing frequency
             */
            EUTelHotPixel(int pixelRow, int pixelColumn);


            /** 'Copy constructor' needed to interpret LCCollection read from file/database.
             */
            EUTelHotPixel(EVENT::LCObject* obj):UTIL::LCFixedObject<EUTelHotPixelNINTVals,
                                                               EUTelHotPixelNFLOATVals,
                                                               EUTelHotPixelNDOUBLEVals>(obj) { }

            /** Important for memory handling*/
            virtual ~EUTelHotPixel();

                          
            int   getHotPixelRow()    const;
            int   getHotPixelColumn() const;
            float getHotPixelFrequency() const;


            void setHotPixelRow   (int   hotpixelRow)  ;
            void setHotPixelColumn(int   hotpixelColumn);
            void setHotPixelFrequency(float hotpixelFrequency);


            static std::string getRevision();

            void print(  std::ostream& os = std::cout ) const;

            // -------- need to implement abstract methods from LCGenericObject
            const std::string getTypeName() const
            {
                return std::string("SiHotPixel");
            }
            const std::string getDataDescription() const
            {
                return std::string("i:PixelRow,PixelColumn:HotPixelFiringFrequency"    );
            }
            
    };// end of class

    std::ostream &operator<<(std::ostream &os, const EUTelHotPixel &p);    

} //end namespace tpcconddata
#endif
