/* EUTelHotPixel.cc
 * $Id$ 
 */

#include "EUTelHotPixel.h"

namespace eutelescope{

    std::string EUTelHotPixel::getRevision()
    {
        return std::string("$Rev: 20110821 $");
    }
         
    EUTelHotPixel::EUTelHotPixel(int pixelRow, int pixelColumn, float pixelFiringFrequency)                
    {
        setHotPixelRow(pixelRow);
        setHotPixelColumn(pixelColumn);
        setHotPixelFrequency(pixelFiringFrequency);
    }
 
    EUTelHotPixel::EUTelHotPixel(int pixelRow, int pixelColumn)                
    {
        setHotPixelRow(pixelRow);
        setHotPixelColumn(pixelColumn);
        setHotPixelFrequency(1.);
    }
    
    EUTelHotPixel::~EUTelHotPixel() { /* no op*/  }
   
  
    int   EUTelHotPixel::getHotPixelRow() const
    {
        return getIntVal( 0 );
    } 
    
    
    int   EUTelHotPixel::getHotPixelColumn() const
    {
        return getIntVal( 1 );
    } 
 
    float EUTelHotPixel::getHotPixelFrequency() const
    {
        return getFloatVal(0);
    } 


        
    void EUTelHotPixel::setHotPixelRow(int pixelRow)     
    {
        obj()->setIntVal( 0, pixelRow);    
    }
    
    
    void EUTelHotPixel::setHotPixelColumn(int pixelColumn)     
    {
        obj()->setIntVal( 1, pixelColumn);    
    }
    
    void EUTelHotPixel::setHotPixelFrequency(float pixelFrequency)     
    {
        obj()->setFloatVal( 0, pixelFrequency);    
    }
    
        
    void EUTelHotPixel::print(  std::ostream& os ) const
    {
        std::streamsize PrecisionBefore = os.precision();
        std::streamsize ColumnBefore = os.width();
        os.precision(3);  os.setf(std::ios_base::fixed);
        os << " HotPixelRow: " << getHotPixelRow() << " HotPixelColumn: " << getHotPixelColumn();
        os << " firing frequency: " << getHotPixelFrequency() << std::endl;
        os.unsetf( std::ios_base::fixed );
        os.precision( PrecisionBefore );
        os.width( ColumnBefore );
    }
    
    
    std::ostream &operator<<(std::ostream &os, const EUTelHotPixel &p) {
	p.print(os);
	return os;
    }

}
