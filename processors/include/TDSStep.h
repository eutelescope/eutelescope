// Version: $Id$
/* 
   Description: Step for Tracker Detailed Simulation (translate to it other (Mokka/Geant)steps)

   Author: Piotr Niezurawski

   Date: 2008-09-25

   Comments:
   Coordinates used:
   L - length coordinate
   W - width coordinate
   H - height coordinate


   ^
   |---------------
   |              |
 W |   sensitive  |
   |              |
 0 +--------------->
   0      L

   ^      L        
 0-+--------------->
   |  sensitive   |
 H |              |   H < 0 for sensitive volume!!! 
   ----------------

*/

#include<cstdlib>

#ifndef TDSSTEP_H 
#define TDSSTEP_H

// this code is built only if CLHEP is available

#ifdef USE_CLHEP
#include <CLHEP/Vector/ThreeVector.h>

namespace TDS {

//! Step for Tracker Detailed Simulation
/*! 
   Step for Tracker Detailed Simulation. Translate to it other (Mokka/Geant) steps.

   <pre>
   Comments:
   Coordinates used:
   L - length coordinate
   W - width coordinate
   H - height coordinate


   ^
   |---------------
   |              |
 W |   sensitive  |
   |              |
 0 +--------------->
   0      L

   ^      L        
 0-+--------------->
   |  sensitive   |
 H |              |   H < 0 for sensitive volume!!! 
   ----------------
   </pre>

   @author Piotr Niezurawski

   Date: 2008-09-25


*/
  class TDSStep {

    friend class TDSPixelsChargeMap;

    public:

    //! Constructor
    inline TDSStep(const double midL, const double midW, const double midH, const double dirL, const double dirW, const double dirH, const double geomLength, const double charge) : 
      midL(midL), midW(midW), midH(midH), dirL(dirL), dirW(dirW), dirH(dirH), geomLength(geomLength), charge(charge) 
      {
	if (midH > 0.)
	  {
	    std::cout << "Middle point (" << midH << ") of the step should have H coordinate < 0!" << std::endl;
	    exit(1);
	  }
      }


    //! Destructor
    inline ~TDSStep() {}


    //! Get L coordinate of the middle point of the step
    inline double getMidL() { return midL; }


    //! Get W coordinate of the middle point of the step
    inline double getMidW() { return midW; }


    //! Get H coordinate of the middle point of the step
    inline double getMidH() { return midH; }


    //! Get vector to the middle point of the step
    inline const CLHEP::Hep3Vector getMidPoint()  { CLHEP::Hep3Vector vec(midL, midW, midH);  return vec; }
 


    //! Set L coordinate of the middle point of the step
    inline void setMidL(const double val) { midL = val; }


    //! Set W coordinate of the middle point of the step
    inline void setMidW(const double val) { midW = val; }


    //! Set H coordinate of the middle point of the step
    inline void setMidH(const double val) 
    { 
      midH = val; 
      if (midH > 0.)
	{
	  std::cout << "Middle point of the step should have H coordinate < 0!" << std::endl;
	  exit(1);
	}
    }

    
    //! Set vector to the middle point of the step
    inline void setMid(const CLHEP::Hep3Vector val) 
    { 
        midL = val.getX(); 
        midW = val.getY(); 
        midH = val.getZ();
        if (midH > 0.)
          {
            std::cout << "Middle point of the step should have H coordinate < 0!" << std::endl;
            exit(1);
          }
    }
    

    //! Set L coordinate of the direction vector of the step (vector should have length 1)
    inline void setDirL(const double val) { dirL = val; }


    //! Set W coordinate of the direction vector of the step (vector should have length 1)
    inline void setDirW(const double val) { dirW = val; }


    //! Set H coordinate of the direction vector of the step (vector should have length 1)
    inline void setDirH(const double val) { dirH = val; }


    //! Set the direction vector (vector should have length 1) of the step
    inline void setDir(const CLHEP::Hep3Vector val) 
      { 
        dirL = val.getX(); 
        dirW = val.getY(); 
        dirH = val.getZ();
      }


    //! Get direction vector of the step
    inline const CLHEP::Hep3Vector getDirVector()  { CLHEP::Hep3Vector vec(dirL, dirW, dirH);  return vec; }


    //! Normalize the direction vector of the step
    inline void normalizeDir() 
      {
	double length = std::sqrt( dirL*dirL + dirW*dirW + dirH*dirH );
	if ( length > 0. )
	  {
	    dirL /= length; 
	    dirW /= length;  
	    dirH /= length;  
	  }
	else
	  {
	    std::cout << "Warning: Step length <= 0." << std::endl;
	  }
      }


    //! Set geometrical length of the step
    inline void setGeomLength(const double val) { geomLength = val; }
 

    //! Set charge of the step
    inline void setCharge(const double val) { charge = val; }
 

    protected:


    //! L coordinate of the middle point of the step
    double midL;


    //! W coordinate of the middle point of the step
    double midW;


    //! H coordinate of the middle point of the step
    double midH;


    //! L coordinate of the direction vector of the step (vector should have length 1)
    double dirL;


    //! W coordinate of the direction vector of the step (vector should have length 1)
    double dirW;


    //! H coordinate of the direction vector of the step (vector should have length 1)
    double dirH;

    // 
    // double trueLength; // Takes scattering into account, unused at the moment

    //! Geometrical length of a step
    double geomLength; // Necessary. Unfortunately at the moment Mokka provides only trueLenght (pathLength) and temporarily this should be used here
    // 

    //! Charge of the step
    double charge;

  };

}

#else 
#WARNING TDSStep will not be built because CLHEP is not available
#endif // USE_CLHEP
#endif //
