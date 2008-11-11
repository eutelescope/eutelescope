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
0 ---------------->
  0      L

  ^      L        
0-|--------------->
  |  sensitive   |
H |              |   H < 0 for sensitive volume!!! 
  ----------------

*/

#include <CLHEP/Vector/ThreeVector.h>

class TDSStep {

  friend class TDSPixelsChargeMap;

 public:
  // Constructor
  inline TDSStep(const double midL, const double midW, const double midH, const double dirL, const double dirW, const double dirH, const double geomLength, const double charge) : 
  midL(midL), midW(midW), midH(midH), dirL(dirL), dirW(dirW), dirH(dirH), geomLength(geomLength), charge(charge) {};
  // Destructor
  inline ~TDSStep() {};

  // Gets
  inline const double getMidL() { return midL; };
  inline const double getMidW() { return midW; };
  inline const double getMidH() { return midH; };

  inline const CLHEP::Hep3Vector getMidPoint()  { CLHEP::Hep3Vector vec(midL, midW, midH);  return vec; } ;
 
  // Sets
  inline void setMidL(const double val) { midL = val; };
  inline void setMidW(const double val) { midW = val; };
  inline void setMidH(const double val) { midH = val; };

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
  };
 
  inline void setDir(const CLHEP::Hep3Vector val) 
  { 
    dirL = val.getX(); 
    dirW = val.getY(); 
    dirH = val.getZ();
  };
 

  protected:
  // Middle point of step
  double midL;
  double midW;
  double midH;
  // Direction of step (normalized to length 1)
  double dirL;
  double dirW;
  double dirH;
  // 
  // double trueLength; // Takes scattering into account, unused at the moment
  double geomLength; // Necessary, unfortunately at the moment Mokka provides only trueLenght and temporarily this is used here
  // 
  double charge;

};
