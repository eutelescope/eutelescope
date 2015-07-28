#ifndef FEI4SINGLE_H
#define	FEI4SINGLE_H

  /** @class FEI4Single
	* This class is the implementation of  @class EUTelGenericPixGeoDescr
	* for a FEI4 layout with edge pixels which are 400 microns long, the
	* other properties are: 80 x 336 pixels, 250 x 50 microns**2 size
	* with exception of the edge pixels (X=0,79) which are longer.
    */

//EUTELESCOPE
#include "EUTelGenericPixGeoDescr.h"

//ROOT
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"

namespace eutelescope {
    namespace geo {
        
        class FEI4Single : public EUTelGenericPixGeoDescr {
            
        public:
            FEI4Single();
            ~FEI4Single();
            
            void createRootDescr(char const *);
            std::string getPixName(int, int);
            std::pair<int, int> getPixIndex(char const *);
            
        protected:
            TGeoMaterial* matSi;
            TGeoMedium* Si;
            TGeoVolume* plane;
            
        };
        
        class FEI4Single250x50 : public EUTelGenericPixGeoDescr {
            
        public:
            FEI4Single250x50();
            ~FEI4Single250x50();
            
            void createRootDescr(char const *);
            std::string getPixName(int, int);
            std::pair<int, int> getPixIndex(char const *);
            
        protected:
            TGeoMaterial* matSi;
            TGeoMedium* Si;
            TGeoVolume* plane;
            
        };
        
        class FEI4Single125x100 : public EUTelGenericPixGeoDescr {
            
        public:
            FEI4Single125x100();
            ~FEI4Single125x100();
            
            void createRootDescr(char const *);
            std::string getPixName(int, int);
            std::pair<int, int> getPixIndex(char const *);
            
        protected:
            TGeoMaterial* matSi;
            TGeoMedium* Si;
            TGeoVolume* plane;
            
        };
        
        class FEI4Single167x125 : public EUTelGenericPixGeoDescr {
            
        public:
            FEI4Single167x125();
            ~FEI4Single167x125();
            
            void createRootDescr(char const *);
            std::string getPixName(int, int);
            std::pair<int, int> getPixIndex(char const *);
            
        protected:
            TGeoMaterial* matSi;
            TGeoMedium* Si;
            TGeoVolume* plane;
            
        };
        
        class FEI4Single167x125_super : public EUTelGenericPixGeoDescr {
            
        public:
            FEI4Single167x125_super();
            ~FEI4Single167x125_super();
            
            void createRootDescr(char const *);
            std::string getPixName(int, int);
            std::pair<int, int> getPixIndex(char const *);
            
        protected:
            TGeoMaterial* matSi;
            TGeoMedium* Si;
            TGeoVolume* plane;
            
        };
        
        class FEI4Single2000x25 : public EUTelGenericPixGeoDescr {
            
        public:
            FEI4Single2000x25();
            ~FEI4Single2000x25();
            
            void createRootDescr(char const *);
            std::string getPixName(int, int);
            std::pair<int, int> getPixIndex(char const *);
            
        protected:
            TGeoMaterial* matSi;
            TGeoMedium* Si;
            TGeoVolume* plane;
            
        };
        
        class FEI4Single500x25 : public EUTelGenericPixGeoDescr {
            
        public:
            FEI4Single500x25();
            ~FEI4Single500x25();
            
            void createRootDescr(char const *);
            std::string getPixName(int, int);
            std::pair<int, int> getPixIndex(char const *);
            
        protected:
            TGeoMaterial* matSi;
            TGeoMedium* Si;
            TGeoVolume* plane;
            
        };
        
        class FEI4Single500x25_alt : public EUTelGenericPixGeoDescr {
            
        public:
            FEI4Single500x25_alt();
            ~FEI4Single500x25_alt();
            
            void createRootDescr(char const *);
            std::string getPixName(int, int);
            std::pair<int, int> getPixIndex(char const *);
            
        protected:
            TGeoMaterial* matSi;
            TGeoMedium* Si;
            TGeoVolume* plane;
            
        };

extern "C"
{
EUTelGenericPixGeoDescr* maker();
}

} //namespace geo
} //namespace eutelescope

#endif	//FEI4SINGLE_H
