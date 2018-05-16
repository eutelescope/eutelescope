/*
 * Created by Mykyta Haranko
 *  (2018 DESY)
 *
 *  email:mykyta.haranko@desy.de
 */

#ifndef CBCHITRECOVERY_H
#define CBCHITRECOVERY_H 1

// marlin includes ".h"
#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// system includes <>
#include <string>
#include <map>

namespace eutelescope
{

    //! CMS CBC Hit Recovery processor for Marlin (recovers real track hits from the given position on the dut)

    class CBCHitRecovery : public marlin::Processor
    {

        public:

            virtual Processor * newProcessor ( )
            {
                return new CBCHitRecovery;
            }

            //! Default constructor
            CBCHitRecovery ();

            virtual void init ();

            virtual void processRunHeader (LCRunHeader * run);

            virtual void processEvent (LCEvent * evt);

            virtual void check (LCEvent * evt);

            void bookHistos();

            void fillHistos();

            void vector_set_length(double *& vec, double length);
            double vector_get_length(const double *vec);
            double cos_alpha(const double *vec1,const double *vec2);
            double* hit_pos(const double track[],const double normal[],const double virtual_hit[], double distance);

            virtual void end();

            std::string _InputTrackCollectionName;

            std::string _InputFitHitsCollectionName;

            std::string _cbcInputCollectionName;

            std::string _cbcDataOutputCollectionName;

            int _cbcVirtualDUTId;

            std::vector<int> _cbcRealDUTsVec;

            double *_VirtualDutNormal_zplus;
            double *_VirtualDutNormal_zminus;
            double *_VirtualDutPos;

            gear::SiPlanesParameters * _siPlanesParameters;
            gear::SiPlanesLayerLayout * _siPlanesLayerLayout;

            // map of the vectors
            std::map < int, double* > _dutPosMap;
            std::map < int, double* > _dutNormalMap;

        protected:

            std::map < std::string, AIDA::IBaseHistogram * > _aidaHistoMap;

    };

    //! A global instance of the processor
    CBCHitRecovery gCBCHitRecovery;

}

#endif
