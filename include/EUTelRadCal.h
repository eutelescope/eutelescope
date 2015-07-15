#ifndef EUTELRADCAL_H
#define	EUTELRADCAL_H

#include "EUTelUtility.h"
#include <iostream>
#include "EUTelTrack.h"
#include "EUTelState.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
#include "EUTelGBLFitter.h"
#include "TMath.h"

#include <map>

namespace eutelescope {
  //! This class performs all radiation length calculations. 
  /*! This class will use EUTelTrack objects directly. 
   *  However the analysis code simply uses points in global space to propgate through.
   *  Radiation length is stored in terms of blocks. All sensors included in fit are assumed a thin scatterer.
   *  The radiation length in front of each sensor is assumed thick. Integration over the propagated trajectory for thick scatterers is needed.
   *  After this integration the weighted position and variance of the arclength is determined. This is used to model the thick scatterer by two scattering points.
   *
   */ 


	class  EUTelRadCal {
        public:
            //! Struct object
            /*! This stores the basic information for each block of radiation length.
             *  Each block is either a sensor or the material in front of it.
             *
             */

            struct Block{
                bool isSen;
                double totalRad;
                double weigMean;
                double weigVar;
                std::vector<double> startPos; 
                std::vector<double> endPos;
            };
            //! Set block homogeneous information.
            /*! This function will take the start and end points just outside the detector system.
             *  It will then propagate through the full detector system and get the total radiation length for each block. 
             *  Will only calculate the radiation length if in the system. If you propagate from outside then this is ignored.
             *  @param[in] start start position of propagation 
             *  @param[in] end end position of propagation
             *  @param[out] blocks This is a vector of blocks.  
             */

            float setHomoBlocks(TVector3& start , TVector3& end ,std::vector<int>& sen,std::map<int ,Block>& blocks);
            //! Set radiation information to blocks 
            /*! 
             *  This will use the first and last hit on a track to parameterise a straight line. 
             *  This straight line is then used to propagate through each detector volume and determine the radiation length
             *  @param[in] track EUTelTrack object 
             */
            /// \todo Can make this propagation through the volumes curved by updating the propagtion in TGeo. Simple but not needed now.

            //! Get weighted moments for thick scatterers.
            /*! This function will take each block which is medium. These all will have negative ID.
             *  The start and end positions of each block are used to propagate in steps allowing the radiation length at each position
             *  to weight the position and variance of the arclength
             *  With these the inhomogenous medium can be modelled with two scattering points. 
             *  
             *  @param[in] blocks Filled vector of blocks 
             *  @param[in] beamEnergy This is needed for the radiation length to be converted to a variance.
             */


            float setInHomoBlocks(std::map<int ,Block>& blocks, double beamEnergy);

            void setRad(EUTelTrack& track);

            //! Moves node and calculates the radiation length from that motion. 
            /*! This will take a node and try to find the next node which it can propagate to. 
             *  If it does this over a sensible range then take the node radiation length and distance travelled and move the node. 
             *  If the distance to the next node was too small then propagate manually and update usign this.
             *  @param [in out ] node TGeoNode
             *  @param [out] rad radiation length from node moving to new location
             *  @param [out] distance travelled 
             */
            void move(TGeoNode *node, double& rad, double& dist);

        private:
            //! Fill the blocks in the correct order with homogeneous information 
            /*! This function will update the map from sensor ID to block. Material infront of each sensor is mapped to -ID to block.   
             *  @param [in] geo EUTelGeometryTelescopeGeoDescription passed to get the node location and global position
             *  @param [in] rad pass radition length to update block with 
             *  @param [in] dist Distance travelled by block
             *  @param [out] blocks updated blocks.
             */

            void fillBlocks(geo::EUTelGeometryTelescopeGeoDescription& geo,double & rad, double& dist,std::map<int ,Block>& blocks);
            //! Calculate the meanWeight for the scatter 
            /*! Function will update block using the start and end positions of the Block object.
             *  The arclength weighted with the variance at each point is added to the point
             *  
             *  @param[in] beamEnergy This is needed for the radiation length to be converted to a variance.
             *  @param[in out] itBl This is an iterator of a map from ID to block
             */
            void getDistWeig(double& beamEnergy, std::map<int, Block>::iterator itBl);
            //! Calculate the varWeight for the scatter 
            /*! The function will not get the variance of the arclength weighted with the variance due to scattering. 
             *  
             *  @param[in] beamEnergy This is needed for the radiation length to be converted to a variance.
             *  @param[in out] itBl This is an iterator of a map from ID to block
             */

            void getVarWeig(double& beamEnergy, std::map<int, Block>::iterator itBl);

        protected:
            ///Planes included in fit.
            std::vector<int> _planes;
            std::map<int,Block> _blocks;


	};

}
#endif	/* EUTELRADCAL */
