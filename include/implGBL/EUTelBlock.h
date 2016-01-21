#ifndef EUTELBLOCK_H
#define	EUTELBLOCK_H
namespace eutelescope {
    //! Struct object
    /*! Each block contains the information to model the scattering for the sensor and the material infront of it.
     */

    struct Block{
        ///Radiation length of the included sensor.
        double senRadPer;
        ///Radiation length of the thick scatterer. This will be split between two scattering points.
        double medRadPer;
        ///The arclength weighted with the radiation length of each point
        double weigMean;
        ///Variance of arc length weighted with the radiation length at each point
        double weigVar;
        ///Variance of the included sensor
        double senVar;
        ///The position the thick scatter should be placed relative to the included plane.
        std::vector<double> scatPos;
        ///This is the scattering variance which must be give to each scattering point. Must be determined from the total variance and then split between the scatterers.
        std::vector<double> scatVar;
        ///This is needed to work out the scatVar and scatPos
        std::vector<std::pair<double,double> > thicknessAndRad; 

    };
}
#endif
