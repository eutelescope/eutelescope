// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *
 * Description:
 * Method for charge distribution in pixels for Tracker Detailed Simulation
 *
 * Author: Piotr Niezurawski
 *
 * First release: 2008-11-07
 *
 * Licence:
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */


#ifndef TDSPixelsChargeMap_H
#define TDSPixelsChargeMap_H 1

// this code is used only if GSL is available
#ifdef USE_GSL

#include <iostream>
#include <fstream>
#include <map>
#include <iterator>
#include <string>
#include <cmath>

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>

#include <TDSStep.h>
#include <TDSIntegrationStorage.h>
#include <TDSPixel.h>
#include <TDSPrecluster.h>


namespace TDS {

  //! Main map type used for pixel storing
  typedef std::map<unsigned long long int, double> type_PixelsChargeMap;

  //! Const used in pixel index coding
  const unsigned long long int tenTo10 = 10000000000ULL;


  //! Pixel Charge Map for Track Detaile Simulation
  /*! Track Detaile Simulation (TDS) tools are developed for realistic
   *  simulation of charged distribution in pixel
   *  detectors. TDSPixelChargeMap is the core of this package.
   *
   *  Each input (Mokka) hit is divided into a number of smaller steps
   *  (based on energy deposit and path length). For each step corresponding
   *  charge deposit in sensor pixels is calculated by 2D integration
   *  of the expected charge density over the pixel surface. Charge
   *  capture in the silicon (signal attenuation) and charge
   *  reflection from the epitaxial layer boundary are taken into
   *  account.
   *
   *  Additional methods are available for taking into account charge
   *  fluctuations, amplification gain, pedestal and noise. A threshold
   *  cut can also be applied before simple clustering algorithm is
   *  used.
   *
   * Comment:
   * Local coordinate system used in calculations
   *
   *   L - length coordinate
   *   W - width coordinate
   *   H - height coordinate
   *
   *
   *    ^
   *    |---------------
   *    |              |
   *  W |   sensitive  |
   *    |              |
   *  0 +--------------->
   *    0      L
   *
   *    ^      L
   *  0-+--------------->
   *    |  sensitive   |
   *  H |              |   H < 0 for sensitive volume!!!
   *    ----------------
   *
   *  @author Piotr Niezurawski, University of Warsaw <mailto:pniez@fuw.edu.pl>
   *  @version $Id: TDSPixelsChargeMap.h,v 1.2 2008-11-12 14:36:57 bulgheroni Exp $
   *
   */
  class TDSPixelsChargeMap {


  public:


    //! Constructor
    TDSPixelsChargeMap(const double length, const double width, const double height, const double firstPixelCoordL=0, const double firstPixelCoordW=0);


    //! Destructor
    ~TDSPixelsChargeMap();


    //! Set pixels pitch
    /*! Set pixels pitch in L and W and check consistency with defined
     * sensor geometry. It is also checked if number of pixels does not
     * exceed maximum allowed value
     */

    void setPixelLength(const double val);
    void setPixelWidth (const double val);


    //! Set local coordinates of the corner of the first pixel (1,1)

    inline void setFirstPixelCoordL(const double val = 0.) { firstPixelCoordL = val; };
    inline void setFirstPixelCoordW(const double val = 0.) { firstPixelCoordW = val; };


    //! Set parameters of charge distribution parametrization
    /*! Parametrization of the charge distribution is based on the MAPS
     *  test results. It is determined by one parameter which is the
     *  assumed charge attenuation length in silicon. Efficiency for
     *  charge reflection off the epitaxial layer boundary is the second
     *  parameter, however its default value is 1.
     */

    inline void setLambda(const double val = 38.)
      { theParamsOfFunChargeDistribution.lambda = val; };

    inline void setReflectedContribution(const double val = 1.0)
      {
        theParamsOfFunChargeDistribution.reflectedContribution = val;
        theParamsOfFunChargeDistribution.addReflectedContribution = true;
      };


    //! Return number of pixels in the sensor (in Length and Width)

    inline const unsigned long int getNumberPixelsAlongL() { return numberPixelsAlongL; };
    inline const unsigned long int getNumberPixelsAlongW() { return numberPixelsAlongW; };


    //! Return parameters of charge distribution parametrization
    inline double getLambda() { return theParamsOfFunChargeDistribution.lambda; };
    inline double getReflectedContribution()
      { return theParamsOfFunChargeDistribution.reflectedContribution; };


    //! Initialization of integration
    /*! To initialize the integration procedure maximum length and
     * maximum charge deposit per track segment have to be defined, as well
     * as number of pixels to be considered in integration (in L and W)
     * and maximum number of calls in numerical integration.
     */
    void initializeIntegration(const double val_integMaxStepInDistance,
                               const double val_integMaxStepInCharge,
                               const unsigned int val_integMaxNumberPixelsAlongL=7,
                               const unsigned int val_integMaxNumberPixelsAlongW=7,
                               const unsigned int val_gsl_calls=10000);

    //! Define integration storage
    /*! To speed up the digitization results of numerical integration
     *  should be stored in dedicated storage. If many sensors of the
     *  same type are digitized, same storage should be used, so
     *  integration results can be shared.
     */

    void setPointerToIntegrationStorage(TDSIntegrationStorage * val_integrationStorage);


    //! Set maximal range of considered pixels during integration
    /*! Considered are integMaxNumberPixelsAlongL/2 down, the same up,
     *  integMaxNumberPixelsAlongW/2 left, the same right from the pixel
     *  under which there is the current point considered. Actual range
     *  can be smaller if the 'core' pixel is near to the layer border.
     */
    void setIntegMaxNumberPixelsAlongL(const unsigned int val = 7);
    void setIntegMaxNumberPixelsAlongW(const unsigned int val = 7);


    //! Set maximum integration step length
    inline void setIntegMaxStepInDistance(const double val = 5.)
      {
        integMaxStepInDistance = val;
        if (integMaxStepInDistance <= 0.)
          {
            std::cout << "Maximal step in distance for integration should be > 0!" << std::endl;
            exit(1);
          };
      };


    //! Set maximum charge deposit in single integration step
    inline void setIntegMaxStepInCharge(const double val = -100.)
      {
        integMaxStepInCharge = val;
        if (integMaxStepInCharge == 0.)
          {
            std::cout << "Maximal step in charge for integration should be non-zero!" << std::endl;
            exit(1);
          };
      };


    //! Add new charge contribution to pixels
    /*! Diffusion of charge deposited in single Geant/Mokka step is
     * calculated and expected charge collected in single pixels is
     * calculated by numerical integration. This is the main 'engine' of
     * the TDS package
     */
    void update(const TDSStep & step);


    //! Print collected charges to ASCII file (for debuging)

    void print(std::string filename = "pixelsChargeMap.out");


    //! Empty pixel charge map

    inline void clear()
      {
        pixelsChargeMap.clear();
      };


    //! Scale charge deposited in map
    /*! New value of total charge is returned
     */
    double scaleCharge(double scaleFactor);

    //! Apply Poisson fluctuations to the charge deposited in single pixels
    /*! It is assumed that the expected stored charge is in the units of
     * elementary charge. If doCleaning flag is set, pixels with no
     * deposited charge are removed from map.
     */
    void applyPoissonFluctuations(bool doCleaning);


    //! Convert deposited charge to expected output charge
    /*! Deposited charge is scaled according to the specified amplifier
     *   gain. Additional fluctuations are added due to gain variations
     *   and noise (gaussian distribution is used in both
     *   cases). Additional offset of the output signal (pedestal) can be
     *   added.
     */

    void applyGain(double gain,double  gainVariation,double  noise,double offset);

    //! Apply threshold cut to all deposits in the map
    /*! Pixels with deposit smaller than given threshold are removed from
     *  the map
     */
    void  applyThresholdCut(double  adcThreshold);


    //! Get charge deposited in single pixel

    double getPixelCharge(unsigned long int indexAlongL, unsigned long int indexAlongW);


    //! Get total charge collected in the whole sensor (map)
    /*! Returned in the total charge stored in map, in units of
     *  elementary charge or in ADC counts, if  applyADC method was
     *  applied.
     */

    double getTotalCharge();


    //! Number of stored (fired) pixels
    /* Returns number of pixels stored in pixel charge map
     */

    inline unsigned int getPixelNumber()
      {
        return pixelsChargeMap.size();
      }


    //! Get vector of pixels stored in the map
    /*! Vector of TDSPixel  is the main output method of the
     *  TDSPixelsMap tool. For each pixel fired its coordinates (indexes
     * along L and W), local coordinates and deposited charge are stored.
     */
    std::vector<TDSPixel > getVectorOfPixels();


    //! Get rectangular precluster (for clustering cheater)
    /*! Create simple rectangular precluster from the charge deposits
     * stored in map. maxNumberOfPixels is the maximal number of pixels
     * with deposit larger than ThresholdCharge. If this number is
     * exceeded (claster too large or many clusters in the map) no
     * cluster is returned.
     */

    TDSPrecluster getPrecluster(unsigned int rectLength, unsigned int rectWidth, unsigned int maxNumberOfPixels, double ThresholdCharge);


  private:


    const double length, width, height;


    double pixelLength, pixelWidth;


    unsigned long int numberPixelsAlongL, numberPixelsAlongW;


    double firstPixelCoordL, firstPixelCoordW;


    // Map of pixels <pixID, pixCharge> for given part of given layer
    type_PixelsChargeMap pixelsChargeMap;


    unsigned int integMaxNumberPixelsAlongL, integMaxNumberPixelsAlongW;


    // Maximal integration steps -- the smallest is chosen
    double integMaxStepInDistance, integMaxStepInCharge;


    // Pointer to storage of the results of integration
    TDSIntegrationStorage * integrationStorage;
    bool useIntegrationStorage;


    // Integration part variables (GSL - C library)
    const gsl_rng_type *gsl_T;
    gsl_rng *gsl_r;
    // Setup of function to integrate
    gsl_monte_function gsl_funToIntegrate;
    // Memory allocation for GSL
    gsl_monte_miser_state *gsl_s;
    // Number of MISER calls per one integration step
    size_t gsl_calls;


    // Parameters of charge distribution
    struct ParamsOfFunChargeDistribution
    {
      double height; // Height of the layer.
      double H; // Height coordinate of the point (H<0 for sensitive volume)
      double lambda;
      double reflectedContribution; // Weight for reflected-charges contribution
      bool   addReflectedContribution;
    } theParamsOfFunChargeDistribution;


    // Function describing charge distribution in the silicon
    // k[0] = L, k[1] = W, dim - dimension, *p - struct with parameters
    static double funChargeDistribution (double *k, size_t dim, void * p)
      {
        ParamsOfFunChargeDistribution * gp = (ParamsOfFunChargeDistribution *)p;
        // Calculate charge distribution
        // MAPS
        double r2 = k[0]*k[0] + k[1]*k[1] + (gp->H)*(gp->H);
        double r  = sqrt(r2);
        // Function = exp(-r/\lambda) |H| / r^3 / (4\pi)
        double A  = exp( - r / gp->lambda ) * std::abs(gp->H) / (r2*r) / (4*M_PI);
        if (gp->addReflectedContribution)
          {
            // Reflection with the same angle (Remember: gp->H < 0  AND gp->height < 0)
            double rimage2 = k[0]*k[0] + k[1]*k[1] + (2*gp->height - gp->H)*(2*gp->height - gp->H);
            double rimage  = sqrt(rimage2);
            A += gp->reflectedContribution * exp( - rimage / gp->lambda ) * std::abs(2*gp->height - gp->H) / (rimage2*rimage) / (4*M_PI);
          };
        //    cout << " r = " << r << " lambda = " << gp->lambda << " H = " << gp->H << endl;
        // Other detectors...
        // Here implementation + conditions at initialization
        // ...

        //    cout << "A= " << A << endl;
        // Because integration hangs up if function result is 0 everywhere in the integration region:
        return (A <= 0.0) ? 1e-100 : A;
        // return A;
        //
        // Tests only:
        // return 1.;
        // return 0.5;
        // return k[0] > 1.5 ? 0. : 1.;
        // return 0.0;
        //
        // WARNING: A dangerous GSL 'feature' found - integration hangs up if function result is 0 everywhere
        // (probably convergence condition does not foresee this case)!!! With __DBL_MIN__ instead of 0.0
        // the integration does not work either.
      };

  };

} // end of TDS namespace


#endif // USE_GSL
#endif // 
