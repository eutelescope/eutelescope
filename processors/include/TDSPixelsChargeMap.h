// Version: $Id$
/*
 *
 * Description:
 * Determination of charge distribution in pixels for tracker
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
#include <algorithm>

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>

#include <TDSStep.h>
#include <TDSIntegrationStorage.h>
#include <TDSPixel.h>
#include <TDSPrecluster.h>

//! Namespace
/*!
    Namespace of Tracker Detailed Simulation
*/
namespace TDS {


  //! Type of the pixel ID (identification number for a pixel)
  typedef unsigned long long int type_PixelID;

  //! Main map type used for pixel storing
  typedef std::map<type_PixelID, double> type_PixelsChargeMap;

  //! Const used in pixel index coding
  const unsigned long long int tenTo10 = 10000000000ULL;


  //! Pixels Charge Map for Tracker Detailed Simulation
  /*!
   * <pre>
   *  Track Detailed Simulation (TDS) tools are developed for realistic
   *  simulation of charged distribution in pixel detectors.
   *  TDSPixelChargeMap is the core of this package.
   *
   *  Each input hit (Geant step) is divided into a number of smaller steps
   *  (based on energy deposit and path length). For each step corresponding
   *  charge deposit in sensor pixels is calculated by 2D integration
   *  of the expected charge density over the pixel surface. Charge
   *  capture in the silicon (signal attenuation) and charge
   *  reflection from the epitaxial layer boundary are taken into
   *  account.
   *
   *  Warning: In Mokka simulation the options
   *             /Mokka/init/detailedHitsStoring [detector_name]
   *             /Mokka/init/lcioDetailedTRKHitMode [detector_name]Collection
   *           should be present. For example:
   *             /Mokka/init/detailedHitsStoring VXD
   *             /Mokka/init/lcioDetailedTRKHitMode VXDCollection
   *
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
   *  0 +---------------> By default corner of the pixel (0,0) has coordinates (0.,0.)
   *    0      L
   *
   *    ^      L
   *  0-+--------------->
   *    |  sensitive   |
   *  H |              |   H < 0 for sensitive volume!!!
   *    ----------------
   * </pre>
   *
   *  @author Piotr Niezurawski, University of Warsaw <mailto:pniez@fuw.edu.pl>
   *  @author Contributors: Aleksander Filip Zarnecki, Lukasz Maczewski, Antonio Bulgheroni
   *
   */


  class TDSPixelsChargeMap {


  public:


    //! Constructor.
    /*!
      Constructor takes dimensions (obligatory) and coordinates of the corner (optional).
      Convention: height < 0, as the H axis direction is layer -> readout -> outside, with readout plane at H = 0.
    */
    TDSPixelsChargeMap(const double length, const double width, const double height, const double firstPixelCornerCoordL=0, const double firstPixelCornerCoordW=0);


    //! Destructor
    ~TDSPixelsChargeMap();


    //! Set pixels pitch along Length coordinate (L)
    /*! Set pixels pitch along L and check consistency with defined
     * sensor geometry. It is also checked if number of pixels does not
     * exceed maximum allowed value
     */

    void setPixelLength(const double val);


    //! Set pixels pitch along Width coordinate (W)
    /*! Set pixels pitch along W and check consistency with defined
     * sensor geometry. It is also checked if number of pixels does not
     * exceed maximum allowed value
     */

    void setPixelWidth (const double val);


    //! Set local L coordinate of the corner of the first pixel (0,0)

    inline void setFirstPixelCornerCoordL(const double val = 0.) { firstPixelCornerCoordL = val; };


    //! Set local W coordinate of the corner of the first pixel (0,0)

    inline void setFirstPixelCornerCoordW(const double val = 0.) { firstPixelCornerCoordW = val; };


    //! Choose detector type
    /*! Currently parametrization only for the MAPS-detector type is implemented
     */
    inline void setDetectorType(const std::string val = "MAPS")
    { theParamsOfFunChargeDistribution.detectorType = val; }


    //! Set lambda parameter of charge distribution parametrization - charge attenuation length
    /*! Parametrization of the charge distribution is based on the MAPS
     *  test results. It is determined by one parameter which is the
     *  assumed charge attenuation length in silicon. Efficiency for
     *  charge reflection off the epitaxial layer boundary is the second
     *  parameter, however its default value is 1.
     */

    inline void setLambda(const double val = 38.)
      { theParamsOfFunChargeDistribution.lambda = val; }


    //! Set factor for reflected-charge contribution
    /*! Parametrization of the charge distribution is based on the MAPS
     *  test results. It is determined by one parameter which is the
     *  assumed charge attenuation length in silicon (lambda parameter). Efficiency for
     *  charge reflection off the epitaxial layer boundary is the second
     *  parameter, however its default value is 1.
     */

    inline void setReflectedContribution(const double val = 1.0)
      {
        if ( val != 0.0 )
          {
            theParamsOfFunChargeDistribution.reflectedContribution = val;
            theParamsOfFunChargeDistribution.addReflectedContribution = true;
          }
        else
          {
            theParamsOfFunChargeDistribution.addReflectedContribution = false;
          }
      }


    //! Return number of pixels in the sensor along Length

    inline unsigned long int getNumberPixelsAlongL() { return numberPixelsAlongL; };


    //! Return number of pixels in the sensor along Width

    inline unsigned long int getNumberPixelsAlongW() { return numberPixelsAlongW; };


    //! Return parameter of charge distribution parametrization

    inline double getLambda() { return theParamsOfFunChargeDistribution.lambda; };


    //! Return parameter of charge distribution parametrization

    inline double getReflectedContribution()
      { return theParamsOfFunChargeDistribution.reflectedContribution; };


    //! Initialization of integration (obligatory)
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
    /*! To speed up the digitization, results of numerical integration
     *  should be stored in dedicated storage. If many sensors of the
     *  same type are digitized, same storage should be used, so
     *  integration results can be shared.
     */

    void setPointerToIntegrationStorage(TDSIntegrationStorage * val_integrationStorage);


    //! Set maximal range along L of considered pixels during integration
    /*! Considered are integMaxNumberPixelsAlongL/2 left, the same right,
     *  integMaxNumberPixelsAlongW/2 down, the same up from the pixel
     *  under which there is the current point considered. Actual range
     *  can be smaller if the 'core' pixel is near to the layer border.
     */

    inline void setIntegMaxNumberPixelsAlongL(const unsigned int val = 7);


    //! Set maximal range along W of considered pixels during integration
    /*! Considered are integMaxNumberPixelsAlongL/2 left, the same right,
     *  integMaxNumberPixelsAlongW/2 down, the same up from the pixel
     *  under which there is the current point considered. Actual range
     *  can be smaller if the 'core' pixel is near to the layer border.
     */

    inline void setIntegMaxNumberPixelsAlongW(const unsigned int val = 7);


    //! Set maximum integration step (distance)
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
    inline void setIntegMaxStepInCharge(const double val = 40.)
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


    //! Print collected charges to ASCII file

    void print(std::string filename = "pixelsChargeMap.out");


    //! Empty the pixels' charge map

    inline void clear()
      {
        pixelsChargeMap.clear();
      };


    //! Scale charge deposited in the map
    /*! New value of total charge is returned
     */

    double scaleCharge(double scaleFactor);


    //! Apply Poisson fluctuations to the charge deposited in single pixels
    /*! It is assumed that the expected stored charge is in the units of
     * elementary charge.
     <br>
     * If doCleaning flag is set, pixels with no
     * deposited charge are removed from map.
     */

    void applyPoissonFluctuations(bool doCleaning);


    //! Convert deposited charge to expected output charge
    /*! Deposited charge is scaled according to the specified amplifier
     *   gain. Additional fluctuations are added due to gain variations
     *   and noise (gaussian distribution is used in both
     *   cases). Additional offset of the output signal (pedestal) can be
     *   added. <br>
     Procedure: <br> <pre>

     varGain  is from Gauss(gain, gainVariation)

     varNoise is from Gauss(offset, noise)

     new_charge = varGain * charge + varNoise </pre>
     */

    void applyGain(double gain, double gainVariation, double noise = 0.0, double offset = 0.0);


    //! Apply threshold cut to all deposits in the map
    /*! Pixels with deposit smaller than given threshold are removed from
     *  the map
     */

    void  applyThresholdCut(double threshold);


    //! Get total charge collected in the whole sensor (map)
    /*! Returned in the total charge stored in map, in units of
     *  elementary charge or in ADC counts, if  applyGain method was
     *  applied.
     */

    double getTotalCharge();


    //! Number of stored (fired) pixels
    /*! Returns number of pixels stored in pixels charge map
     */

    inline unsigned int getPixelsNumber()
      {
        return pixelsChargeMap.size();
      }


    //! Is any pixel stored in the map?

    inline bool isPixelStored()
    {
      return pixelsChargeMap.size() != 0;
    }


    //! Is pixel stored in the map?

    inline bool isPixelStored(type_PixelID pixID)
    {
      return pixelsChargeMap.find(pixID) != pixelsChargeMap.end();
    }


    //! Is pixel stored in the map?

    inline bool isPixelStored(unsigned long int indexAlongL, unsigned long int indexAlongW)
    {
      return pixelsChargeMap.find( getPixelID(indexAlongL, indexAlongW) ) != pixelsChargeMap.end();
    }


    //! Translate L and W pixel's indexes to the internal pixel ID

    inline type_PixelID getPixelID(unsigned long int indexAlongL, unsigned long int indexAlongW)
    {
      return tenTo10*indexAlongL + indexAlongW;
    }


    //! Get charge deposited in single pixel
    /*! indexAlongL is a number between 0 and numberPixelsAlongL-1;
        indexAlongW is a number between 0 and numberPixelsAlongW-1
     */

    double getPixelCharge(unsigned long int indexAlongL, unsigned long int indexAlongW);


    //! Get charge deposited in single pixel

    double getPixelCharge(type_PixelID pixID);


    //! Get pixel index along L
    /*! pixelIndexAlongL is a number between 0 and numberPixelsAlongL-1
     */

    unsigned long int getPixelIndexAlongL(type_PixelID pixID);


    //! Get pixel index along W
    /*! pixelIndexAlongW is a number between 0 and numberPixelsAlongW-1
     */

    unsigned long int getPixelIndexAlongW(type_PixelID pixID);


    //! Get pixel coordinate along L (in the middle of the pixel)
    /*! indexAlongL is a number between 0 and numberPixelsAlongL-1
     */

    double getPixelCoordL(unsigned long int indexAlongL);


    //! Get pixel coordinate along L (in the middle of the pixel)

    double getPixelCoordL(type_PixelID pixID);


    //! Get pixel coordinate along W (in the middle of the pixel)
    /*! indexAlongW is a number between 0 and numberPixelsAlongW-1
     */

    double getPixelCoordW(unsigned long int indexAlongW);


    //! Get pixel coordinate along W (in the middle of the pixel)

    double getPixelCoordW(type_PixelID pixID);


    //! Get ID of the pixel with maximal deposit
    /*! Returns ID (internal identification code) of the pixel with the greatest deposit (i.e. greatest absolute value of charge: |10| < |-100|)
     *  or the first pixel of many with the same, biggest deposit.
     */

    type_PixelID getPixelID_maxDeposit();


    //! Comparison function - used for maximum/minimum finding
    static bool smallerDeposit( std::pair<type_PixelID, double>  a, std::pair<type_PixelID, double> b)
    {
      return std::abs(a.second) < std::abs(b.second);
    }


    //! Get ID of the pixel with maximal charge
    /*! Returns ID (internal identification code) of the pixel with the greatest charge (here 10 > -100)
     *  or the first pixel of many with the same, maximal charge.
     */

    type_PixelID getPixelID_maxCharge();


    //! Comparison function - used for maximum/minimum finding
    static bool smallerCharge( std::pair<type_PixelID, double>  a, std::pair<type_PixelID, double> b)
    {
      return a.second < b.second;
    }


    //! Get ID of the pixel with minimal charge
    /*! Returns ID (internal identification code) of the pixel with the smallest charge (here 10 > -100)
     *  or the first pixel of many with the same, minimal charge.
     */

    type_PixelID getPixelID_minCharge();


    //! Erase a pixel
    /*! indexAlongL is a number between 0 and numberPixelsAlongL-1;
        indexAlongW is a number between 0 and numberPixelsAlongW-1
     */

    void erasePixel(unsigned long int indexAlongL, unsigned long int indexAlongW);


    //! Erase a pixel

    void erasePixel(type_PixelID pixID);


    //! Get vector of pixels IDs
    /*! Vector contains current IDs of pixels.
     */

    std::vector<type_PixelID> getVectorOfPixelsIDs();


    //! Get vector of pixels stored in the map
    /*! Vector of TDSPixel  is the main output method of the
     *  TDSPixelsChargeMap tool. For each pixel fired its coordinates (indexes
     * along L and W), local coordinates and deposited charge are stored.
     * Pixels are sorted in charge in descending order (e.g.: +10., +5., -1., -20.).
     * Use for example reverse(v.begin(), v.end()) to have ascending ordering.
     */

    std::vector<TDSPixel> getVectorOfPixels();



    //! Get rectangular precluster
    /*! Create simple rectangular precluster from the charge deposits
     * stored in map.
     * Precluster is created around the pixel chosen by the user.
     * Pixels in the internal vector are sorted in charge in descending order (e.g.: +10., +5., -1., -20.).
     * Use for example reverse(v.begin(), v.end()) to have ascending ordering.
     * Pixels belonging to the cluster are removed from the map if flag 'removePixels' is set to 'true'.
     */

    TDSPrecluster getPrecluster(unsigned long int seedIndexAlongL, unsigned long int seedIndexAlongW, unsigned int rectLength, unsigned int rectWidth, bool removePixels = false);


    //! Get rectangular precluster
    /*! Create simple rectangular precluster from the charge deposits
     * stored in map.
     * Precluster is created around the pixel chosen by the user.
     * Pixels in the internal vector are sorted in charge in descending order (e.g.: +10., +5., -1., -20.).
     * Use for example reverse(v.begin(), v.end()) to have ascending ordering.
     * Pixels belonging to the cluster are removed from the map if flag 'removePixels' is set to 'true'.
     */

    TDSPrecluster getPrecluster(type_PixelID seedPixelID, unsigned int rectLength, unsigned int rectWidth, bool removePixels = false);


    //! Get vector of rectangular preclusters
    /*! Create vector of simple rectangular preclusters from the charge deposits
     * stored in map. Very simple algorithm is implemented at the moment.
     * Precluster is created around the pixel with greatest deposit; pixels belonging to this precluster are then removed from the map; the next hottest pixel is found etc.
     * Preclusters are sorted in charge in descending order (e.g.: +10., +5., -1., -20.).
     * Use for example reverse(v.begin(), v.end()) to have ascending ordering.
     */

    std::vector<TDSPrecluster> getVectorOfPreclusters(unsigned int rectLength, unsigned int rectWidth);


  private:


    const double length, width, height;


    double pixelLength, pixelWidth;

    bool isPixelLengthSet, isPixelWidthSet;


    unsigned long int numberPixelsAlongL, numberPixelsAlongW;


    double firstPixelCornerCoordL, firstPixelCornerCoordW;


    // Map of pixels <pixID, pixCharge> for given part of given layer
    type_PixelsChargeMap pixelsChargeMap;

    bool isIntegrationInitialized;

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


    // Do not insert this in DOXYGEN-generated documentation
    /// @cond
    // Internal stuct - to provide parameters of charge distribution to function through GSL machinery
    struct ParamsOfFunChargeDistribution
    {
      double height; // Height of the layer. Used for reflection contribution calculation.
      double H; // Height coordinate of the point (H<0 for sensitive volume)
      double lambda;
      double reflectedContribution; // Weight for reflected-charges contribution
      bool   addReflectedContribution;
      std::string detectorType;
    } theParamsOfFunChargeDistribution;
    /// @endcond


    // Function describing charge distribution in the silicon
    // k[0] = L, k[1] = W, dim - dimension (not used), *p - struct with parameters
    static double funChargeDistribution (double *k, size_t /* dim */ , void * p)
    {
      ParamsOfFunChargeDistribution * gp = static_cast< ParamsOfFunChargeDistribution* >(p);

      // Variable for result storing
      double A;
      int debug = 0 ;
      // Calculate charge distribution

      if ( gp->detectorType == "MAPS" )
        {
          // MAPS
          double r2 = k[0]*k[0] + k[1]*k[1] + (gp->H)*(gp->H);
          double r  = std::sqrt(r2);
          if(debug)std::cout << "MAPS r = " << r << " lambda = " << gp->lambda << " H = " << gp->H <<std::endl;
          // Function = exp(-r/\lambda) |H| / r^3 / (4\pi)
          A  = std::exp( - r / gp->lambda ) * std::abs(gp->H) / (r2*r) / (4*M_PI);
          if(debug)std::cout << " function A: " << A << "    (-r/la) " <<         ( - r / gp->lambda ) << " " <<          (gp->H) << " " << r2*r << " " << 4*M_PI << std::endl;
          if(debug)std::cout << " function A: " << A << " exp(-r/la) " << std::exp( - r / gp->lambda ) << " " <<  std::abs(gp->H) << " " << r2*r << " " << 4*M_PI << std::endl;
          if (gp->addReflectedContribution)
            {
              // Reflection with the same angle (Remember: gp->H < 0  AND gp->height < 0)
              double rimage2 = k[0]*k[0] + k[1]*k[1] + (2*gp->height - gp->H)*(2*gp->height - gp->H);
              double rimage  = std::sqrt(rimage2);
              A += gp->reflectedContribution * std::exp( - rimage / gp->lambda ) * std::abs(2*gp->height - gp->H) / (rimage2*rimage) / (4*M_PI);
            };
        }
      //
      // Other detectors...
      else if ( gp->detectorType == "FEI4" )
        {
//          std::cout << "Detector Type : " << gp->detectorType << "!" <<std::endl;

          // a copy of MAPS at the moment 
          double r2 = k[0]*k[0] + k[1]*k[1] + (gp->H)*(gp->H);
          double r  = std::sqrt(r2);
          if(debug)std::cout << "MAPS r = " << r << " lambda = " << gp->lambda << " H = " << gp->H <<std::endl;
          // Function = exp(-r/\lambda) |H| / r^3 / (4\pi)
          A  = std::exp( - r / gp->lambda ) * std::abs(gp->H) / (r2*r) / (4*M_PI);
          if(debug)std::cout << " function A: " << A << "    (-r/la) " <<         ( - r / gp->lambda ) << " " <<          (gp->H) << " " << r2*r << " " << 4*M_PI << std::endl;
          if(debug)std::cout << " function A: " << A << " exp(-r/la) " << std::exp( - r / gp->lambda ) << " " <<  std::abs(gp->H) << " " << r2*r << " " << 4*M_PI << std::endl;
          if (gp->addReflectedContribution)
            {
              // Reflection with the same angle (Remember: gp->H < 0  AND gp->height < 0)
              double rimage2 = k[0]*k[0] + k[1]*k[1] + (2*gp->height - gp->H)*(2*gp->height - gp->H);
              double rimage  = std::sqrt(rimage2);
              A += gp->reflectedContribution * std::exp( - rimage / gp->lambda ) * std::abs(2*gp->height - gp->H) / (rimage2*rimage) / (4*M_PI);
            };
           
        } 
      else
        {
          std::cout << "Error: Detector type " << gp->detectorType << " not found!" <<std::endl;
          exit(1);
        }
      //   std::cout << "A= " << A << std::endl;
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

#else
#warning WARNING: GSL not available! TDS will not compile.
#endif // USE_GSL
#endif //
