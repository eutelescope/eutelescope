// Version: $Id$
/*!

Description: Map of charge distribution in pixels for Tracker Detailed Simulation

Author: Piotr Niezurawski

Date: 2008-11-07

*/

#include <TDSPixelsChargeMap.h>

#ifdef USE_CLHEP
#include <CLHEP/Random/RandGauss.h>
#include <CLHEP/Random/RandPoisson.h>

#include "marlin/Processor.h"


using namespace TDS;
using namespace std;

// Constructor
TDSPixelsChargeMap::TDSPixelsChargeMap(const double length, const double width, const double height, const double firstPixelCornerCoordL, const double firstPixelCornerCoordW) :
  length(length), width(width), height(height), firstPixelCornerCoordL(firstPixelCornerCoordL), firstPixelCornerCoordW(firstPixelCornerCoordW)
{
  std::cout << " booking width="<< width << " and length=" << length << std::endl;
  if( height > 0. )
    {
      cout << "Error: height must be negative!" << endl;
      exit(1);
    };
 
  clear();

  setDetectorType();
  theParamsOfFunChargeDistribution.height = height;
  setLambda();
  setReflectedContribution();

  // By default no integration storage is used
  useIntegrationStorage = false;

  // Integration should be initialized by user
  isIntegrationInitialized = false;

  // Pixels' dimensions not set
  isPixelLengthSet = isPixelWidthSet = false;
}

// Destructor
TDSPixelsChargeMap::~TDSPixelsChargeMap()
{
  if (isIntegrationInitialized)
    {
      gsl_monte_miser_free (gsl_s); // Free allocated memory (GSL)
    }
}


// Layer of pixels
void TDSPixelsChargeMap::setPixelLength(const double val)
{
  pixelLength = val;
  numberPixelsAlongL = static_cast< unsigned long int >(length/pixelLength);
  if( static_cast< double >(numberPixelsAlongL) < length/pixelLength )
    {
      cout << "Warning: Layer length is not an integer multiple of pixel length." << endl;
    };
  // Check if the system limits are not too low
  if( numberPixelsAlongL > 999999999ULL )
    {
      cout << "Too many pixels to consider (more than 999999999 at least along length)!" << endl;
      exit(1);
    }
  isPixelLengthSet = true;
}


void TDSPixelsChargeMap::setPixelWidth (const double val) {
  pixelWidth = val;
  numberPixelsAlongW = static_cast< unsigned long int >(width/pixelWidth);
  if( static_cast< double >(numberPixelsAlongW) < width/pixelWidth )
    {
      cout << "Warning: Layer width (" << val << ")is not an integer multiple of pixel width ("<< width<<") width/pixelWidth="<<(width/pixelWidth)<<". numberPixelsAlongW="<< static_cast< double >(numberPixelsAlongW) << endl;
    }
  // Check if the system limits are not too low
  if( numberPixelsAlongW > 999999999ULL )
    {
      cout << "Too many pixels to consider (more than 999999999 at least along width)!" << endl;
      exit(1);
    }
  isPixelWidthSet = true;
}


// Initialization of integration
void TDSPixelsChargeMap::initializeIntegration(const double val_integMaxStepInDistance, const double val_integMaxStepInCharge, const unsigned int val_integMaxNumberPixelsAlongL, const unsigned int val_integMaxNumberPixelsAlongW, const unsigned int val_gsl_calls)
{
  // Integration parameters
  integMaxStepInDistance = val_integMaxStepInDistance;
  if (integMaxStepInDistance <= 0.)
    {
      cout << "Maximal step in distance for integration should be > 0!" << endl;
      exit(1);
    };
  integMaxStepInCharge   = val_integMaxStepInCharge;
  if (integMaxStepInCharge == 0.)
    {
      cout << "Maximal step in charge for integration should be non-zero!" << endl;
      exit(1);
    };

  integMaxNumberPixelsAlongL = ((val_integMaxNumberPixelsAlongL-1)/2)*2+1;
  if (integMaxNumberPixelsAlongL > 1999)
    {
      cout << "Too many pixels to integrate over along length!" << endl;
      exit(1);
    }
  integMaxNumberPixelsAlongW = ((val_integMaxNumberPixelsAlongW-1)/2)*2+1;
  if (integMaxNumberPixelsAlongW > 1999)
    {
      cout << "Too many pixels to integrate over along width!" << endl;
      exit(1);
    }
  // GSL random number generator
  gsl_rng_env_setup ();
  gsl_T = gsl_rng_default;
  gsl_r = gsl_rng_alloc (gsl_T);
  // Setup of function to integrate (2-dim integration)
  gsl_funToIntegrate.f = &funChargeDistribution;
  gsl_funToIntegrate.dim = 2;
  gsl_funToIntegrate.params = &theParamsOfFunChargeDistribution;
  // Memory allocation for GSL
  gsl_s = gsl_monte_miser_alloc (2);
  // Number of MISER calls per one integration
  gsl_calls = val_gsl_calls;

  isIntegrationInitialized = true;
}


// Integration storage
void TDSPixelsChargeMap::setPointerToIntegrationStorage(TDSIntegrationStorage * val_integrationStorage)
{
  if (val_integrationStorage != NULL)
    {
      integrationStorage = val_integrationStorage;
      useIntegrationStorage = true;
    }
  else
    {
      cout << "setPointerToIntegrationStorage: Provide non-NULL pointer to integration storage" << endl;
      exit(1);
    }
}

// Maximal range of considered pixels during integration (integMaxNumberPixelsAlongL/2 down, the same up, integMaxNumberPixelsAlongW/2 left, the same right from the pixel under which there is the current point considered). Range can be smaller if the 'core' pixel is near to the layer border.
void TDSPixelsChargeMap::setIntegMaxNumberPixelsAlongL(const unsigned int val)
{
  integMaxNumberPixelsAlongL = ((val-1)/2)*2+1;
  if (integMaxNumberPixelsAlongL > 1999)
    {
      cout << "Too many pixels to integrate over along length!" << endl;
      exit(1);
    }
}

void TDSPixelsChargeMap::setIntegMaxNumberPixelsAlongW(const unsigned int val)
{
  integMaxNumberPixelsAlongW = ((val-1)/2)*2+1;
  if (integMaxNumberPixelsAlongW > 1999)
    {
      cout << "Too many pixels to integrate over along width!" << endl;
      exit(1);
    }
}


// Function which adds charge contribution to pixels
void TDSPixelsChargeMap::update(const TDSStep & step)
{
  int debug = 0;

  if ( ( ! isPixelLengthSet ) || ( ! isPixelWidthSet ) )
    {
      cout << "Error: Pixels' dimensions are not set!" << endl;
      exit (1);
    }

  if ( ! isIntegrationInitialized )
    {
      cout << "Error: Integration is not initialized!" << endl;
      exit(1);
    }

  if (step.geomLength == 0.)    // Return (and take next step)
    {
//      cout << "Warning: Step length = 0." << endl;
//      return;
    }

  if (step.geomLength < 0.)
    {
      cout << "Error: Step length less than 0!" << endl;
      exit(1);
    }

  // Choose the smallest step -- the greatest number of steps
  unsigned int integStepsNumber;
  unsigned int temp1 = static_cast< unsigned int >( step.geomLength / integMaxStepInDistance + 0.5);
  unsigned int temp2 = static_cast< unsigned int >( abs(step.charge / integMaxStepInCharge)  + 0.5);
  unsigned int temp  = max(temp1,temp2);
  integStepsNumber = ( temp > 0 ? temp : 1 );

  double integStep = step.geomLength / integStepsNumber;
  if(debug>1) cout << "integ: length: " << step.geomLength << " stepnumber: " << integStepsNumber << endl;
  // Charge per integration step
  // Charge per integration step
  double integChargePerStep = step.charge / integStepsNumber;
  if(debug>1) cout << "integChargePerStep= " << integChargePerStep << ";  integStep= " << integStep << endl;


  // Initialize position before integration loop (one integration point back)
  double currentPoint[3];
  currentPoint[0] = step.midL - step.dirL*(step.geomLength + integStep)/2.;
  currentPoint[1] = step.midW - step.dirW*(step.geomLength + integStep)/2.;
  currentPoint[2] = step.midH - step.dirH*(step.geomLength + integStep)/2.;

  // Go through points - integration along the step
  for (unsigned int is = 0; is < integStepsNumber ; is++ )
    {

      // new position
      currentPoint[0] += step.dirL*integStep;
      currentPoint[1] += step.dirW*integStep;
      currentPoint[2] += step.dirH*integStep;

      if(debug>1) cout << " currentPoint[0] = " << currentPoint[0] << "  currentPoint[1] = " << currentPoint[1] << " currentPoint[2] = " << currentPoint[2] << endl;

      // Determine integer coordinates of the main (core) pixel (under which the current point is placed)
      unsigned long int iL, iW;
      iL = static_cast< unsigned long int >((currentPoint[0]-firstPixelCornerCoordL)/pixelLength);
      iW = static_cast< unsigned long int >((currentPoint[1]-firstPixelCornerCoordW)/pixelWidth);
      if(debug>1) cout << "iL: " << iL << " " << static_cast< long unsigned int >((currentPoint[0]-firstPixelCornerCoordL)/pixelLength)  << " " ;
      if(debug>1) cout << "iW: " << iW << " " << static_cast< long unsigned int >((currentPoint[1]-firstPixelCornerCoordW)/pixelWidth)  << endl;
      if ( iL >= numberPixelsAlongL  || iW >= numberPixelsAlongW )
        {
          cout << "Error: Core pixel (and step) outside the boundary of Length-Width plane!" << endl;
          //      exit(1);
          break;
        }


      if (currentPoint[2] > 0.)
        {
          cout << "Error: Point outside sensitive volume (Height > 0)!" << endl;
          //      exit(1);
          streamlog_out(ERROR4) << 
                                  " currentPoint[0] " << currentPoint[0] <<
                                  " currentPoint[1] " << currentPoint[1] <<
                                  " currentPoint[2] " << currentPoint[2] <<
                                    endl;
          break;
        }

      // Set H for funChargeDistribution
      theParamsOfFunChargeDistribution.H=currentPoint[2];

      // Pixel segment for integration storage
      unsigned int segmentL=0, segmentW=0, segmentH=0;
      bool segmentL_reduced = false, segmentW_reduced = false;
      if (useIntegrationStorage)
        {
          // Determine pixel segment for integration results storage
          // Unreduced segments
          segmentL = static_cast< unsigned int >( integrationStorage->integPixelSegmentsAlongL * ((currentPoint[0]-firstPixelCornerCoordL-pixelLength*iL ) / pixelLength) );
          segmentW = static_cast< unsigned int >( integrationStorage->integPixelSegmentsAlongW * ((currentPoint[1]-firstPixelCornerCoordW-pixelWidth *iW ) / pixelWidth ) ) ;
          segmentH = static_cast< unsigned int >( integrationStorage->integPixelSegmentsAlongH * (abs(currentPoint[2]) / abs(height) ) );
          if(debug>1) std::cout << "segmentH: " << segmentH  << " point2: " << abs(currentPoint[2]) <<
                    " " << integrationStorage->integPixelSegmentsAlongH << " 1./height:" << 1./height <<  std::endl;
          // Thanks to symmetry we can reduce L and W segments (we have to reduce pixels then, too!)
          if (segmentL >= integrationStorage->integPixelSegmentsAlongL/2)
            {
              segmentL = integrationStorage->integPixelSegmentsAlongL - segmentL - 1;
              segmentL_reduced = true;
            };
          if (segmentW >= integrationStorage->integPixelSegmentsAlongW/2)
            {
              segmentW = integrationStorage->integPixelSegmentsAlongW - segmentW - 1;
              segmentW_reduced = true;
            }
        }

      // Result of integration and its error
      double gsl_res, gsl_err;

      /* Tables for 2-dim limits of integration on the Length-Width plane*/
      double limitsLow[2];
      double limitsUp[2];

      unsigned long int i,j, imin, imax, jmin, jmax; // Indexes of pixels for which contributions will be calculated
      // Limits take into account borders of the pixel plane (rectangle)
      long int temp;
      temp = iL - integMaxNumberPixelsAlongL / 2;
      temp < 0 ? imin = 0 : imin = temp;
      temp = iL + integMaxNumberPixelsAlongL / 2;
      temp >= static_cast< long int >(numberPixelsAlongL) ? imax = numberPixelsAlongL - 1 : imax = temp;
      temp = iW - integMaxNumberPixelsAlongW / 2;
      temp < 0 ? jmin = 0 : jmin = temp;
      temp = iW + integMaxNumberPixelsAlongW / 2;
      temp >= static_cast< long int >(numberPixelsAlongW) ? jmax = numberPixelsAlongW - 1 : jmax = temp;

      // Loops over important pixels
      // (borders of a layer part taken into account - see above)
      if(debug>2) cout << "iL = " << iL << " iW = " << iW << " imin = " << imin << " imax = " << imax << " jmin = " << jmin << " jmax = " << jmax << endl;

      // I use map<unsigned long long int pixId, double pixCharge> to keep charges collected in pixels.
      // Machine limit for unsigned long int:
      // ULLONG_MAX = 18 446 744 073 709 551 615
      //                           ^           ^
      // pixID = 10^10*i + j - key for the pixel (i,j) which is used in map<> container [(i,j) <-> (L,W)]
      type_PixelID pixID;

      for (i = imin ; i <= imax ; i++ )
        {
          if(debug>2) cout << "i = " << i << endl;
          // L limits of integral
          limitsLow[0] = firstPixelCornerCoordL + i*pixelLength - currentPoint[0];
          limitsUp[0]  = limitsLow[0] + pixelLength;

          for (j = jmin ; j <= jmax ; j++ )
            {
              if(debug>2) cout << "j = " << j << endl;
              // W limits of integral
              limitsLow[1] = firstPixelCornerCoordW + j*pixelWidth - currentPoint[1];
              limitsUp[1]  = limitsLow[1] + pixelWidth;

              if(debug>2) std::cout << "limitsLow " << limitsLow[0] << " " << limitsLow[1] << std::endl; 
              if(debug>2) std::cout << "limitsUp  " << limitsUp[0] << " " << limitsUp[1] << std::endl; 
 
              // Should we use integration-results?
              if (useIntegrationStorage)
                {

                  if(debug>2) cout << "Integration-results storage is used!" << endl;

                  // Relative integer coordinates of pixel from main pixel
                  unsigned int pixelL, pixelW;
                  // Thanks to symmetry we can reduce L and W pixels indexes. We have to reduce segments simultaneously!

                  if(debug>2) std::cout << " " << i << " " << iL << " " << integMaxNumberPixelsAlongL / 2 << std::endl;
                  pixelL = static_cast< long int >(i) - static_cast< long int >(iL) + integMaxNumberPixelsAlongL / 2;
                  if ( segmentL_reduced   &&  i != iL )
                    {
                      pixelL = integMaxNumberPixelsAlongL - pixelL - 1;
                    }

                  if(debug>2) std::cout << " " << j << " " << iW << " " << integMaxNumberPixelsAlongW / 2 << std::endl;
                  pixelW = static_cast< long int >(j) - static_cast< long int >(iW) + integMaxNumberPixelsAlongW / 2;
                  if ( segmentW_reduced   &&  j != iW )
                    {
                      pixelW = integMaxNumberPixelsAlongW - pixelW - 1;
                    }

                  if(debug>2) std::cout << "L: " << segmentL << " W: " << segmentW << " H: " << segmentH << " pixelL:" << pixelL << " pixelW:" << pixelW << std::endl;
                  unsigned long long int integSegID = integSegmentID(segmentL, segmentW, segmentH, pixelL, pixelW);

                  if(debug>2) std::cout << " " << pixelL << " " << pixelW << " " << integSegID << std::endl;
                  if ( integrationStorage->isResultStored(integSegID) )
                    {
                      if(debug>2) std::cout << "1here no storage " << std::endl; 
                      gsl_res = integrationStorage->getResult(integSegID);
                    }
                  else
                    {
                      if(debug>2) std::cout << "2here no storage " << std::endl; 
                      // Integrate
                      gsl_monte_miser_integrate (&gsl_funToIntegrate, limitsLow, limitsUp, 2, gsl_calls, gsl_r, gsl_s, &gsl_res, &gsl_err);
                      if(debug>2) std::cout << " integration done    gsl_calls: " << gsl_calls <<
                                "   r: " << gsl_r << "   s: " << gsl_s << "   res: " << gsl_res << "   err:" << gsl_err << std::endl;
                      // Store integration result
                      integrationStorage->rememberResult(integSegID,gsl_res);
                    };
                }
              else
                {
                  if(debug>2) std::cout << " here no storage " << std::endl; 
                  // Integrate (here no storage)
                  gsl_monte_miser_integrate (&gsl_funToIntegrate, limitsLow, limitsUp, 2, gsl_calls, gsl_r, gsl_s, &gsl_res, &gsl_err);
                };

              if(debug>2)cout << "integChargePerStep " << integChargePerStep << " " << gsl_res << std::endl;
               
              // Pixels Charge Map
              // "code" of the pixel - it serves as a key in the map container of pixels (relations: i <-> L, j <-> W)
              pixID = 0UL + tenTo10*i + j ;
              // Add contribution to pixelsChargeMap
              pixelsChargeMap[ pixID ] = pixelsChargeMap[ pixID ] + ( gsl_res * integChargePerStep );
              if(debug>2)cout << "charge collected in the pixel " << pixID << " = " << pixelsChargeMap[ pixID ] << endl;
              if(debug>2)cout << "miser= " << gsl_res << " +- " << gsl_err << endl;
            }
        }
    }
}


void TDSPixelsChargeMap::print(string filename)
{
  ofstream fout(filename.c_str());

  type_PixelsChargeMap::iterator i;
  for( i = pixelsChargeMap.begin(); i != pixelsChargeMap.end(); ++i )
    {
      fout << (i->first)/tenTo10 << "\t" << (i->first)%tenTo10 << "\t" << (i->second) << endl;
    }

}


double TDSPixelsChargeMap::getPixelCharge(unsigned long int indexAlongL, unsigned long int indexAlongW)
{
  type_PixelID pixID;
  pixID = tenTo10*indexAlongL + indexAlongW ;
  if (pixelsChargeMap.find(pixID) == pixelsChargeMap.end() )
    {
      return 0.;
    }
  else
    {
      return pixelsChargeMap[pixID];
    }
}


double TDSPixelsChargeMap::getPixelCharge(type_PixelID pixID)
{
  if (pixelsChargeMap.find(pixID) == pixelsChargeMap.end() )
    {
      cout << "Error: pixID not found in the map!" << endl;
      exit(1);
    }
  else
    {
      return pixelsChargeMap[pixID];
    }
}


unsigned long int TDSPixelsChargeMap::getPixelIndexAlongL(type_PixelID pixID)
{
  if (pixelsChargeMap.find(pixID) == pixelsChargeMap.end() )
    {
      cout << "Error: pixID not found in the map!" << endl;
      exit(1);
    }
  else
    {
      return pixID/tenTo10;
    }
}  


unsigned long int TDSPixelsChargeMap::getPixelIndexAlongW(type_PixelID pixID)
{
  if (pixelsChargeMap.find(pixID) == pixelsChargeMap.end() )
    {
      cout << "Error: pixID not found in the map!" << endl;
      exit(1);
    }
  else
    {
      return pixID%tenTo10;
    }
}  


double TDSPixelsChargeMap::getPixelCoordL(unsigned long int indexAlongL)
{
  if ( indexAlongL >= numberPixelsAlongL )
    {
      cout << "Error: Pixel index outside the allowed range!" << endl;
      exit(1);
    }
  else
    {
      return (static_cast< double >(indexAlongL)+0.5)*pixelLength + firstPixelCornerCoordL;
    }
}  


double TDSPixelsChargeMap::getPixelCoordL(type_PixelID pixID)
{
  if (pixelsChargeMap.find(pixID) == pixelsChargeMap.end() )
    {
      cout << "Error: pixID not found in the map!" << endl;
      exit(1);
    }
  else
    {
      return (static_cast< double >(pixID/tenTo10)+0.5)*pixelLength + firstPixelCornerCoordL;
    }
}  


double TDSPixelsChargeMap::getPixelCoordW(unsigned long int indexAlongW)
{
  if ( indexAlongW >= numberPixelsAlongW )
    {
      cout << "Error: Pixel index outside the allowed range!" << endl;
      exit(1);
    }
  else
    {
      return (static_cast< double >(indexAlongW)+0.5)*pixelWidth + firstPixelCornerCoordW;
    }
}  


double TDSPixelsChargeMap::getPixelCoordW(type_PixelID pixID)
{
  if (pixelsChargeMap.find(pixID) == pixelsChargeMap.end() )
    {
      cout << "Error: pixID not found in the map!" << endl;
      exit(1);
    }
  else
    {
      return (static_cast< double >(pixID%tenTo10)+0.5)*pixelWidth + firstPixelCornerCoordW;
    }
}  


type_PixelID TDSPixelsChargeMap::getPixelID_maxDeposit()
{
  type_PixelsChargeMap::iterator i_maxDep;

  i_maxDep = max_element(pixelsChargeMap.begin(), pixelsChargeMap.end(), TDSPixelsChargeMap::smallerDeposit);

  return i_maxDep->first;
}


type_PixelID TDSPixelsChargeMap::getPixelID_maxCharge()
{
  type_PixelsChargeMap::iterator i_maxCharge;

  i_maxCharge = max_element(pixelsChargeMap.begin(), pixelsChargeMap.end(), TDSPixelsChargeMap::smallerCharge);

  return i_maxCharge->first;
}


type_PixelID TDSPixelsChargeMap::getPixelID_minCharge()
{
  type_PixelsChargeMap::iterator i_minCharge;

  i_minCharge = min_element(pixelsChargeMap.begin(), pixelsChargeMap.end(), TDSPixelsChargeMap::smallerCharge);

  return i_minCharge->first;
}


void TDSPixelsChargeMap::erasePixel(type_PixelID pixID)
{
  pixelsChargeMap.erase(pixID);
}


void TDSPixelsChargeMap::erasePixel(unsigned long int indexAlongL, unsigned long int indexAlongW)
{
  type_PixelID pixID = tenTo10*indexAlongL + indexAlongW;
  erasePixel(pixID);

}


std::vector<type_PixelID> TDSPixelsChargeMap::getVectorOfPixelsIDs()
{
  type_PixelsChargeMap::iterator i;

  vector<type_PixelID> vectorOfPixelsIDs;

  // Speed-up vector filling
  vectorOfPixelsIDs.reserve(pixelsChargeMap.size());

  for( i = pixelsChargeMap.begin(); i != pixelsChargeMap.end(); ++i )
    {
      vectorOfPixelsIDs.push_back(i->first);
    }
  return vectorOfPixelsIDs;
}



// Get total charge collected in the whole map
double TDSPixelsChargeMap::getTotalCharge()
{
  type_PixelsChargeMap::iterator i;

  double totalCharge = 0.;
  int debug = 0;

  int ipixel=0;
  if(debug) streamlog_out ( MESSAGE5 ) << " pixelsChargeMap : " << pixelsChargeMap.size() << endl;

  for( i = pixelsChargeMap.begin(); i != pixelsChargeMap.end(); ++i )
    {
      if(debug) 
           {
             std::cout << "ipixel " << ipixel << " charge" << i->second << std::endl; 
           }
      totalCharge += i->second;
      ipixel++;
    }
  return totalCharge;
}


// Scale charge deposited in map
double TDSPixelsChargeMap::scaleCharge(double scaleFactor)
{
  type_PixelsChargeMap::iterator i;

  double totalCharge = 0.;

  for( i = pixelsChargeMap.begin(); i != pixelsChargeMap.end(); ++i )
    {
      i->second *= scaleFactor;
      totalCharge += i->second;
    }
  return totalCharge;
}

  // Apply Poisson fluctuations to the charge deposited in single pixels

  void TDSPixelsChargeMap::applyPoissonFluctuations(bool doCleaning)
{
  // As pixels will be removed we can not do simple for loop

  type_PixelsChargeMap::iterator i = pixelsChargeMap.begin();

  while( i != pixelsChargeMap.end() )
     {
     double charge =  abs(i->second);
     double varCharge;
     if (charge > 1000.)
       { // assume Gaussian
	 double sigma = std::sqrt(charge);
	 varCharge = double(CLHEP::RandGauss::shoot(charge,sigma));
       }
     else
       { // assume Poisson
	 varCharge = double(CLHEP::RandPoisson::shoot(charge));
       }

     if ( i->second < 0.) varCharge = -varCharge;

      if (varCharge == 0. && doCleaning)
        pixelsChargeMap.erase(i++);
      else
        {
	  i->second = varCharge;
	  ++i;
        }
    }

  return;
}

  // Convert deposited charge to expected charge output

void TDSPixelsChargeMap::applyGain(double gain, double gainVariation, double noise, double offset)
{
  type_PixelsChargeMap::iterator i;

  for( i = pixelsChargeMap.begin(); i != pixelsChargeMap.end(); ++i )
    {
      double charge =  i->second;
      
      double varGain = double(CLHEP::RandGauss::shoot(gain,gainVariation));
      
      double varNoise = double(CLHEP::RandGauss::shoot(offset,noise));
      
      i->second = varGain*charge + varNoise;
    }

  return;
}


// Apply threshold cut to all deposits in the map

void TDSPixelsChargeMap::applyThresholdCut(double threshold)
{
  // As pixels will be removed we can not do simple for loop

  type_PixelsChargeMap::iterator i = pixelsChargeMap.begin();

  while( i != pixelsChargeMap.end() )
    {
    double charge =  i->second;

    if ( (threshold > 0 && charge < threshold) ||
         (threshold < 0 && charge > threshold) )
      {
	pixelsChargeMap.erase(i++);
      }
    else
      {
	++i;
      }
    }

  return;
}


// Get vector of pixels
vector<TDSPixel> TDSPixelsChargeMap::getVectorOfPixels()
{
  type_PixelsChargeMap::iterator i;

  TDSPixel thePixel;
  vector<TDSPixel> vectorOfPixels;

  // Speed-up vector filling
  vectorOfPixels.reserve(pixelsChargeMap.size());

  int ipixel = 0;
  for( i = pixelsChargeMap.begin(); i != pixelsChargeMap.end(); ++i )
    {
      thePixel.indexAlongL = i->first/tenTo10;
      thePixel.indexAlongW = i->first%tenTo10;
      thePixel.coordL = (static_cast< double >(thePixel.indexAlongL)+0.5)*pixelLength + firstPixelCornerCoordL;
      thePixel.coordW = (static_cast< double >(thePixel.indexAlongW)+0.5)*pixelWidth  + firstPixelCornerCoordW;
      thePixel.charge = i->second;
      vectorOfPixels.push_back(thePixel);
      //streamlog_out( MESSAGE5 ) << "ipixel= " << ipixel << " " << thePixel.coordL << " " << thePixel.coordW << endl;
      ipixel++;
    }

  // Sort pixels in charge in descending order.
  // Use for example reverse(v.begin(), v.end()) to have ascending ordering.
  sort(vectorOfPixels.begin(), vectorOfPixels.end(), TDSPixel::greaterCharge);
  
  return vectorOfPixels;
}



// Get rectangular precluster
TDSPrecluster TDSPixelsChargeMap::getPrecluster(unsigned long int seedIndexAlongL, unsigned long int seedIndexAlongW, unsigned int rectLength, unsigned int rectWidth, bool removePixels)
{
  // Checks
  if (rectLength < 1)
    {
      cout << "Error: rectLength < 1" << endl;
      exit(1);
    }
  if (rectWidth < 1)
    {
      cout << "Error: rectWidth < 1" << endl;
      exit(1);
    }
  if ( seedIndexAlongL >= numberPixelsAlongL )
    {
      cout << "Error: Seed-pixel-index outside the allowed range!" << endl;
      exit(1);
    }
  if ( seedIndexAlongW >= numberPixelsAlongW )
    {
      cout << "Error: Seed-pixel-index outside the allowed range!" << endl;
      exit(1);
    }

  TDSPrecluster thePrecluster;
  
  // Fill precluster
  thePrecluster.pixelL = seedIndexAlongL;
  thePrecluster.pixelW = seedIndexAlongW;
  thePrecluster.coordL = getPixelCoordL(seedIndexAlongL);
  thePrecluster.coordW = getPixelCoordW(seedIndexAlongW);
  thePrecluster.rectLength = rectLength;
  thePrecluster.rectWidth  = rectWidth;

  unsigned  long int l, lmin, lmax, w, wmin, wmax;
  long int temp;
  temp = thePrecluster.pixelL-rectLength/2;
  temp < 0 ? lmin = 0 : lmin = temp;
  temp = thePrecluster.pixelL+rectLength/2;
  temp >= static_cast< long int >(numberPixelsAlongL) ? lmax = numberPixelsAlongL - 1 : lmax = temp;
  temp = thePrecluster.pixelW-rectWidth/2;
  temp < 0 ? wmin = 0 : wmin = temp;
  temp = thePrecluster.pixelW+rectWidth/2;
  temp >= static_cast< long int >(numberPixelsAlongW) ? wmax = numberPixelsAlongL - 1 : wmax = temp;

  thePrecluster.rectLmin = lmin;
  thePrecluster.rectLmax = lmax;
  thePrecluster.rectWmin = wmin;
  thePrecluster.rectWmax = wmax;

  // Speed-up vector filling
  thePrecluster.vectorOfPixels.reserve((lmax-lmin+1)*(wmax-wmin+1));

  // For center of charge calculations
  double tempL = 0., tempW = 0.;
  double preclusterCharge = 0., tempCharge;
      
  for (l=lmin; l<=lmax; l++)
    {
      for (w=wmin; w<=wmax; w++)
	{
	  if ( isPixelStored(l,w) )
	    {
	      tempCharge = getPixelCharge(l,w);
	      preclusterCharge += tempCharge;
	      tempL += tempCharge * getPixelCoordL(l);
	      tempW += tempCharge * getPixelCoordW(w);
	      
	      // Fill vector of pixels
	      thePrecluster.vectorOfPixels.push_back( TDSPixel( l, w, getPixelCoordL(l), getPixelCoordW(w), getPixelCharge(l,w) ) );

	      if ( removePixels ) erasePixel(l,w);
	      
	    }
	} // for w
    } // for l

  thePrecluster.charge = preclusterCharge;

  if ( !thePrecluster.vectorOfPixels.empty() )
    {
      thePrecluster.empty = false;
  
      thePrecluster.coordL_chargeCenter = tempL / preclusterCharge;
      thePrecluster.coordW_chargeCenter = tempW / preclusterCharge;
  
      // Sort pixels in descending order.
      // Use for example reverse(v.begin(), v.end()) to have ascending ordering.
      sort(thePrecluster.vectorOfPixels.begin(), thePrecluster.vectorOfPixels.end(), TDSPixel::greaterCharge);
    }

  return thePrecluster;
}


// Get rectangular precluster
TDSPrecluster TDSPixelsChargeMap::getPrecluster(type_PixelID seedPixelID, unsigned int rectLength, unsigned int rectWidth, bool removePixels)
{
  return getPrecluster(getPixelIndexAlongL(seedPixelID), getPixelIndexAlongW(seedPixelID), rectLength, rectWidth, removePixels);
}


//! Get vector of rectangular preclusters
/*! Create vector of simple rectangular preclusters from the charge deposits
 * stored in map. Very simple algorithm is implemented at the moment.
 * Precluster is created around the pixel with greatest deposit; pixels belonging to this precluster are then removed from the map; the next hottest pixel is found etc.
 * Preclusters are sorted in charge in descending order (e.g.: +10., +5., -1., -20.).
 * Use for example reverse(v.begin(), v.end()) to have ascending ordering.
 */

vector<TDSPrecluster> TDSPixelsChargeMap::getVectorOfPreclusters(unsigned int rectLength, unsigned int rectWidth)
{

  type_PixelsChargeMap::iterator i;

  TDSPrecluster thePrecluster;
  vector<TDSPrecluster> theVectorOfPreclusters;

  // Continue if there are still pixels in the map
  while ( getPixelsNumber() > 0 )
    {
      // Find hottest pixel's ID
      type_PixelID hottestPixID = getPixelID_maxDeposit();

      // Create precluster, remove pixels
      thePrecluster = getPrecluster(hottestPixID, rectLength, rectWidth, true);
      
      // Fill vector of preclusters
      theVectorOfPreclusters.push_back(thePrecluster);
    }

  // Sort preclusters in charge in descending order.
  // Use for example reverse(v.begin(), v.end()) to have ascending ordering.
  sort(theVectorOfPreclusters.begin(), theVectorOfPreclusters.end(), TDSPrecluster::greaterCharge);

  return theVectorOfPreclusters;
}


#else
#warning TDSPixelsChargeMap will not be built because CLHEP is not available
#endif
