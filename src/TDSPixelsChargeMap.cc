/* 
Description: Map of charge distribution in pixels for Tracker Detailed Simulation 

Author: Piotr Niezurawski

Date: 2008-11-07
*/

#include "TDSPixelsChargeMap.h"

#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h" 

// Constructor
TDSPixelsChargeMap::TDSPixelsChargeMap(const double length, const double width, const double height, const double firstPixelCoordL, const double firstPixelCoordW) : 
  length(length), width(width), height(height), firstPixelCoordL(firstPixelCoordL), firstPixelCoordW(firstPixelCoordW) 
{ 
  if( height > 0. ) 
    { 
      cout << "Error: height must be negative!" << endl; 
      exit(1); 
    };
    
  // By default no integration storage is used
  useIntegrationStorage = false;

  theParamsOfFunChargeDistribution.height = height;
  theParamsOfFunChargeDistribution.lambda = 30.;
  theParamsOfFunChargeDistribution.reflectedContribution = 0.;
  theParamsOfFunChargeDistribution.addReflectedContribution = false;
}

// Destructor
TDSPixelsChargeMap::~TDSPixelsChargeMap() 
{ 
  gsl_monte_miser_free (gsl_s); // Free allocated memory (GSL)
}


// Layer of pixels
void TDSPixelsChargeMap::setPixelLength(const double val) 
{ 
  pixelLength = val; 
  numberPixelsAlongL = (unsigned long int)(length/pixelLength); 
  if( (double)numberPixelsAlongL < length/pixelLength ) 
    {
      cout << "Warning: Layer length is not an integer multiple of pixel length." << endl; 
    };
  // Check if the system limits are not too low
  if( numberPixelsAlongL > 999999999ULL ) 
    {
      cout << "Too many pixels to consider (more than 999999999 at least along length)!" << endl;
      exit(1);
    }
}

  
void TDSPixelsChargeMap::setPixelWidth (const double val) { 
  pixelWidth = val; 
  numberPixelsAlongW = (unsigned long int)(width/pixelWidth); 
  if( (double)numberPixelsAlongW < width/pixelWidth ) 
    {
      cout << "Warning: Layer width is not an integer multiple of pixel width." << endl; 
    }
  // Check if the system limits are not too low
  if( numberPixelsAlongW > 999999999ULL ) 
    {
      cout << "Too many pixels to consider (more than 999999999 at least along width)!" << endl;
      exit(1);
    }
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
  if (step.geomLength == 0.) 	// Return (and take next step)
    {
      cout << "Warning: Step length = 0." << endl;
      return;
    }
    
  if (step.geomLength < 0.)
    {
      cout << "Error: Step length less than 0!" << endl;
      exit(1);
    }

  // Choose the smallest step -- the greatest number of steps
  unsigned int integStepsNumber;
  unsigned int temp1 = (unsigned int) ( step.geomLength / integMaxStepInDistance + 0.5); 
  unsigned int temp2 = (unsigned int) ( abs(step.charge / integMaxStepInCharge)  + 0.5); 
  unsigned int temp  = max(temp1,temp2);
  integStepsNumber = ( temp > 0 ? temp : 1 );

  double integStep = step.geomLength / integStepsNumber;
  // cout << integStepsNumber << endl;
  // Charge per integration step
  // Charge per integration step
  double integChargePerStep = step.charge / integStepsNumber;
  // cout << "integChargePerStep= " << integChargePerStep << ";  integStep= " << integStep << endl;

  
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

      // cout << " currentPoint[0] = " << currentPoint[0] << "  currentPoint[1] = " << currentPoint[1] << " currentPoint[2] = " << currentPoint[2] << endl;

      // Determine integer coordinates of the main (core) pixel (under which the current point is placed)
      unsigned long int iL, iW; 
      iL = (unsigned long int)((currentPoint[0]-firstPixelCoordL)/pixelLength) + 1; 
      iW = (unsigned long int)((currentPoint[1]-firstPixelCoordW)/pixelWidth) + 1;
      // cout << "iL: " << iL << " " << (long unsigned int)((currentPoint[0]-firstPixelCoordL)/pixelLength)  << endl;
      // cout << "iW: " << iW << " " << (long unsigned int)((currentPoint[1]-firstPixelCoordW)/pixelWidth)  << endl;
      if ( iL < 1  ||  iL > numberPixelsAlongL  ||  iW < 1  ||  iW > numberPixelsAlongW )
	{
	  cout << "Error: Pixel outside the boundary of Length-Width plane!" << endl;
	  //	  exit(1);
          break;
	}
	

      if (currentPoint[2] > 0.)
	{
	  cout << "Error: Point outside sensitive volume (Height > 0)!" << endl;
	  //	  exit(1);
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
	  segmentL = (unsigned int)( integrationStorage->integPixelSegmentsAlongL * ((currentPoint[0]-firstPixelCoordL-pixelLength*(iL-1) ) / pixelLength) ) + 1;
	  segmentW = (unsigned int)( integrationStorage->integPixelSegmentsAlongW * ((currentPoint[1]-firstPixelCoordW-pixelWidth *(iW-1) ) / pixelWidth ) ) + 1;
	  segmentH = (unsigned int)( integrationStorage->integPixelSegmentsAlongH * (abs(currentPoint[2]) / height ) ) + 1;
	  // Thanks to symmetry we can reduce L and W segments (we have to reduce pixels then, too!)
	  if (segmentL > integrationStorage->integPixelSegmentsAlongL/2) 
	    {
	      segmentL = integrationStorage->integPixelSegmentsAlongL - segmentL + 1;  
	      segmentL_reduced = true;
	    };
	  if (segmentW > integrationStorage->integPixelSegmentsAlongW/2)
	    {
	      segmentW = integrationStorage->integPixelSegmentsAlongW - segmentW + 1;  
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
      temp < 1 ? imin = 1 : imin = temp;
      temp = iL + integMaxNumberPixelsAlongL / 2;
      temp > (long int)numberPixelsAlongL ? imax = numberPixelsAlongL : imax = temp;
      temp = iW - integMaxNumberPixelsAlongW / 2;
      temp < 1 ? jmin = 1 : jmin = temp;
      temp = iW + integMaxNumberPixelsAlongW / 2;
      temp > (long int)numberPixelsAlongW ? jmax = numberPixelsAlongW : jmax = temp;
      
      // Loops over important pixels
      // (borders of a layer part taken into account - see above)
      // cout << "iL = " << iL << " iW = " << iW << " imin = " << imin << " imax = " << imax << " jmin = " << jmin << " jmax = " << jmax << endl;

      // I use map<unsigned long long int pixId, double pixCharge> to keep charges collected in pixels.
      // Machine limit for unsigned long int:
      // ULLONG_MAX = 18 446 744 073 709 551 615
      //                           ^           ^
      // pixID = 10^10*i + j - key for the pixel (i,j) which is used in map<> container [(i,j) <-> (L,W)]
      unsigned long long int pixID;

      for (i = imin ; i <= imax ; i++ )
	{
	  // cout << "i = " << i << endl; 
	  // L limits of integral
	  limitsLow[0] = firstPixelCoordL + (i-1)*pixelLength - currentPoint[0];
	  limitsUp[0]  = limitsLow[0] + pixelLength;     

	  for (j = jmin ; j <= jmax ; j++ )
	    {
	      // cout << "j = " << j << endl; 
	      // W limits of integral
	      limitsLow[1] = firstPixelCoordW + (j-1)*pixelWidth - currentPoint[1];
	      limitsUp[1]  = limitsLow[1] + pixelWidth;     

		
	      // Should we use integration-results?
	      if (useIntegrationStorage)
		{

		  // cout << "Integration-results storage is used!" << endl;
		  // Relative integer coordinates of pixel from main pixel
		  unsigned int pixelL, pixelW;
		  // Thanks to symmetry we can reduce L and W pixels indexes. We have to reduce segments simultaneously!
		  pixelL = (long int)(i) - (long int)(iL) + integMaxNumberPixelsAlongL / 2 + 1;
		  if ( segmentL_reduced   &&  i != iL )
		    {
		      pixelL = integMaxNumberPixelsAlongL - pixelL + 1; 
		    }
		  pixelW = (long int)(j) - (long int)(iW) + integMaxNumberPixelsAlongW / 2 + 1;
		  if ( segmentW_reduced   &&  j != iW )
		    {
		      pixelW = integMaxNumberPixelsAlongW - pixelW + 1; 
		    }

		  unsigned long long int integSegID = integSegmentID(segmentL, segmentW, segmentH, pixelL, pixelW);

		  if ( integrationStorage->isResultStored(integSegID) )
		    {
		      gsl_res = integrationStorage->getResult(integSegID);
		    }
		  else
		    {
		      // Integrate
		      gsl_monte_miser_integrate (&gsl_funToIntegrate, limitsLow, limitsUp, 2, gsl_calls, gsl_r, gsl_s, &gsl_res, &gsl_err);
		      // Store integration result
		      integrationStorage->rememberResult(integSegID,gsl_res);
		    };
		}
	      else
		{
		  // Integrate (here no storage)
		  gsl_monte_miser_integrate (&gsl_funToIntegrate, limitsLow, limitsUp, 2, gsl_calls, gsl_r, gsl_s, &gsl_res, &gsl_err);
		};

	      // Pixels Charge Map
	      // "code" of the pixel - it serves as a key in the map container of pixels (relations: i <-> L, j <-> W)
	      pixID = tenTo10*i + j ;
	      // Add contribution to pixelsChargeMap
	      pixelsChargeMap[ pixID ] = pixelsChargeMap[ pixID ] + ( gsl_res * integChargePerStep );
	      //cout << "charge collected in the pixel " << pixID << " = " << pixelsChargeMap[ pixID ] << endl;
	      //cout << "miser= " << gsl_res << " +- " << gsl_err << endl;
	    }
	}
    }


}


void TDSPixelsChargeMap::print(string filename)
{
  ofstream fout(filename.c_str()); 
    
  type_PixelsChargeMap::iterator i;
  for( i = pixelsChargeMap.begin(); i != pixelsChargeMap.end(); i++ )
    {
      fout << (i->first)/tenTo10 << "\t" << (i->first)%tenTo10 << "\t" << (i->second) << endl;
    }
    
}


double TDSPixelsChargeMap::getPixelCharge(unsigned long int indexAlongL, unsigned long int indexAlongW)
{
  unsigned long long int pixID;
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


// Get total charge collected in the whole map
double TDSPixelsChargeMap::getTotalCharge()
{
  type_PixelsChargeMap::iterator i;

  double totalCharge = 0.;
    
  for( i = pixelsChargeMap.begin(); i != pixelsChargeMap.end(); i++ )
    {
      totalCharge += i->second;
    }
  return totalCharge;
}

// Scale charge deposited in map
double TDSPixelsChargeMap::scaleCharge(double scaleFactor)
{
  type_PixelsChargeMap::iterator i;

  double totalCharge = 0.;
    
  for( i = pixelsChargeMap.begin(); i != pixelsChargeMap.end(); i++ )
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

  while( i != pixelsChargeMap.end())
     {
     double charge =  abs(i->second);
     double varCharge;
     if (charge > 1000.) 
       { // assume Gaussian 
       double sigma = sqrt(charge);
       varCharge = double(CLHEP::RandGauss::shoot(charge,sigma));
       }
     else 
       { // assume Poisson
       varCharge = double(CLHEP::RandPoisson::shoot(charge));
       }
    
     if( i->second < 0.)varCharge = -varCharge;

      if (varCharge == 0. && doCleaning) 
        pixelsChargeMap.erase(i++);
      else    
        {
        i->second = varCharge;
        i++;
        }
    }

  return;
}

  // Convert deposited charge to expected charge output

  void TDSPixelsChargeMap::applyGain(double gain, double gainVariation, double noise, double offset)
{
  type_PixelsChargeMap::iterator i;

  for( i = pixelsChargeMap.begin(); i != pixelsChargeMap.end(); i++ )
    {
    double charge =  i->second;

    double varGain = double(CLHEP::RandGauss::shoot(gain,gainVariation));

    double varNoise = double(CLHEP::RandGauss::shoot(offset,noise));

    i->second = varGain*charge + varNoise;
    }

  return;
}


  // Apply threshold cut to all deposits in the map

  void TDSPixelsChargeMap::applyThresholdCut(double  adcThreshold)
{
  // As pixels will be removed we can not do simple for loop

  type_PixelsChargeMap::iterator i = pixelsChargeMap.begin();

  while( i != pixelsChargeMap.end())
    {
    double charge =  i->second;

    if ( (adcThreshold > 0 && charge < adcThreshold) ||
         (adcThreshold < 0 && charge > adcThreshold) )
      pixelsChargeMap.erase(i++);
    else
      i++;
    }

  return;
}


// Get vector of pixels
vector<TDSPixel> TDSPixelsChargeMap::getVectorOfPixels()
{
  type_PixelsChargeMap::iterator i;

  TDSPixel thePixel;
  vector<TDSPixel> vectorOfPixels;
    
  for( i = pixelsChargeMap.begin(); i != pixelsChargeMap.end(); i++ )
    {
      thePixel.indexAlongL = i->first/tenTo10;
      thePixel.indexAlongW = i->first%tenTo10;
      thePixel.coordL = ((double)(thePixel.indexAlongL)-0.5)*pixelLength + firstPixelCoordL;
      thePixel.coordW = ((double)(thePixel.indexAlongW)-0.5)*pixelWidth  + firstPixelCoordW;
      thePixel.charge = i->second;
      vectorOfPixels.push_back(thePixel);
    }
  return vectorOfPixels;
}

// Get rectangular precluster (for clustering cheater)
// maxNumberOfPixels is the maximal number of pixels with deposit larger than ThresholdCharge.
// ThresholdCharge is only for maxNumberOfPixels cut. In precluster can be included pixels with smaller deposits.
TDSPrecluster TDSPixelsChargeMap::getPrecluster(unsigned int rectLength, unsigned int rectWidth, unsigned int maxNumberOfPixels, double ThresholdCharge)
{
  type_PixelsChargeMap::iterator i, icore;

  unsigned int n=0;
    
  // First pixel
  icore = pixelsChargeMap.begin();
  double greatestDeposit = abs(icore->second);

  for( i = pixelsChargeMap.begin(); i != pixelsChargeMap.end(); i++ )
    {
      if (abs(i->second) > abs(ThresholdCharge))
	{
	  n++;
	    
	  if (abs(i->second) > abs(greatestDeposit))
	    {
	      greatestDeposit = i->second;
	      icore = i;
	    }
	}
    }

  TDSPrecluster thePrecluster;

  if ( n == 0  ||  n > maxNumberOfPixels ) 
    {
      // Empty precluster
      return thePrecluster; 
    }
  else
    {
      // Fill precluster
      thePrecluster.empty = false;
      thePrecluster.pixelL = (icore->first)/tenTo10;
      thePrecluster.pixelW = (icore->first)%tenTo10;
      thePrecluster.coordL = ((double)(thePrecluster.pixelL)-0.5)*pixelLength + firstPixelCoordL;
      thePrecluster.coordW = ((double)(thePrecluster.pixelW)-0.5)*pixelWidth  + firstPixelCoordW;
      thePrecluster.rectLength = rectLength;
      thePrecluster.rectWidth  = rectWidth;

      unsigned  long int l, lmin, lmax, w, wmin, wmax;
      long int temp;
      temp = thePrecluster.pixelL-rectLength/2;
      temp < 1 ? lmin = 1 : lmin = temp;
      temp = thePrecluster.pixelL+rectLength/2;
      temp > (long int)numberPixelsAlongL ? lmax = numberPixelsAlongL : lmax = temp;
      temp = thePrecluster.pixelW-rectWidth/2;
      temp < 1 ? wmin = 1 : wmin = temp;
      temp = thePrecluster.pixelW+rectWidth/2;
      temp > (long int)numberPixelsAlongW ? wmax = numberPixelsAlongL : wmax = temp;

      thePrecluster.rectLmin = lmin;
      thePrecluster.rectLmax = lmax;
      thePrecluster.rectWmin = wmin;
      thePrecluster.rectWmax = wmax;

      for (l=lmin; l<=lmax; l++)
	{
	  for (w=wmin; w<=wmax; w++)
	    {
	      thePrecluster.pixelsCharges.push_back(getPixelCharge(l,w));
	    }
	}
	
      return thePrecluster;
    }
}
