// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
// root includes
#include <TFile.h>
#include <TList.h>
#include <TKey.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TPaveLabel.h>
#include <TPaveStats.h>
#include <TDirectoryFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TProfile.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include <TCutG.h>
#include <TSystem.h>
#include <TAxis.h>

// system includes
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <set>

#ifdef DOXY
//! The histogram namespace
/*! This namespace contains all the global variables and function used
 *  to display and print standard reference histograms from the ROOT
 *  files genereated by the AIDAProcessor. For more information see
 *  \ref histogramPackage
 */
namespace eutelHistogram {
#endif

  /*! \page histogramPackage The histogram package
   *
   *  The great majority of EUTelescope processors are generating
   *  control histograms if the AIDAProcessor is added in the steering
   *  file and if a suitable AIDA implementation is installed into the
   *  system.
   *
   *  In case RAIDA (The ROOT based implementation of AIDA) is used,
   *  then all histograms will be saved into a ROOT file divided into
   *  subfolder named after the processor which generated them. For
   *  example if in the steering file there is the
   *  EUTelEtaCalculatorProcessor and the user named it "ETA", then all
   *  the histograms will be save into the ROOT file within a subfolder
   *  called "ETA".
   *
   *  To check these histograms, the user has to open the histogram file
   *  in an interactive ROOT session, browse the file and plot
   *  histograms according to his/her needs. To make the user life
   *  easier, we provided also a series of macros contained in the \a
   *  eutelHistogram namespace to display the standard and most commonly
   *  used reference histograms. All the canvases displayed on the video
   *  are also saved as pictures into a dedicated folder.
   *
   *  \section histoHowTo How to use the package
   *
   *  Using the eutelescopeHistos package should be as simple as
   *  starting an interactive ROOT session. The only thing the user has
   *  to do is to copy the following three files: \a eutelescopeHistos.h
   *  \a eutelescopeHistos.C and \a rootlogon.C in the folder where is
   *  doing the analysis.
   *  When starting ROOT will execute the ROOT logon and on the fly
   *  generating a shared library (\a eutelescopeHistos.so) that is soon
   *  after loaded in ROOT as a plugin.
   *  A list of all the available functions is shown together with a
   *  very brief help.
   *
   *  To use the library, the first thing one has to do is to set the
   *  run name. This is absolutely <b>optional</b> and the only result
   *  it will produce is that the run name/number will appear on all the
   *  canvases so that they can be directly used for a presentation or
   *  even a paper. In case the user doesn't set any name, the Unknown
   *  Run label will be used instead.
   *
   *  After that the run name is set, the user can play with \e show
   *  function to display histograms. All these functions have been
   *  coded in order to be as general as possible, in the meaning that
   *  they will automatically guess how many sensors there are in the
   *  the run and consequently the number of canvases that should be
   *  produced. The package is assuming that the histogram the user
   *  wants to see are stored in a sub folder with a standard name. If
   *  the user in the Marlin steering file used a different name, then
   *  the package will complain about the not found folder. Since the
   *  folder name is a parameter, the user can type, for example,
   *  \code
   *  root [0] correlatorFolderName = "MyCorrelatorFolder"
   *  (class TString)"MyCorrelatorFolder"
   *  root [1] showCorrelationPlot("correlation-histo-file.root")
   *  \endcode
   *  to change the default value.
   *  The same is true for the output format of the canvas pictures. The
   *  default value is ".png" and can be changed to any other supported
   *  format typing:
   *  \code
   *  root [3] pictureOutputFormat = ".gif"
   *  \endcode
   *
   *  \section histoIssue Know Issue
   *  A possible problem that can arise from the use of the package is
   *  caused by a bug of RAIDA that in some circumstances is generating
   *  two instances of the same object. This is severely affecting the
   *  capability of the package to properly guess the number of sensors
   *  in the setup. Unfortunately there is no way to know in advance if
   *  this issue will affect your system. In case you see a segmentation
   *  fault and an unexpected number of canvases, then this could be an
   *  indication that something is corrupted in the input ROOT file. To
   *  fix this issue, you have to follow this very simple procedure:
   *
   *  <b>Create an empty ROOT file</b>. Use the commands below in a ROOT
   *  shell:
   *  \code
   *  TFile * file = TFile::Open("empty.root","RECREATE");
   *  file->Write();
   *  file->Close();
   *  \endcode
   *
   *  <b>Add this empty file to your corrupted ROOT file</b>. This
   *  procedure will remove the duplicated instances of objects and fix
   *  the number of sensor guess. Use the following commands from a
   *  command shell (outside ROOT).
   *  \code
   *  hadd myFixedFile.root myCorruptedFile.root empty.root
   *  \endcode
   *
   *  Now try again to execute the show command from ROOT
   */

// global variables
  //! Global variable containing the run name
  TString runName = "Unknown run";

  //! Global variable with the name of pedestal processor folder
  TString pedeProcessorFolderName = "PedestalAndNoiseCalculator";

  //! Global variable with the cluster histogramming processor folder
  TString clusterHistoFolderName = "FillHisto";

  //! Global variable with the Eta calculator folder
  TString etaFolderName = "ETA";

  //! Global variable with the Correlator folder
  TString correlatorFolderName = "Correlator";

  //! Global variable with the TestFitter folder
  TString trackerFolderName = "TestFitter";

  //! Global variable for the Alignment folder
  TString milleFolderName = "Align";

  //! Global variable with the picture output format
  TString pictureOutputFormat = ".png";

// real function prototypes
  //! Show cluster plot
  /*! This function can be used to plot in a nice way the standard
   *  control histograms concerning clusters.
   *
   *  The function will open the ROOT file specified by \a filename and
   *  look for the folder containing the cluster histograms. The name of
   *  this folder can be changed by the user using the \a
   *  clusterHistoFolderName global variable.
   *
   *  The following canvases will be plotted and saved as picture files:
   *
   *  \li <b>SNRCanvas</b>. This canvas contains two histograms for each
   *  detector up to three detectors. If the total number of sensors is
   *  exceeding 3, another canvas will be created. The histogram on the
   *  left hand side is Seed pixel signal to noise ratio distribution,
   *  the one on the right hand side is the cluster 3x3 SNR. Both
   *  distribution are drawn along with their Landau fit functions.
   *
   *  \li <b>HitMapCanvas</b>. This canvas contains one 2D histogram for
   *  each detector up to six detectors. If the total number of sensors is
   *  exceeding 6, another canvas will be created. Each histogram
   *  represents the seed pixel position accumulated overall the
   *  events. The color code is also draw.
   *
   *  \li <b>NoiseCanvas</b>. This canvas has not to be confused with
   *  the single pixel noise canvas generated by showPedeNoisePlot. It
   *  draws one histogram per sensor, up to 6 sensors per canvas. Each
   *  histogram represents the cluster noise distribution calculated as
   *  the sum in quadrature of all the pixels belonging to the
   *  clusters. Since not all clusters are made by the same amount of
   *  pixels, having a multi-peaked distribution is normal.
   *
   *  @param filename The ROOT file name.
   *
   */
  void showClusterPlot( const char * filename );



  //! Show pedestal and noise plot
  /*! This function can be used to plot in a nice way the standard
   *  control histograms concerning pedestal and noise.
   *
   *  The function will open the ROOT file specified by \a filename and
   *  look for the folder containing the pedestal, noise and status
   *  histograms. The name of this folder can be changed by the user
   *  using the pedeProcessorFolderName global variable.
   *
   *  The following canvases will be plotted and saved as picture files:
   *
   *  \li <b>PedeCanvas</b> This canvas contains two histograms for each
   *  sensor, up to three sensors per canvas. The 2D histogram on the
   *  left represents the pedestal distribution on the sensors surface,
   *  while the 1D plot is the pedestal distribution.
   *
   *  \li <b>NoiseCanvas</b> Same as the PedeCanvas but now representing
   *  noise figures instead of pedestal values.
   *
   *  \li <b>StatusCanvas</b> This canvas contains one 2D histogram for
   *  each detector up to six sensors per canvas. The red dots
   *  corresponds to masked pixels.
   *
   *  In case the detector is a known one (e.g. a mimotel) then the user
   *  can specified it in the function call to get a more detailed
   *  analysis with a breakdown of the noise channel by channel. To
   *  switch off the advanced analysis simply set \a detector to "".
   *
   *  @bug If the telescope setup is a mixed one (mimotel and mimosa18),
   *  the advance noise analysis will fail.
   *
   *  @param filename The ROOT file name.
   *  @param detector The name of the detector plane
   */
  void showPedeNoisePlot( const char * filename,  const char * detector = "mimotel"  );

  //! Show Eta function plot
  /*! This function can be used to plot in a nice way the standard
   *  control histograms concerning Eta function.
   *
   *  The function will open the ROOT file specified by \a filename and
   *  look for the folder containing the Eta histograms. The name of
   *  this folder can be changed by the user using the
   *  etaProcessorFolderName global variable.
   *
   *  The following canvases will be plotted and saved as picture files:
   *
   *  \li <b>EtaXCanvas</b>. This canvas contains two histograms for
   *  each sensor, up to three sensors per canvas. The plot on left is
   *  the cluster center position as calculated by the Center Of Gravity
   *  algorithm. The plot on the right hand side is the corresponding
   *  Eta function. All these plots are related to the X direction.
   *
   *  \li <b>EtaYCanvas</b>. As the EtaXCanvas, but related to the Y
   *  direction.
   *
   *  \li <b>CoGCanvas</b>. This canvas contains one 2D histograms for
   *  each sensor, up to six sensors per canvas. The histogram
   *  represents the center of gravity position of the cluster center
   *  relative to the central pixel.
   *
   *  \param filename The ROOT file name.
   */
  void showEtaPlot( const char * filename );

  //! Show correlation plot
  /*! This function can be used to plot in a nice way the standard
   *  control histograms concerning Correlator processor.
   *
   *  Depending if the correlation is done on clusters or hits or both,
   *  the function will display one or two series of canvases. For each
   *  canvases, the user can see the 2D correlation plot and the
   *  straight line fit resulting from the projection of this 2D plot
   *  onto the x axis.
   *
   *  In order not to be overwhelmed by the combinatorial background,
   *  the function is automatically guessing a region of interest around
   *  the correlation line and projecting only these values. The ROI is
   *  identified on the 2D plot with a dashed spline.
   *
   *  The correlation histograms are shown for following telescope
   *  plane, it is to say plane N vs plane N+1. So if there are 6 planes
   *  in the setup there will be only 5 correlation histograms with the
   *  corresponding fits.
   */
  void showCorrelationPlot( const char * filename );

  //! Show tracker histogram
  /*! This function is plotting several histograms taken from
   *  TestFitter output. Those are quite useful to understand the
   *  tracker performance and to verify the tracker alignment.
   */
  void showTrackerPlot( const char * filename );

  //! Show Mille histogram
  /*! This function is plotting the output histograms of the Mille
   *  processor.
   */
  void showMillePlot( const char * filename );

  // utility prototypes
  //! Open a ROOT file
  /*! This utility function is used to open a ROOT file, if this is not
   *  yet opened. In case, the file is opened but a Zombie, than it is
   *  closed and re-opened.
   */
  TFile * closeAndReopenFile( const char * filename );

  //! Set the run name for canvas title
  /*! This function can be used to set the global variable \a
   *  runName. This string is automatically added at the canvas title
   *  pad.
   */
  void setRunName(const char * name );

  //! Close canvases
  /*! Before opening a new canvas, the user should call this function in
   *  order to close all the other canvases having the same name to
   *  avoid confusion.
   */
  void closeCanvases( const char * canvasBaseName );

  //! Prepare the output folder
  /*! Each show plot function is saving the canvases as figure. This
   *  function is used to create a folder named after the \a runName and
   *  another folder named according to the kind of plots.
   */
  const char * prepareOutputFolder( const char * type );

  //! Usage
  /*! Simple function to print out which are the functions available.
   */
  void usage();

  //! Set default size for histogram axis
  void setDefaultAxis(TAxis * axis );

  //! Check the folder existence
  /*! The AIDA processor is creating in every ROOT file a folder named
   *  after the name of the processor generating the histograms. Every
   *  users can in principle change the name of the processor without
   *  changing the functionality. This function is checking among a
   *  list of possible names for a folder in the input file. If none
   *  of them is found, then an error message is shown and it returns
   *  NULL. If more than one is found, a list of good
   *  folders is presented to the user who has to take a decision.
   *
   *  @param candidateNames A list of all names to look for
   *  @param inputFile A pointer to the ROOT input file name
   *  @return A pointer to the folder
   */
  TDirectoryFile * checkFolder( std::vector< string >& candidateNames, TFile * inputFile );

  //! Convert to string
  /*! This template function is used to convert to string any type
   *  having a streamer. This utility can be used to "sum" numbers to
   *  string, for example.
   */
  template<typename T>
  std::string toString(T data ) {
    std::stringstream ss;
    ss << data ;
    return ss.str();

  }

#ifdef DOXY
}
#endif

//  LocalWords:  filename pedeProcessorFolderName PedeCanvas NoiseCanvas png
//  LocalWords:  StatusCanvas mimotel etaProcessorFolderName iostream sstream
//  LocalWords:  PedestalAndNoiseCalculator FillHisto histogramming Correlator
//  LocalWords:  TestFitter
