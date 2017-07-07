########################################################################
# cmake module for finding G4
#
# returns:
#   G4_FOUND         : set to TRUE or FALSE
#   G4_INCLUDE_DIRS  : paths to G4 includes
#   G4_LIBRARY_DIRS  : paths to G4 libraries
#   G4_LIBRARIES     : list of G4 libraries (either set to
#       libCGAPack - if found  - OR - G4 libraries in $G4/tmp )
#   G4_CGA_FOUND     : set to TRUE or FALSE (checks for libCGAPack)
#
# @author Jan Engels, DESY
########################################################################

  SET ( G4_DIR "$ENV{G4INSTALL}" )
  
  find_path(G4_INCLUDE_DIRS NAMES G4Cache.hh 
    PATHS
    ${G4_DIR}/include/Geant4/
  )  

 message ( "${G4_INCLUDE_DIRS}" )  
 FILE( GLOB G4_LIBRARIES "${G4_DIR}/lib*/libG4*.so" )

# ---------- final checking ---------------------------------------------------

  INCLUDE( FindPackageHandleStandardArgs )

# set G4_FOUND to TRUE if all listed variables are TRUE and not empty
  FIND_PACKAGE_HANDLE_STANDARD_ARGS( G4 DEFAULT_MSG G4_INCLUDE_DIRS G4_LIBRARIES )
 


