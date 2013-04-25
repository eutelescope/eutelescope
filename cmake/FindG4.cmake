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

  EXECUTE_PROCESS(
      COMMAND test -d ${G4_DIR}/include/Geant4/
      OUTPUT_VARIABLE G4_INCLUDE_DIRS
      RESULT_VARIABLE _exit_code
  )

  FILE( GLOB G4_LIBRARIES "${G4_DIR}/lib/libG4*.so" )


# ---------- final checking ---------------------------------------------------

  INCLUDE( FindPackageHandleStandardArgs )

# set G4_FOUND to TRUE if all listed variables are TRUE and not empty
  FIND_PACKAGE_HANDLE_STANDARD_ARGS( G4 DEFAULT_MSG G4_INCLUDE_DIRS G4_LIBRARIES )
 


