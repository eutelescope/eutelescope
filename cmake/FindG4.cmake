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

#    MESSAGE ( STATUS "looking for G4: [${G4_DIR}] " )    

    EXECUTE_PROCESS(
        COMMAND find ${G4_DIR}/include/Geant4/ -type d # -printf "%p;" # not recognized on mac osx
        OUTPUT_VARIABLE G4_INCLUDE_DIRS
        RESULT_VARIABLE _exit_code
    )
 
#   EXECUTE_PROCESS(
#       COMMAND "ls ${G4_DIR}/lib/*so"
#        OUTPUT_VARIABLE G4_LIBRARIES
#        RESULT_VARIABLE _exit_code
#    )

    SET( G4_LIBRARIES 
"/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4analysis.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4physicslists.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4run.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4persistency.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4readout.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4clhep.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4digits_hits.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4error_propagation.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4event.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4FR.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4geometry.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4gl2ps.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4global.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4GMocren.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4graphics_reps.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4intercoms.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4interfaces.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4materials.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4modeling.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4OpenGL.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4OpenInventor.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4parmodels.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4particles.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4processes.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4RayTracer.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4tracking.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4track.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4Tree.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4visHepRep.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4vis_management.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4visXXX.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4VRML.so;/home/ilcsoft/v01-15/geant4/9.5.p01/lib/libG4zlib.so")


    # libraries
    FIND_LIBRARY( G4_LIBRARIES
        NAMES G4
        PATHS "${G4_DIR}" 
        PATH_SUFFIXES lib lib/$ENV{G4SYSTEM}
    )





# ---------- final checking ---------------------------------------------------
INCLUDE( FindPackageHandleStandardArgs )
# set G4_FOUND to TRUE if all listed variables are TRUE and not empty
FIND_PACKAGE_HANDLE_STANDARD_ARGS( G4 DEFAULT_MSG G4_INCLUDE_DIRS G4_LIBRARIES )


SET( G4_FOUND ${G4_FOUND} )

#MESSAGE( STATUS "G4_FOUND: ${G4_FOUND}")
#MESSAGE( STATUS "G4_DIR: ${G4_DIR}")
#MESSAGE( STATUS "G4_INCS: ${G4_INCLUDE_DIRS}")
#MESSAGE( STATUS "G4_LIB: ${G4_LIBRARIES}")





