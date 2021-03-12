#include(${CMAKE_MODULE_PATH}/FindPackageHandleStandardArgs.cmake)

#=============================================================================
# If the user has provided ``Minuit_ROOT_DIR``, use it!  Choose items found
# at this location over system locations.
if( EXISTS "$ENV{Minuit_ROOT_DIR}" )
  file( TO_CMAKE_PATH "$ENV{Minuit_ROOT_DIR}" Minuit_ROOT_DIR )
  set( Minuit_ROOT_DIR "${Minuit_ROOT_DIR}" CACHE PATH "Prefix for Minuit installation." )
endif()
if( NOT EXISTS "${Minuit_ROOT_DIR}" )
  set( Minuit_USE_PKGCONFIG ON )
endif()

#=============================================================================

# As a first try, use the PkgConfig module.  This will work on many
# *NIX systems.  See :module:`findpkgconfig`
# This will return ``Minuit_INCLUDEDIR`` and ``Minuit_LIBDIR`` used below.
if( Minuit_USE_PKGCONFIG )
  find_package(PkgConfig)
  pkg_check_modules( Minuit QUIET Minuit )

  if( EXISTS "${Minuit_INCLUDEDIR}" )
    get_filename_component( Minuit_ROOT_DIR "${Minuit_INCLUDEDIR}" DIRECTORY CACHE)
  endif()
endif()

#=============================================================================
# Set Minuit_INCLUDE_DIRS and Minuit_LIBRARIES. If we skipped the PkgConfig step, try
# to find the libraries at $Minuit_ROOT_DIR (if provided) or in standard system
# locations.  These find_library and find_path calls will prefer custom
# locations over standard locations (HINTS).  If the requested file is not found
# at the HINTS location, standard system locations will be still be searched
# (/usr/lib64 (Redhat), lib/i386-linux-gnu (Debian)).

find_path( Minuit_INCLUDE_DIR
  NAMES Minuit2/MnMinos.h
  HINTS ${Minuit_ROOT_DIR}/include ${Minuit_INCLUDEDIR}
)
find_library( Minuit_LIBRARY
  NAMES Minuit2
  HINTS ${Minuit_ROOT_DIR}/lib ${Minuit_LIBDIR}
  PATH_SUFFIXES Release Debug
)

#=============================================================================
# handle the QUIETLY and REQUIRED arguments and set Minuit_FOUND to TRUE if all
# listed variables are TRUE
find_package_handle_standard_args( Minuit
  FOUND_VAR
    Minuit_FOUND
  REQUIRED_VARS
    Minuit_INCLUDE_DIR
    Minuit_LIBRARY
    )

mark_as_advanced( Minuit_ROOT_DIR Minuit_LIBRARY Minuit_INCLUDE_DIR
 )

if(Minuit_FOUND)
    set( Minuit_INCLUDE_DIRS ${Minuit_INCLUDE_DIR} )
    set( Minuit_LIBRARIES ${Minuit_LIBRARY})
    if(NOT TARGET Minuit)
        add_library(Minuit UNKNOWN IMPORTED)
        set_target_properties(Minuit PROPERTIES
            IMPORTED_LOCATION        "${Minuit_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${Minuit_INCLUDE_DIR}")
    endif()
endif()
