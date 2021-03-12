#include(${CMAKE_MODULE_PATH}/FindPackageHandleStandardArgs.cmake)

#=============================================================================
# If the user has provided ``Cfitsio_ROOT_DIR``, use it!  Choose items found
# at this location over system locations.
if( EXISTS "$ENV{Cfitsio_ROOT_DIR}" )
  file( TO_CMAKE_PATH "$ENV{Cfitsio_ROOT_DIR}" Cfitsio_ROOT_DIR )
  set( Cfitsio_ROOT_DIR "${Cfitsio_ROOT_DIR}" CACHE PATH "Prefix for Cfitsio installation." )
endif()
if( NOT EXISTS "${Cfitsio_ROOT_DIR}" )
  set( Cfitsio_USE_PKGCONFIG ON )
endif()

#=============================================================================

# As a first try, use the PkgConfig module.  This will work on many
# *NIX systems.  See :module:`findpkgconfig`
# This will return ``Cfitsio_INCLUDEDIR`` and ``Cfitsio_LIBDIR`` used below.
if( Cfitsio_USE_PKGCONFIG )
  find_package(PkgConfig)
  pkg_check_modules( Cfitsio QUIET Cfitsio )

  if( EXISTS "${Cfitsio_INCLUDEDIR}" )
    get_filename_component( Cfitsio_ROOT_DIR "${Cfitsio_INCLUDEDIR}" DIRECTORY CACHE)
  endif()
endif()

#=============================================================================
# Set Cfitsio_INCLUDE_DIRS and Cfitsio_LIBRARIES. If we skipped the PkgConfig step, try
# to find the libraries at $Cfitsio_ROOT_DIR (if provided) or in standard system
# locations.  These find_library and find_path calls will prefer custom
# locations over standard locations (HINTS).  If the requested file is not found
# at the HINTS location, standard system locations will be still be searched
# (/usr/lib64 (Redhat), lib/i386-linux-gnu (Debian)).

find_path( Cfitsio_INCLUDE_DIR
  NAMES fitsio.h
  HINTS ${Cfitsio_ROOT_DIR}/include ${Cfitsio_INCLUDEDIR}
)
find_library( Cfitsio_LIBRARY
  NAMES cfitsio
  HINTS ${Cfitsio_ROOT_DIR}/lib ${Cfitsio_LIBDIR}
  PATH_SUFFIXES Release Debug
)

#=============================================================================
# handle the QUIETLY and REQUIRED arguments and set Cfitsio_FOUND to TRUE if all
# listed variables are TRUE
find_package_handle_standard_args( Cfitsio
  FOUND_VAR
    Cfitsio_FOUND
  REQUIRED_VARS
    Cfitsio_INCLUDE_DIR
    Cfitsio_LIBRARY
    )

mark_as_advanced( Cfitsio_ROOT_DIR Cfitsio_LIBRARY Cfitsio_INCLUDE_DIR
 )

if(Cfitsio_FOUND)
    set( Cfitsio_INCLUDE_DIRS ${Cfitsio_INCLUDE_DIR} )
    set( Cfitsio_LIBRARIES ${Cfitsio_LIBRARY})
    if(NOT TARGET Cfitsio)
        add_library(Cfitsio UNKNOWN IMPORTED)
        set_target_properties(Cfitsio PROPERTIES
            IMPORTED_LOCATION        "${Cfitsio_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${Cfitsio_INCLUDE_DIR}")
    endif()
endif()
