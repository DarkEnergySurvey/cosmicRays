#include(FindPackageHandleStandardArgs.cmake)

#=============================================================================
# If the user has provided ``NDARRAY_ROOT_DIR``, use it!  Choose items found
# at this location over system locations.
if( EXISTS "$ENV{NDARRAY_ROOT_DIR}" )
  file( TO_CMAKE_PATH "$ENV{NDARRAY_ROOT_DIR}" NDARRAY_ROOT_DIR )
  set( NDARRAY_ROOT_DIR "${NDARRAY_ROOT_DIR}" CACHE PATH "Prefix for NDARRAY installation." )
endif()
if( NOT EXISTS "${NDARRAY_ROOT_DIR}" )
  set( NDARRAY_USE_PKGCONFIG ON )
endif()

#=============================================================================

# As a first try, use the PkgConfig module.  This will work on many
# *NIX systems.  See :module:`findpkgconfig`
# This will return ``NDARRAY_INCLUDEDIR`` and ``NDARRAY_LIBDIR`` used below.
if( NDARRAY_USE_PKGCONFIG )
  find_package(PkgConfig)
  pkg_check_modules( NDARRAY3 QUIET NDARRAY3 )

  if( EXISTS "${NDARRAY_INCLUDEDIR}" )
    get_filename_component( NDARRAY_ROOT_DIR "${NDARRAY_INCLUDEDIR}" DIRECTORY CACHE)
  endif()
endif()

#=============================================================================
# Set NDARRAY_INCLUDE_DIRS and NDARRAY_LIBRARIES. If we skipped the PkgConfig step, try
# to find the libraries at $NDARRAY_ROOT_DIR (if provided) or in standard system
# locations.  These find_library and find_path calls will prefer custom
# locations over standard locations (HINTS).  If the requested file is not found
# at the HINTS location, standard system locations will be still be searched
# (/usr/lib64 (Redhat), lib/i386-linux-gnu (Debian)).

find_path( NDARRAY_INCLUDE_DIR
  NAMES ndarray.h
  HINTS ${NDARRAY_ROOT_DIR}/include ${NDARRAY_INCLUDEDIR}
)

#=============================================================================
# handle the QUIETLY and REQUIRED arguments and set NDARRAY_FOUND to TRUE if all
# listed variables are TRUE
find_package_handle_standard_args( NDARRAY
  FOUND_VAR
    NDARRAY_FOUND
  REQUIRED_VARS
    NDARRAY_INCLUDE_DIR
    )

mark_as_advanced( NDARRAY_ROOT_DIR NDARRAY_INCLUDE_DIR )

if(NDARRAY_FOUND)
    set( NDARRAY_INCLUDE_DIRS ${NDARRAY_INCLUDE_DIR} )
    if(NOT TARGET NDARRAY)
        add_library(NDARRAY UNKNOWN IMPORTED)
        set_target_properties(NDARRAY PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${NDARRAY_INCLUDE_DIR}")
    endif()
endif()
