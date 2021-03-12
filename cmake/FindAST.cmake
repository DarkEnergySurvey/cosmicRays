#include(${CMAKE_MODULE_PATH}/FindPackageHandleStandardArgs.cmake)

#=============================================================================
# If the user has provided ``AST_ROOT_DIR``, use it!  Choose items found
# at this location over system locations.
if( EXISTS "$ENV{AST_ROOT_DIR}" )
  file( TO_CMAKE_PATH "$ENV{AST_ROOT_DIR}" AST_ROOT_DIR )
  set( AST_ROOT_DIR "${AST_ROOT_DIR}" CACHE PATH "Prefix for AST installation." )
endif()
if( NOT EXISTS "${AST_ROOT_DIR}" )
  set( AST_USE_PKGCONFIG ON )
endif()

#=============================================================================

# As a first try, use the PkgConfig module.  This will work on many
# *NIX systems.  See :module:`findpkgconfig`
# This will return ``AST_INCLUDEDIR`` and ``AST_LIBDIR`` used below.
if( AST_USE_PKGCONFIG )
  find_package(PkgConfig)
  pkg_check_modules( AST QUIET ast )

  if( EXISTS "${AST_INCLUDEDIR}" )
    get_filename_component( AST_ROOT_DIR "${AST_INCLUDEDIR}" DIRECTORY CACHE)
  endif()
endif()

#=============================================================================
# Set AST_INCLUDE_DIRS and AST_LIBRARIES. If we skipped the PkgConfig step, try
# to find the libraries at $AST_ROOT_DIR (if provided) or in standard system
# locations.  These find_library and find_path calls will prefer custom
# locations over standard locations (HINTS).  If the requested file is not found
# at the HINTS location, standard system locations will be still be searched
# (/usr/lib64 (Redhat), lib/i386-linux-gnu (Debian)).

find_path( AST_INCLUDE_DIR
  NAMES star/ast.h
  HINTS ${AST_ROOT_DIR}/include ${AST_INCLUDEDIR}
)
find_library( AST_LIBRARY
  NAMES ast
  HINTS ${AST_ROOT_DIR}/lib ${AST_LIBDIR}
  PATH_SUFFIXES Release Debug
)
find_library( AST_ERR_LIBRARY
  NAMES ast_err
  HINTS ${AST_ROOT_DIR}/lib ${AST_LIBDIR}
  PATH_SUFFIXES Release Debug
)
find_library( AST_GRF3D_LIBRARY
  NAMES ast_grf3d
  HINTS ${AST_ROOT_DIR}/lib ${AST_LIBDIR}
  PATH_SUFFIXES Release Debug
)
find_library( AST_GRF2_LIBRARY
  NAMES ast_grf_2.0
  HINTS ${AST_ROOT_DIR}/lib ${AST_LIBDIR}
  PATH_SUFFIXES Release Debug
)
find_library( AST_GRF3_LIBRARY
  NAMES ast_grf_3.2
  HINTS ${AST_ROOT_DIR}/lib ${AST_LIBDIR}
  PATH_SUFFIXES Release Debug
)
find_library( AST_GRF5_LIBRARY
  NAMES ast_grf_5.6
  HINTS ${AST_ROOT_DIR}/lib ${AST_LIBDIR}
  PATH_SUFFIXES Release Debug
)

#=============================================================================
# handle the QUIETLY and REQUIRED arguments and set AST_FOUND to TRUE if all
# listed variables are TRUE
find_package_handle_standard_args( AST
  FOUND_VAR
    AST_FOUND
  REQUIRED_VARS
    AST_INCLUDE_DIR
    AST_LIBRARY
    AST_ERR_LIBRARY
    AST_GRF3D_LIBRARY
    AST_GRF2_LIBRARY
    AST_GRF3_LIBRARY
    AST_GRF5_LIBRARY
    )

mark_as_advanced( AST_ROOT_DIR AST_LIBRARY AST_INCLUDE_DIR
  AST_ERR_LIBRARY AST_GRF3D_LIBRARY AST_GRF2_LIBRARY AST_GRF3_LIBRARY
    AST_GRF5_LIBRARY
 )

if(AST_FOUND)
    set( AST_INCLUDE_DIRS ${AST_INCLUDE_DIR} )
    set( AST_LIBRARIES ${AST_LIBRARY} ${AST_CBLAS_LIBRARY} ${AST_ERR_LIBRARY} ${AST_GRF3D_LIBRARY} ${AST_GRF2_LIBRARY} ${AST_GRF3_LIBRARY} ${AST_GRF5_LIBRARY})
    if(NOT TARGET Star::AST)
        add_library(Star::AST UNKNOWN IMPORTED)
        set_target_properties(Star::AST PROPERTIES
            IMPORTED_LOCATION        "${AST_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${AST_INCLUDE_DIR}")
    endif()
endif()
