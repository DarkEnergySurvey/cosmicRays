cmake_minimum_required(VERSION 3.5)

project(cosmicRays LANGUAGES CXX)

set(PACKAGE_NAME "cosmicRays")
set(PACKAGE_STRING "cosmicRays 1.0.0")

set(PACKAGE_VERSION_MAJOR "1")
set(PACKAGE_VERSION_MINOR "0")
set(PACKAGE_VERSION_PATCH "0")
set(PACKAGE_VERSION_LABEL "")

if(PACKAGE_VERSION_LABEL)
  set(PACKAGE_VERSION "${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}-${PACKAGE_VERSION_LABEL}")
else()
  set(PACKAGE_VERSION "${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}")
endif()


set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Prohibit in-source build
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "In-source build prohibited.")
endif("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/cmake")
set(Boost_USE_RELEASE_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
find_package(Boost 1.36.0 REQUIRED COMPONENTS regex iostreams filesystem)
find_package(OpenMP REQUIRED)
find_package(GSL REQUIRED)
find_package(Cfitsio REQUIRED)
find_package(Minuit REQUIRED)
find_package(AST REQUIRED)
find_package(FFTW REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Python3 3.7.0 REQUIRED Development Interpreter)
find_package(pybind11 CONFIG REQUIRED HINTS "${Python3_SITELIB}/pybind11/share/cmake/pybind11")
find_package(NDARRAY REQUIRED)

include_directories(cosmicRays "${PROJECT_SOURCE_DIR}/include")
include_directories(cosmicRays ${Boost_INCLUDE_DIR})
include_directories(cosmicRays ${OpenMP_CXX_INCLUDE_DIR})
include_directories(cosmicRays ${GSL_INCLUDE_DIR})
include_directories(cosmicRays ${Cfitsio_INCLUDE_DIR})
include_directories(cosmicRays ${Minuit_INCLUDE_DIR})
include_directories(cosmicRays ${AST_INCLUDE_DIR})
include_directories(cosmicRays ${FFTW_INCLUDES})
include_directories(cosmicRays ${EIGEN3_INCLUDE_DIR})
include_directories(cosmicRays ${Python3_INCLUDE_DIRS})
include_directories(cosmicRays ${NDARRAY_INCLUDE_DIR})

file(GLOB_RECURSE SOURCES "src/*.cc")
file(GLOB_RECURSE INCLUDES "include/*.h*")

add_library(cosmicRays SHARED
    ${SOURCES}
    ${INCLUDES}
)
target_link_libraries(cosmicRays ${Boost_LIBRARIES})
message(STATUS "BO ${Boost_LIBRARIES}")
target_link_libraries(cosmicRays ${OpenMP_CXX_LIBRARIES})
message(STATUS "OM ${OpenMP_CXX_LIBRARIES}")
target_link_libraries(cosmicRays ${GSL_LIBRARIES})
message(STATUS "GS ${GSL_LIBRARIES}")
target_link_libraries(cosmicRays ${Minuit_LIBRARIES})
message(STATUS "MI ${Minuit_LIBRARIES}")
target_link_libraries(cosmicRays ${Cfitsio_LIBRARIES})
message(STATUS "CF ${Cfitsio_LIBRARIES}")
target_link_libraries(cosmicRays ${AST_LIBRARIES})
message(STATUS "AS ${AST_LIBRARIES}")
target_link_libraries(cosmicRays ${FFTW_LIBRARIES})
message(STATUS "FF ${FFTW_LIBRARIES}")
target_link_libraries(cosmicRays pybind11::module)
target_compile_options(cosmicRays PUBLIC ${OpenMP_CXX_FLAGS})
message(STATUS "")
target_link_options(cosmicRays PUBLIC  "-fopenmp")
target_compile_definitions(cosmicRays PRIVATE COSMICRAYS_LIBRARY)

install(TARGETS cosmicRays DESTINATION lib)
