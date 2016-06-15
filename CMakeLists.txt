project(nurbslib)
cmake_minimum_required(VERSION 2.8)

# C++11 required to build this library
set(CMAKE_CXX_FLAGS "-Wall -std=c++11")

# Compiler-specific C++11 activation.
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(
        COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (NOT (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7))
        message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.7 or greater.")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
else ()
    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()

# Find the Eigen library
INCLUDE_DIRECTORIES("${CMAKE_SOURCE_DIR}/libs/eigen/")

# look for the boost package
find_package(Boost REQUIRED)
if( Boost_FOUND )
  include_directories(${Boost_INCLUDE_DIRS})
  add_definitions( "-DHAS_BOOST" )
endif()

# find the gsl library
FIND_PACKAGE(GSL REQUIRED)
if(GSL_FOUND)
  include_directories(${GSL_INCLUDE_DIRS})
  set(THIS_LINK_LIBS ${THIS_bLINK_LIBS} ${GSL_LIBRARIES})
endif()

FIND_PACKAGE( Trilinos )
IF( Trilinos_FOUND )
    include_directories( ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})
    link_directories(${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})
    set(THIS_LINK_LIBS ${THIS_LINK_LIBS} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
ENDIF()

# look for VTK package
find_package( VTK REQUIRED )
include( ${VTK_USE_FILE} )
set(THIS_LINK_LIBS ${THIS_LINK_LIBS} ${VTK_LIBRARIES})

SET(LIB_INCLUDE_DIR ${nurbslib_SOURCE_DIR}/include/)
INCLUDE_DIRECTORIES(${LIB_INCLUDE_DIR}
                    ${LIB_INCLUDE_DIR}/2D/
                    ${LIB_INCLUDE_DIR}/quadrature/)
SET(libheaders
        ${LIB_INCLUDE_DIR}/algebra.h
        ${LIB_INCLUDE_DIR}/2D/AnalysisElement2D.h
        ${LIB_INCLUDE_DIR}/2D/BSplineSpace2D.h
        ${LIB_INCLUDE_DIR}/2D/Forest2D.h
        ${LIB_INCLUDE_DIR}/2D/Geometry2D.h
        ${LIB_INCLUDE_DIR}/2D/GeometryElement2D.h
        ${LIB_INCLUDE_DIR}/base.h
        ${LIB_INCLUDE_DIR}/AnalysisElement.h
        ${LIB_INCLUDE_DIR}/BSplineSpace.h
        ${LIB_INCLUDE_DIR}/DOFManager.h
	${LIB_INCLUDE_DIR}/emag.h
        ${LIB_INCLUDE_DIR}/Forest.h
        ${LIB_INCLUDE_DIR}/Geometry.h
        ${LIB_INCLUDE_DIR}/GeometryElement.h
        ${LIB_INCLUDE_DIR}/HCElement.h
        ${LIB_INCLUDE_DIR}/HConformingForest.h
        ${LIB_INCLUDE_DIR}/IAssembly.h
        ${LIB_INCLUDE_DIR}/IAssemblyEmag.h
        ${LIB_INCLUDE_DIR}/IAssemblyFastEmag.h
        ${LIB_INCLUDE_DIR}/IAssemblyHDivNorm.h
        ${LIB_INCLUDE_DIR}/IElem.h
        ${LIB_INCLUDE_DIR}/InputDataStructures.h
        ${LIB_INCLUDE_DIR}/quadrature/IBaseIntegrate.h
        ${LIB_INCLUDE_DIR}/quadrature/IEdgePolarIntegrate.h
        ${LIB_INCLUDE_DIR}/quadrature/IEdgeQuadrature.h
        ${LIB_INCLUDE_DIR}/quadrature/IElemIntegrate.h
        ${LIB_INCLUDE_DIR}/quadrature/IEqualQuadrature.h
        ${LIB_INCLUDE_DIR}/quadrature/IGalerkinIntegrate.h
        ${LIB_INCLUDE_DIR}/quadrature/IGaussQuadrature1D.h
        ${LIB_INCLUDE_DIR}/quadrature/IGaussQuadrature1D.h
        ${LIB_INCLUDE_DIR}/quadrature/IPolarIntegrate.h
        ${LIB_INCLUDE_DIR}/quadrature/ISubElem.h
        ${LIB_INCLUDE_DIR}/quadrature/IVertexQuadrature.h
        ${LIB_INCLUDE_DIR}/MaterialParams.h
        ${LIB_INCLUDE_DIR}/MultiForest.h
        ${LIB_INCLUDE_DIR}/MultithreadBEM.h
        ${LIB_INCLUDE_DIR}/NodalElement.h
        ${LIB_INCLUDE_DIR}/NURBSCommon.h
        ${LIB_INCLUDE_DIR}/NURBSCurve.h
        ${LIB_INCLUDE_DIR}/NURBSSurface.h
        ${LIB_INCLUDE_DIR}/Point.h
        ${LIB_INCLUDE_DIR}/Point3D.h
        ${LIB_INCLUDE_DIR}/Point4D.h
        ${LIB_INCLUDE_DIR}/SurfaceFunctor.h
        ${LIB_INCLUDE_DIR}/VectorSolution.h)

set(LIB_SOURCE_DIR ${nurbslib_SOURCE_DIR}/src/)
set(libsrcs
        ${LIB_SOURCE_DIR}/algebra.cpp
	${LIB_SOURCE_DIR}/emag.cpp	
        ${LIB_SOURCE_DIR}/2D/AnalysisElement2D.cpp
        ${LIB_SOURCE_DIR}/2D/BSplineSpace2D.cpp
        ${LIB_SOURCE_DIR}/2D/Forest2D.cpp
        ${LIB_SOURCE_DIR}/2D/Geometry2D.cpp
        ${LIB_SOURCE_DIR}/2D/GeometryElement2D.cpp
        ${LIB_SOURCE_DIR}/GeometryElement.cpp
        ${LIB_SOURCE_DIR}/Geometry.cpp
        ${LIB_SOURCE_DIR}/HCElement.cpp
        ${LIB_SOURCE_DIR}/HConformingForest.cpp
        ${LIB_SOURCE_DIR}/IAssemblyEmag.cpp
        ${LIB_SOURCE_DIR}/IAssemblyFastEmag.cpp
        ${LIB_SOURCE_DIR}/IAssemblyHDivNorm.cpp
        ${LIB_SOURCE_DIR}/quadrature/IEdgePolarIntegrate.cpp
        ${LIB_SOURCE_DIR}/quadrature/IEdgeQuadrature.cpp
        ${LIB_SOURCE_DIR}/quadrature/IElemIntegrate.cpp
        ${LIB_SOURCE_DIR}/quadrature/IEqualQuadrature.cpp
        ${LIB_SOURCE_DIR}/quadrature/IGalerkinIntegrate.cpp
        ${LIB_SOURCE_DIR}/quadrature/IGaussQuadrature1D.cpp
        ${LIB_SOURCE_DIR}/quadrature/IRegularQuadrature.cpp
        ${LIB_SOURCE_DIR}/quadrature/IPolarIntegrate.cpp
        ${LIB_SOURCE_DIR}/quadrature/ISubElem.cpp
        ${LIB_SOURCE_DIR}/quadrature/IVertexQuadrature.cpp
        ${LIB_SOURCE_DIR}/Forest.cpp
	    ${LIB_SOURCE_DIR}/Point3D.cpp
	    ${LIB_SOURCE_DIR}/Point4D.cpp	
	    ${LIB_SOURCE_DIR}/Point.cpp
	    ${LIB_SOURCE_DIR}/MaterialParams.cpp
	    ${LIB_SOURCE_DIR}/MultiForest.cpp
	    ${LIB_SOURCE_DIR}/NURBSCurve.cpp
	    ${LIB_SOURCE_DIR}/NURBSCommon.cpp
	    ${LIB_SOURCE_DIR}/BSplineSpace.cpp
	    ${LIB_SOURCE_DIR}/NURBSSurface.cpp
	    ${LIB_SOURCE_DIR}/base.cpp
	    ${LIB_SOURCE_DIR}/VectorSolution.cpp
            ${LIB_SOURCE_DIR}/SurfaceFunctor.cpp
        ${libheaders})

# add the external libraries that are contained in the submodule directory
ADD_SUBDIRECTORY(${nurbslib_SOURCE_DIR}/lib/bemquad)
INCLUDE_DIRECTORIES(${nurbslib_SOURCE_DIR}/lib/bemquad/include)
set(THIS_LINK_LIBS ${THIS_LINK_LIBS} bemquad)

add_library(nurbs ${libsrcs})
target_link_libraries(nurbs ${THIS_LINK_LIBS})

INCLUDE(CTest)
ADD_SUBDIRECTORY(test)
ENABLE_TESTING()


