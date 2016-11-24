project(confnurbs)
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
    add_definitions ("-Wno-unused-local-typedefs")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
else ()
    message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif ()


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
set(THIS_LINK_LIBS ${THIS_LINK_LIBS} ${GSL_LIBRARIES})
endif()

# look for VTK package
find_package( VTK REQUIRED )
include( ${VTK_USE_FILE} )
set(THIS_LINK_LIBS ${THIS_LINK_LIBS} ${VTK_LIBRARIES})

SET(LIB_INCLUDE_DIR ${confnurbs_SOURCE_DIR}/include/)
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
        ${LIB_INCLUDE_DIR}/BezierNodalElement.h
        ${LIB_INCLUDE_DIR}/BezierVectorElement.h
        ${LIB_INCLUDE_DIR}/BoundingBoxIterator.h
        ${LIB_INCLUDE_DIR}/BSplineSpace.h
        ${LIB_INCLUDE_DIR}/DOFManager.h
        ${LIB_INCLUDE_DIR}/Forest.h
        ${LIB_INCLUDE_DIR}/Geometry.h
        ${LIB_INCLUDE_DIR}/GeometryElement.h
        ${LIB_INCLUDE_DIR}/HCElement.h
        ${LIB_INCLUDE_DIR}/HConformingForest.h
        ${LIB_INCLUDE_DIR}/IElem.h
        ${LIB_INCLUDE_DIR}/InputDataStructures.h
        ${LIB_INCLUDE_DIR}/quadrature/IBaseIntegrate.h
        ${LIB_INCLUDE_DIR}/quadrature/IEdgePolarIntegrate.h
        ${LIB_INCLUDE_DIR}/quadrature/IEdgeQuadrature.h
        ${LIB_INCLUDE_DIR}/quadrature/IElemIntegrate.h
        ${LIB_INCLUDE_DIR}/quadrature/ISubElemIntegrate.h
        ${LIB_INCLUDE_DIR}/quadrature/IEqualQuadrature.h
        ${LIB_INCLUDE_DIR}/quadrature/IEqualQuadratureTri.h
        ${LIB_INCLUDE_DIR}/quadrature/IGalerkinIntegrate.h
        ${LIB_INCLUDE_DIR}/quadrature/IGaussQuadrature1D.h
        ${LIB_INCLUDE_DIR}/quadrature/IGaussQuadrature1D.h
        ${LIB_INCLUDE_DIR}/quadrature/IRegularQuadrature.h
        ${LIB_INCLUDE_DIR}/quadrature/IPolarIntegrate.h
        ${LIB_INCLUDE_DIR}/quadrature/IPolarDegenerate.h
        ${LIB_INCLUDE_DIR}/quadrature/ISubElem.h
        ${LIB_INCLUDE_DIR}/quadrature/IVertexQuadrature.h
        ${LIB_INCLUDE_DIR}/MultiForest.h
        ${LIB_INCLUDE_DIR}/NedelecVectorElement.h
        ${LIB_INCLUDE_DIR}/NodalElement.h
        ${LIB_INCLUDE_DIR}/Norm.h
        ${LIB_INCLUDE_DIR}/NURBSCommon.h
        ${LIB_INCLUDE_DIR}/NURBSCurve.h
        ${LIB_INCLUDE_DIR}/NURBSSurface.h
        ${LIB_INCLUDE_DIR}/OutputVTK.h
        ${LIB_INCLUDE_DIR}/Point.h
        ${LIB_INCLUDE_DIR}/Point3D.h
        ${LIB_INCLUDE_DIR}/Point4D.h)

set(LIB_SOURCE_DIR ${confnurbs_SOURCE_DIR}/src/)
set(libsrcs
        ${LIB_SOURCE_DIR}/algebra.cpp
        ${LIB_SOURCE_DIR}/2D/AnalysisElement2D.cpp
        ${LIB_SOURCE_DIR}/2D/BSplineSpace2D.cpp
        ${LIB_SOURCE_DIR}/2D/Forest2D.cpp
        ${LIB_SOURCE_DIR}/2D/Geometry2D.cpp
        ${LIB_SOURCE_DIR}/2D/GeometryElement2D.cpp
	    ${LIB_SOURCE_DIR}/base.cpp
        ${LIB_SOURCE_DIR}/BezierNodalElement.cpp
        ${LIB_SOURCE_DIR}/BezierVectorElement.cpp
        ${LIB_SOURCE_DIR}/BoundingBoxIterator.cpp
        ${LIB_SOURCE_DIR}/GeometryElement.cpp
        ${LIB_SOURCE_DIR}/Geometry.cpp
        ${LIB_SOURCE_DIR}/HCElement.cpp
        ${LIB_SOURCE_DIR}/HConformingForest.cpp
        ${LIB_SOURCE_DIR}/Forest.cpp
	    ${LIB_SOURCE_DIR}/Point3D.cpp
	    ${LIB_SOURCE_DIR}/Point4D.cpp
	    ${LIB_SOURCE_DIR}/Point.cpp
        ${LIB_SOURCE_DIR}/NedelecVectorElement.cpp
	    ${LIB_SOURCE_DIR}/MultiForest.cpp
	    ${LIB_SOURCE_DIR}/Norm.cpp
	    ${LIB_SOURCE_DIR}/NURBSCurve.cpp
	    ${LIB_SOURCE_DIR}/NURBSCommon.cpp
	    ${LIB_SOURCE_DIR}/BSplineSpace.cpp
	    ${LIB_SOURCE_DIR}/NURBSSurface.cpp
	    ${LIB_SOURCE_DIR}/NURBSSurface.cpp
	    ${LIB_SOURCE_DIR}/OutputVTK.cpp
        ${LIB_SOURCE_DIR}/quadrature/IEdgePolarIntegrate.cpp
        ${LIB_SOURCE_DIR}/quadrature/IEdgeQuadrature.cpp
        ${LIB_SOURCE_DIR}/quadrature/IElemIntegrate.cpp
        ${LIB_SOURCE_DIR}/quadrature/IEqualQuadrature.cpp
        ${LIB_SOURCE_DIR}/quadrature/IEqualQuadratureTri.cpp
        ${LIB_SOURCE_DIR}/quadrature/IGalerkinIntegrate.cpp
        ${LIB_SOURCE_DIR}/quadrature/IGaussQuadrature1D.cpp
        ${LIB_SOURCE_DIR}/quadrature/IRegularQuadrature.cpp
        ${LIB_SOURCE_DIR}/quadrature/IPolarIntegrate.cpp
        ${LIB_SOURCE_DIR}/quadrature/IPolarDegenerate.cpp
        ${LIB_SOURCE_DIR}/quadrature/ISubElem.cpp
        ${LIB_SOURCE_DIR}/quadrature/IVertexQuadrature.cpp
        ${libheaders})

add_library(confnurbs ${libsrcs})
target_link_libraries(confnurbs ${THIS_LINK_LIBS})

INCLUDE(CTest)
ADD_SUBDIRECTORY(test)
ENABLE_TESTING()
