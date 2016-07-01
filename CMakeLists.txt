PROJECT(NuConvex)

#Required in DGtal
cmake_minimum_required (VERSION 2.6) 

FIND_PACKAGE(DGtal 0.6 REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

# SET(Boost_USE_STATIC_LIBS   ON)
# SET(Boost_USE_MULTITHREADED ON)
# SET(Boost_USE_STATIC_RUNTIME OFF)
SET(Boost_FOUND FALSE)
FIND_PACKAGE(Boost 1.46.0 REQUIRED COMPONENTS program_options)
if ( Boost_FOUND )
  message(STATUS "Boost and boost_program_options found.")
  include_directories( ${Boost_INCLUDE_DIRS} )
  SET(BOOST_LIBRARIES
     ${Boost_LIBRAIRIES}  
     ${Boost_PROGRAM_OPTIONS_LIBRARY})
   SET(DGtalLibInc ${Boost_INCLUDE_DIRS})
endif( Boost_FOUND )
message(STATUS ${BOOST_LIBRARIES})

ADD_EXECUTABLE(convexcover convex-cover)
TARGET_LINK_LIBRARIES(convexcover ${DGTAL_LIBRARIES})

ADD_EXECUTABLE(convexnormals convex-normals)
TARGET_LINK_LIBRARIES(convexnormals ${DGTAL_LIBRARIES})

ADD_EXECUTABLE(testMetricCluster testMetricCluster)
TARGET_LINK_LIBRARIES(testMetricCluster ${DGTAL_LIBRARIES})

ADD_EXECUTABLE(nu-convex-viewer nu-convex-viewer)
TARGET_LINK_LIBRARIES(nu-convex-viewer ${BOOST_LIBRARIES} ${DGTAL_LIBRARIES})
