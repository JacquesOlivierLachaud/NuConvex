PROJECT(NuConvex)

#Required in DGtal
cmake_minimum_required (VERSION 2.6) 

FIND_PACKAGE(DGtal 0.6 REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

ADD_EXECUTABLE(convexcover convex-cover)
TARGET_LINK_LIBRARIES(convexcover ${DGTAL_LIBRARIES})

ADD_EXECUTABLE(convexnormals convex-normals)
TARGET_LINK_LIBRARIES(convexnormals ${DGTAL_LIBRARIES})

ADD_EXECUTABLE(testMetricCluster testMetricCluster)
TARGET_LINK_LIBRARIES(testMetricCluster ${DGTAL_LIBRARIES})

ADD_EXECUTABLE(nu-convex-viewer nu-convex-viewer)
TARGET_LINK_LIBRARIES(nu-convex-viewer ${DGTAL_LIBRARIES})
