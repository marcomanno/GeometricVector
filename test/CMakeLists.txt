#include (ACGCommon)
 cmake_minimum_required(VERSION 3.5)
include_directories (
  ${CMAKE_CURRENT_SOURCE_DIR}/../include
  ${CMAKE_CURRENT_SOURCE_DIR}/../Catch/include
#  C:/Users/marco/OneDrive/Documents/PROJECTS/ThirdParties/boost_1_61_0
)

# Create new target named unittests_hexmeshing
FILE(GLOB UNITTEST_CC *.cc)
FILE(GLOB UNITTEST_HH *.hh )
FILE(GLOB VECTOR_H "../include/*.h")

# Create unittest executable
add_executable(unittests ${UNITTEST_CC} ${UNITTEST_HH} ${VECTOR_H})

# Set output directory to ${BINARY_DIR}/Unittests
set (OUTPUT_DIR "${CMAKE_BINARY_DIR}/Unittests")
set_target_properties(unittests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIR}) 

#target_link_libraries(unittests Base PolygonTriangularization Geo Boolean Import Topology Utils) 
