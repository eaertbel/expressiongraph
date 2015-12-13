cmake_minimum_required(VERSION 2.8.3)
project(expressiongraph)

find_package(catkin REQUIRED)

find_package(orocos_kdl REQUIRED)
find_package(cmake_modules REQUIRED)
find_package(Eigen REQUIRED)
include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${Eigen_INCLUDE_DIRS}
  ${orocos_kdl_INCLUDE_DIRS})

# BUILDING AND LINKING LIBRARY
FILE( GLOB EXPRESSIONTREE_SRCS src/[^.]*.cpp src/[^.]*.cxx)
add_library(${PROJECT_NAME} ${EXPRESSIONTREE_SRCS})
target_link_libraries(${PROJECT_NAME} 
  ${catkin_LIBRARIES} ${orocos_kdl_LIBRARIES} ${Eigen_LIBRARIES})

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES ${PROJECT_NAME}
  DEPENDS orocos_kdl Eigen)


# BUILDING AND LINKING TESTS
catkin_add_gtest(${PROJECT_NAME}_test tests/expressiongraph_test.cpp) 
target_link_libraries(${PROJECT_NAME}_test
  ${catkin_LIBRARIES} ${orocos_kdl_LIBRARIES} ${PROJECT_NAME} ${Eigen_LIBRARIES})

# POTENTIALLY, BUILDING AND LINKING EXAMPLES
OPTION(ENABLE_EXAMPLES "enable compilation of a series of examples" ON) 
INCLUDE(${PROJECT_SOURCE_DIR}/examples/CMakeLists.txt)

# INSTALLING LIBRARY AND HEADER-FILES
install(TARGETS ${PROJECT_NAME} ${EXAMPLES}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION})

install(DIRECTORY include/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
  FILES_MATCHING PATTERN "*.hpp"
  PATTERN ".svn" EXCLUDE)

