cmake_minimum_required(VERSION 2.8.3)
project(expressiongraph)

find_package(catkin REQUIRED)

find_package(orocos_kdl REQUIRED)
find_package(cmake_modules REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS random)
catkin_package(
  INCLUDE_DIRS include
  LIBRARIES ${PROJECT_NAME}
  DEPENDS orocos_kdl EIGEN3
)


include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
  ${orocos_kdl_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS})

# BUILDING AND LINKING LIBRARY
FILE( GLOB EXPRESSIONTREE_SRCS src/[^.]*.cpp src/[^.]*.cxx)

set(EXPRESSIONTREE_SRCS
    src/expressiontree_expressions.cpp  
    src/expressiontree_motionprofiles.cpp    
    src/expressiontree_rotation.cpp  
    src/expressiontree_wrench.cpp
    src/expressiontree_chain.cpp   
    src/expressiontree_frame.cpp        
    src/expressiontree_twist.cpp     
    src/mptrap.cpp
    src/expressiontree_double.cpp  
    src/expressiontree_mimo.cpp         
    src/expressiontree_vector.cpp    
    )

add_library(${PROJECT_NAME} ${EXPRESSIONTREE_SRCS})
target_link_libraries(${PROJECT_NAME} 
  ${catkin_LIBRARIES} ${orocos_kdl_LIBRARIES} ${EIGEN3_LIBRARIES} ${Boost_LIBRARIES})

# BUILDING AND LINKING TESTS
catkin_add_gtest(${PROJECT_NAME}_test tests/expressiongraph_test.cpp) 
target_link_libraries(${PROJECT_NAME}_test
  ${catkin_LIBRARIES} ${orocos_kdl_LIBRARIES} ${PROJECT_NAME} ${EIGEN3_LIBRARIES} ${Boost_LIBRARIES})

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
)

