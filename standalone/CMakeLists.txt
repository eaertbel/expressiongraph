cmake_minimum_required(VERSION 3.1)
project(expressiongraph
        VERSION 0.0.1
        DESCRIPTION "library for automatic differentiation of robotics-related types such as transformation matrices, rotation matrices, scalars, wrenches and twists")

find_package(orocos_kdl REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS random)
include(GoogleTest)

enable_testing()


include_directories(
  ../include
  ${EIGEN3_INCLUDE_DIRS}
  ${orocos_kdl_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
)

option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)

set(EXPRESSIONGRAPH_SRCS
    ../src/expressiontree_quaternion.cpp  
    ../src/expressiontree_expressions.cpp  
    ../src/expressiontree_motionprofiles.cpp    
    ../src/expressiontree_rotation.cpp  
    ../src/expressiontree_wrench.cpp
    ../src/expressiontree_chain.cpp   
    ../src/expressiontree_frame.cpp        
    ../src/expressiontree_twist.cpp     
    ../src/mptrap.cpp
    ../src/expressiontree_double.cpp  
    ../src/expressiontree_mimo.cpp         
    ../src/expressiontree_vector.cpp    
    )

#IF(BUILD_SHARED_LIBS)
#    add_library(${PROJECT_NAME} SHARED ${EXPRESSIONGRAPH_SRCS})
#    set_target_properties(${PROJECT_NAME} PROPERTIES SOVERSION 0)   # major version of the shared library file.
#    set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${PROJECT_VERSION})
#    target_link_libraries(${PROJECT_NAME} ${orocos_kdl_LIBRARIES} ${EIGEN3_LIBRARIES} ${Boost_LIBRARIES})
#ENDIF(BUILD_SHARED_LIBS)
add_library(${PROJECT_NAME} ${EXPRESSIONGRAPH_SRCS})
target_link_libraries(${PROJECT_NAME} ${orocos_kdl_LIBRARIES} ${EIGEN3_LIBRARIES} ${Boost_LIBRARIES})

configure_file(expressiongraph.pc.in expressiongraph.pc @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/expressiongraph.pc DESTINATION lib/pkgconfig)


# BUILDING AND LINKING TESTS
#add_executable(expressiongraph_test ${PROJECT_SOURCE_DIR}/tests/expressiongraph_test.cpp)
#gtest_add_tests(TARGET expressiongraph_test) 
#target_link_libraries(expressiongraph_test ${orocos_kdl_LIBRARIES} ${PROJECT_NAME} ${EIGEN3_LIBRARIES} ${Boost_LIBRARIES})

# POTENTIALLY, BUILDING AND LINKING EXAMPLES
OPTION(ENABLE_EXAMPLES "enable compilation of a series of examples" ON) 
OPTION(ENABLE_TESTS "enable compilation of unit tests" ON) 
#INCLUDE(../examples/CMakeLists.txt)
#ADD_SUBDIRECTORY(../examples)

IF(ENABLE_EXAMPLES)

    set(EXAMPLES
            conversions_example
            debug1
            expressiontree_cached_named
            expressiontree_conditional
            expressiontree_distance_to_line
            expressiontree_example
            expressiontree_example10
            expressiontree_example11
            expressiontree_example12
            expressiontree_example2
            expressiontree_example3
            expressiontree_example4
            expressiontree_example5
            expressiontree_example6
            expressiontree_example7
            expressiontree_example8
            expressiontree_example9
            expressiontree_mimo_ex
            expressiontree_motion_profile
            expressiontree_motion_profile2
            expressiontree_optimizer
            expressiontree_perpendicular
            expressiontree_sensor
            expressiontree_variabletype
            geometry
            initial_value
            matrix_traits
            mptrap_tst
            resolved_motion_rate_control
            saturate
            solving_and_cloning
            tutorial1
            tutorial2
            tutorial3
            quat
    )

     set(BDIR ..)

     add_executable(solving_and_cloning ${BDIR}/examples/solving_and_cloning.cpp )
     TARGET_LINK_LIBRARIES(solving_and_cloning ${PROJECT_NAME} ${Eigen_LIBRARIES})


     add_executable(mptrap_tst ${BDIR}/examples/mptrap_tst.cpp )
     TARGET_LINK_LIBRARIES(mptrap_tst ${PROJECT_NAME} ${Eigen_LIBRARIES})


     add_executable(matrix_traits ${BDIR}/examples/matrix_traits.cpp )
     TARGET_LINK_LIBRARIES(matrix_traits ${PROJECT_NAME} ${Eigen_LIBRARIES})

     add_executable(geometry ${BDIR}/examples/geometry.cpp )
     TARGET_LINK_LIBRARIES(geometry ${PROJECT_NAME} ${Eigen_LIBRARIES})

     add_executable(expressiontree_variabletype ${BDIR}/examples/expressiontree_variabletype.cpp )
     TARGET_LINK_LIBRARIES(expressiontree_variabletype ${PROJECT_NAME} ${Eigen_LIBRARIES})


     add_executable(expressiontree_example ${BDIR}/examples/expressiontree_example.cpp )
     TARGET_LINK_LIBRARIES(expressiontree_example ${PROJECT_NAME} ${Eigen_LIBRARIES})
     
     add_executable(expressiontree_example2 ${BDIR}/examples/expressiontree_example2.cpp )
     TARGET_LINK_LIBRARIES(expressiontree_example2 ${PROJECT_NAME} ${Eigen_LIBRARIES})
     
     add_executable(expressiontree_example3 ${BDIR}/examples/expressiontree_example3.cpp )
     TARGET_LINK_LIBRARIES(expressiontree_example3 ${PROJECT_NAME} ${Eigen_LIBRARIES})
     
     add_executable(expressiontree_example4 ${BDIR}/examples/expressiontree_example4.cpp )
     TARGET_LINK_LIBRARIES(expressiontree_example4 ${PROJECT_NAME} ${Eigen_LIBRARIES})
     
     add_executable(expressiontree_example5 ${BDIR}/examples/expressiontree_example5.cpp )
     TARGET_LINK_LIBRARIES(expressiontree_example5 ${PROJECT_NAME} ${Eigen_LIBRARIES})
     
     add_executable(expressiontree_example6 ${BDIR}/examples/expressiontree_example6.cpp )
     TARGET_LINK_LIBRARIES(expressiontree_example6 ${PROJECT_NAME} ${Eigen_LIBRARIES})
     
     add_executable(expressiontree_example7 ${BDIR}/examples/expressiontree_example7.cpp )
     TARGET_LINK_LIBRARIES(expressiontree_example7 ${PROJECT_NAME} ${Eigen_LIBRARIES})
     
     add_executable(expressiontree_example8 ${BDIR}/examples/expressiontree_example8.cpp )
     TARGET_LINK_LIBRARIES(expressiontree_example8 ${PROJECT_NAME} ${Eigen_LIBRARIES})
     
     add_executable(expressiontree_example9 ${BDIR}/examples/expressiontree_example9.cpp )
     TARGET_LINK_LIBRARIES(expressiontree_example9 ${PROJECT_NAME} ${Eigen_LIBRARIES}) 
     add_executable(expressiontree_example10 ${BDIR}/examples/expressiontree_example10.cpp )
     TARGET_LINK_LIBRARIES(expressiontree_example10 ${PROJECT_NAME} ${Eigen_LIBRARIES})
     
     add_executable(expressiontree_example11 ${BDIR}/examples/expressiontree_example11.cpp )
     TARGET_LINK_LIBRARIES(expressiontree_example11 ${PROJECT_NAME} ${Eigen_LIBRARIES})

     add_executable(expressiontree_example12 ${BDIR}/examples/expressiontree_example12.cpp )
     TARGET_LINK_LIBRARIES(expressiontree_example12 ${PROJECT_NAME} ${Eigen_LIBRARIES})
     
     #add_executable(expressiontree_example13 examples/expressiontree_example13.cpp )
     #TARGET_LINK_LIBRARIES(expressiontree_example13 ${PROJECT_NAME} ${Eigen_LIBRARIES})
     
     add_executable(resolved_motion_rate_control ${BDIR}/examples/resolved_motion_rate_control.cpp )
     TARGET_LINK_LIBRARIES(resolved_motion_rate_control  ${PROJECT_NAME} ${Eigen_LIBRARIES})
     
     add_executable(expressiontree_perpendicular ${BDIR}/examples/expressiontree_perpendicular.cpp )
     TARGET_LINK_LIBRARIES(expressiontree_perpendicular  ${PROJECT_NAME} ${Eigen_LIBRARIES})
     
     add_executable(conversions_example ${BDIR}/examples/conversions_example.cpp )
     TARGET_LINK_LIBRARIES(conversions_example ${PROJECT_NAME} ${Eigen_LIBRARIES})
     
     add_executable(expressiontree_distance_to_line  ${BDIR}/examples/expressiontree_distance_to_line.cpp)
     TARGET_LINK_LIBRARIES(expressiontree_distance_to_line ${PROJECT_NAME} ${Eigen_LIBRARIES})
     
     add_executable(expressiontree_conditional  ${BDIR}/examples/expressiontree_conditional.cpp)
     TARGET_LINK_LIBRARIES(expressiontree_conditional ${PROJECT_NAME} ${Eigen_LIBRARIES})

     add_executable(expressiontree_cached_named ${BDIR}/examples/expressiontree_cached_named.cpp)
     TARGET_LINK_LIBRARIES(expressiontree_cached_named ${PROJECT_NAME} ${Eigen_LIBRARIES})

     add_executable(expressiontree_optimizer ${BDIR}/examples/expressiontree_optimizer.cpp)
     TARGET_LINK_LIBRARIES(expressiontree_optimizer ${PROJECT_NAME} ${Eigen_LIBRARIES})

     add_executable(expressiontree_mimo_ex ${BDIR}/examples/expressiontree_mimo_ex.cpp)
     TARGET_LINK_LIBRARIES(expressiontree_mimo_ex ${PROJECT_NAME} ${Eigen_LIBRARIES})

     add_executable(expressiontree_sensor ${BDIR}/examples/expressiontree_sensor.cpp)
     TARGET_LINK_LIBRARIES(expressiontree_sensor ${PROJECT_NAME} ${Eigen_LIBRARIES})

     add_executable(expressiontree_motion_profile ${BDIR}/examples/expressiontree_motion_profile.cpp)
     TARGET_LINK_LIBRARIES(expressiontree_motion_profile ${PROJECT_NAME} ${Eigen_LIBRARIES})

     add_executable(expressiontree_motion_profile2 ${BDIR}/examples/expressiontree_motion_profile2.cpp)
     TARGET_LINK_LIBRARIES(expressiontree_motion_profile2 ${PROJECT_NAME} ${Eigen_LIBRARIES})


     add_executable(initial_value ${BDIR}/examples/initial_value.cpp)
     TARGET_LINK_LIBRARIES(initial_value ${PROJECT_NAME} ${Eigen_LIBRARIES})




    # TUTORIAL FILES:
     add_executable(tutorial1 ${BDIR}/examples/tutorial1.cpp)
     TARGET_LINK_LIBRARIES(tutorial1 ${PROJECT_NAME} ${Eigen_LIBRARIES})

     add_executable(tutorial2 ${BDIR}/examples/tutorial2.cpp)
     TARGET_LINK_LIBRARIES(tutorial2 ${PROJECT_NAME} ${Eigen_LIBRARIES})

     add_executable(tutorial3 ${BDIR}/examples/tutorial3.cpp)
     TARGET_LINK_LIBRARIES(tutorial3 ${PROJECT_NAME} ${Eigen_LIBRARIES})

     add_executable(saturate ${BDIR}/examples/saturate.cpp)
     TARGET_LINK_LIBRARIES(saturate ${PROJECT_NAME} ${Eigen_LIBRARIES})

     add_executable(debug1 ${BDIR}/examples/debug1.cpp)
     TARGET_LINK_LIBRARIES(debug1 ${PROJECT_NAME} ${Eigen_LIBRARIES})


     # DEVELOPMENT FILES:

     add_executable(quat ${BDIR}/examples/quat.cpp)
     TARGET_LINK_LIBRARIES(quat ${PROJECT_NAME} ${Eigen_LIBRARIES})


ENDIF(ENABLE_EXAMPLES)  

# INSTALLING LIBRARY AND HEADER-FILES
install(TARGETS ${PROJECT_NAME} DESTINATION lib)

install(TARGETS ${EXAMPLES} DESTINATION share/doc/libexpressiongraph/examples)

install(DIRECTORY ../include/
  DESTINATION include 
  FILES_MATCHING PATTERN "*.hpp"
)

install(DIRECTORY ../examples/
  DESTINATION share/doc/libexpressiongraph/examples
  FILES_MATCHING PATTERN "*.cpp"
)

install(FILES ../LICENSE ../license_EUPL_v1_1_en.pdf ../COPYING ../README.md
  DESTINATION share/doc/libexpressiongraph
)

###########################################################################################################
# Testing
###########################################################################################################
# Download and unpack googletest at configure time

configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories("${gtest_SOURCE_DIR}/include")
endif()



#IF(ENABLE_TESTS)
#   add_executable(expressiongraph_test ../tests/expressiongraph_test.cpp )
#   TARGET_LINK_LIBRARIES(expressiongraph_test ${PROJECT_NAME} ${Eigen_LIBRARIES})
#ENDIF(ENABLE_TESTS)


# Now simply link against gtest or gtest_main as needed. Eg
#add_executable(expressiongraph_test ../tests/expressiongraph_test.cpp )
add_executable(test_expressiongraph ../tests/test_expressiongraph.cpp )
#target_link_libraries(expressiongraph_test gtest_main ${PROJECT_NAME} ${Eigen_LIBRARIES})
target_link_libraries(test_expressiongraph gtest_main ${PROJECT_NAME} ${Eigen_LIBRARIES})
gtest_discover_tests(test_expressiongraph)
#add_test(NAME quat_test COMMAND quat_test)

