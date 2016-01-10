IF(ENABLE_EXAMPLES)

  set(EXAMPLES
    geometry
    expressiontree_example
    expressiontree_example2
    expressiontree_example3
    expressiontree_example4
    expressiontree_example5
    expressiontree_example6
    expressiontree_example7
    expressiontree_example8
    expressiontree_example9
    expressiontree_example10
    expressiontree_example11
    expressiontree_example12
    resolved_motion_rate_control
    expressiontree_perpendicular
    conversions_example
    expressiontree_distance_to_line
    expressiontree_conditional
    expressiontree_cached_named
    expressiontree_optimizer
    expressiontree_mimo_ex
    expressiontree_sensor
    expressiontree_motion_profile
    tutorial1
    tutorial2
    tutorial3
    saturate
    debug1
    matrix_traits
    initial_value
 )

 add_executable(matrix_traits examples/matrix_traits.cpp )
 TARGET_LINK_LIBRARIES(matrix_traits ${PROJECT_NAME} ${Eigen_LIBRARIES})

 add_executable(geometry examples/geometry.cpp )
 TARGET_LINK_LIBRARIES(geometry ${PROJECT_NAME} ${Eigen_LIBRARIES})

 add_executable(expressiontree_example examples/expressiontree_example.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_example ${PROJECT_NAME} ${Eigen_LIBRARIES})
 
 add_executable(expressiontree_example2 examples/expressiontree_example2.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_example2 ${PROJECT_NAME} ${Eigen_LIBRARIES})
 
 add_executable(expressiontree_example3 examples/expressiontree_example3.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_example3 ${PROJECT_NAME} ${Eigen_LIBRARIES})
 
 add_executable(expressiontree_example4 examples/expressiontree_example4.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_example4 ${PROJECT_NAME} ${Eigen_LIBRARIES})
 
 add_executable(expressiontree_example5 examples/expressiontree_example5.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_example5 ${PROJECT_NAME} ${Eigen_LIBRARIES})
 
 add_executable(expressiontree_example6 examples/expressiontree_example6.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_example6 ${PROJECT_NAME} ${Eigen_LIBRARIES})
 
 add_executable(expressiontree_example7 examples/expressiontree_example7.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_example7 ${PROJECT_NAME} ${Eigen_LIBRARIES})
 
 add_executable(expressiontree_example8 examples/expressiontree_example8.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_example8 ${PROJECT_NAME} ${Eigen_LIBRARIES})
 
 add_executable(expressiontree_example9 examples/expressiontree_example9.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_example9 ${PROJECT_NAME} ${Eigen_LIBRARIES}) 
 add_executable(expressiontree_example10 examples/expressiontree_example10.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_example10 ${PROJECT_NAME} ${Eigen_LIBRARIES})
 
 add_executable(expressiontree_example11 examples/expressiontree_example11.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_example11 ${PROJECT_NAME} ${Eigen_LIBRARIES})

 add_executable(expressiontree_example12 examples/expressiontree_example12.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_example12 ${PROJECT_NAME} ${Eigen_LIBRARIES})
 
 #add_executable(expressiontree_example13 examples/expressiontree_example13.cpp )
 #TARGET_LINK_LIBRARIES(expressiontree_example13 ${PROJECT_NAME} ${Eigen_LIBRARIES})
 
 add_executable(resolved_motion_rate_control examples/resolved_motion_rate_control.cpp )
 TARGET_LINK_LIBRARIES(resolved_motion_rate_control  ${PROJECT_NAME} ${Eigen_LIBRARIES})
 
 add_executable(expressiontree_perpendicular examples/expressiontree_perpendicular.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_perpendicular  ${PROJECT_NAME} ${Eigen_LIBRARIES})
 
 add_executable(conversions_example examples/conversions_example.cpp )
 TARGET_LINK_LIBRARIES(conversions_example ${PROJECT_NAME} ${Eigen_LIBRARIES})
 
 add_executable(expressiontree_distance_to_line  examples/expressiontree_distance_to_line.cpp)
 TARGET_LINK_LIBRARIES(expressiontree_distance_to_line ${PROJECT_NAME} ${Eigen_LIBRARIES})
 
 add_executable(expressiontree_conditional  examples/expressiontree_conditional.cpp)
 TARGET_LINK_LIBRARIES(expressiontree_conditional ${PROJECT_NAME} ${Eigen_LIBRARIES})

 add_executable(expressiontree_cached_named examples/expressiontree_cached_named.cpp)
 TARGET_LINK_LIBRARIES(expressiontree_cached_named ${PROJECT_NAME} ${Eigen_LIBRARIES})

 add_executable(expressiontree_optimizer examples/expressiontree_optimizer.cpp)
 TARGET_LINK_LIBRARIES(expressiontree_optimizer ${PROJECT_NAME} ${Eigen_LIBRARIES})

 add_executable(expressiontree_mimo_ex examples/expressiontree_mimo_ex.cpp)
 TARGET_LINK_LIBRARIES(expressiontree_mimo_ex ${PROJECT_NAME} ${Eigen_LIBRARIES})

 add_executable(expressiontree_sensor examples/expressiontree_sensor.cpp)
 TARGET_LINK_LIBRARIES(expressiontree_sensor ${PROJECT_NAME} ${Eigen_LIBRARIES})

 add_executable(expressiontree_motion_profile examples/expressiontree_motion_profile.cpp)
 TARGET_LINK_LIBRARIES(expressiontree_motion_profile ${PROJECT_NAME} ${Eigen_LIBRARIES})

 add_executable(initial_value examples/initial_value.cpp)
 TARGET_LINK_LIBRARIES(initial_value ${PROJECT_NAME} ${Eigen_LIBRARIES})



# TUTORIAL FILES:
 add_executable(tutorial1 examples/tutorial1.cpp)
 TARGET_LINK_LIBRARIES(tutorial1 ${PROJECT_NAME} ${Eigen_LIBRARIES})

 add_executable(tutorial2 examples/tutorial2.cpp)
 TARGET_LINK_LIBRARIES(tutorial2 ${PROJECT_NAME} ${Eigen_LIBRARIES})

 add_executable(tutorial3 examples/tutorial3.cpp)
 TARGET_LINK_LIBRARIES(tutorial3 ${PROJECT_NAME} ${Eigen_LIBRARIES})

 add_executable(saturate examples/saturate.cpp)
 TARGET_LINK_LIBRARIES(saturate ${PROJECT_NAME} ${Eigen_LIBRARIES})

 add_executable(debug1 examples/debug1.cpp)
 TARGET_LINK_LIBRARIES(debug1 ${PROJECT_NAME} ${Eigen_LIBRARIES})


ENDIF(ENABLE_EXAMPLES)  

