IF(ENABLE_EXAMPLES)

 add_executable(geometry examples/geometry.cpp )
 TARGET_LINK_LIBRARIES(geometry ${PROJECT_NAME})
 

 add_executable(expressiontree_example examples/expressiontree_example.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_example ${PROJECT_NAME})
 
 add_executable(expressiontree_example2 examples/expressiontree_example2.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_example2 ${PROJECT_NAME})
 
 add_executable(expressiontree_example3 examples/expressiontree_example3.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_example3 ${PROJECT_NAME})
 
 add_executable(expressiontree_example4 examples/expressiontree_example4.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_example4 ${PROJECT_NAME})
 
 add_executable(expressiontree_example5 examples/expressiontree_example5.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_example5 ${PROJECT_NAME})
 
 add_executable(expressiontree_example6 examples/expressiontree_example6.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_example6 ${PROJECT_NAME})
 
 add_executable(expressiontree_example7 examples/expressiontree_example7.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_example7 ${PROJECT_NAME})
 
 add_executable(expressiontree_example8 examples/expressiontree_example8.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_example8 ${PROJECT_NAME})
 
 add_executable(expressiontree_example9 examples/expressiontree_example9.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_example9 ${PROJECT_NAME})
 
 add_executable(expressiontree_example10 examples/expressiontree_example10.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_example10 ${PROJECT_NAME})
 
 add_executable(expressiontree_example11 examples/expressiontree_example11.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_example11 ${PROJECT_NAME})

 add_executable(expressiontree_example12 examples/expressiontree_example12.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_example12 ${PROJECT_NAME})
 
 #add_executable(expressiontree_example13 examples/expressiontree_example13.cpp )
 #TARGET_LINK_LIBRARIES(expressiontree_example13 ${PROJECT_NAME})
 
 add_executable(resolved_motion_rate_control examples/resolved_motion_rate_control.cpp )
 TARGET_LINK_LIBRARIES(resolved_motion_rate_control  ${PROJECT_NAME})
 
 add_executable(expressiontree_perpendicular examples/expressiontree_perpendicular.cpp )
 TARGET_LINK_LIBRARIES(expressiontree_perpendicular  ${PROJECT_NAME})
 
 add_executable(conversions_example examples/conversions_example.cpp )
 TARGET_LINK_LIBRARIES(conversions_example ${PROJECT_NAME})
 
 add_executable(expressiontree_distance_to_line  examples/expressiontree_distance_to_line.cpp)
 TARGET_LINK_LIBRARIES(expressiontree_distance_to_line ${PROJECT_NAME})
 
 add_executable(expressiontree_conditional  examples/expressiontree_conditional.cpp)
 TARGET_LINK_LIBRARIES(expressiontree_conditional ${PROJECT_NAME})

 add_executable(expressiontree_cached_named examples/expressiontree_cached_named.cpp)
 TARGET_LINK_LIBRARIES(expressiontree_cached_named ${PROJECT_NAME})

 add_executable(expressiontree_optimizer examples/expressiontree_optimizer.cpp)
 TARGET_LINK_LIBRARIES(expressiontree_optimizer ${PROJECT_NAME})

 add_executable(expressiontree_mimo_ex examples/expressiontree_mimo_ex.cpp)
 TARGET_LINK_LIBRARIES(expressiontree_mimo_ex ${PROJECT_NAME})

 add_executable(expressiontree_sensor examples/expressiontree_sensor.cpp)
 TARGET_LINK_LIBRARIES(expressiontree_sensor ${PROJECT_NAME})

 add_executable(expressiontree_motion_profile examples/expressiontree_motion_profile.cpp)
 TARGET_LINK_LIBRARIES(expressiontree_motion_profile ${PROJECT_NAME})



# TUTORIAL FILES:
 add_executable(tutorial1 examples/tutorial1.cpp)
 TARGET_LINK_LIBRARIES(tutorial1 ${PROJECT_NAME})

 add_executable(tutorial2 examples/tutorial2.cpp)
 TARGET_LINK_LIBRARIES(tutorial2 ${PROJECT_NAME})

 add_executable(tutorial3 examples/tutorial3.cpp)
 TARGET_LINK_LIBRARIES(tutorial3 ${PROJECT_NAME})

 add_executable(saturate examples/saturate.cpp)
 TARGET_LINK_LIBRARIES(saturate ${PROJECT_NAME})

 add_executable(debug1 examples/debug1.cpp)
 TARGET_LINK_LIBRARIES(debug1 ${PROJECT_NAME})




ENDIF(ENABLE_EXAMPLES)  

