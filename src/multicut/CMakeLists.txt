SET(SOURCE_FILES 
   multicut_srmp.cpp 
   multicut_higher_order.cpp 
   multicut_higher_order_odd_wheel.cpp 
   multicut_higher_order_bicycle_wheel.cpp 
   multicut_opengm_srmp_cycle.cpp 
   multicut_opengm_srmp_cycle_odd_wheel.cpp
   multicut_opengm_srmp_cycle_odd_wheel_odd_bicycle_wheel.cpp
   lifted_multicut_text.cpp
   lifted_multicut_h5.cpp
   lifted_multicut_h5_grid.cpp
#   potts_via_multiway_cut.cpp
   )


# note: HDF5 is only required for the opengm versions. Do not require HDF5 per default
foreach( source_file ${SOURCE_FILES} )
   string( REPLACE ".cpp" "" executable_file ${source_file} )
   add_executable( ${executable_file} ${source_file} ${headers} ${sources})
   target_link_libraries( ${executable_file} m stdc++ pthread LP_MP-Cut ${HDF5_LIBRARIES} ${HDF5_CXX_LIBRARIES} DD_ILP lingeling)
endforeach( source_file ${SOURCE_FILES} )


