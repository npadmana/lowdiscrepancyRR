set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpicxx")

# Simple mangle tests
add_executable(test_mangle test_mangle.cpp mangle.cpp)

# Copy files
configure_file(bossmask_template.cfg ${CMAKE_CURRENT_BINARY_DIR}/bossmask_template.cfg COPYONLY)

add_library(bosslib SHARED mangle.cpp BossMask.cpp ang2d.cpp)
target_link_libraries(bosslib npgsl m boost_program_options boost_filesystem boost_system)

# Test codes 
add_executable(test_program_options test_program_options.cpp)
target_link_libraries(test_program_options boost_program_options)


# Now start building codes
add_executable(test_bossmask test_bossmask.cpp)
target_link_libraries(test_bossmask bosslib)

add_executable(bossmask_area bossmask_area.cpp)
target_link_libraries(bossmask_area bosslib)

add_executable(bossRR bossRR.cpp)
target_link_libraries(bossRR bosslib)

add_executable(test_partition test_partition.cpp)