
# Test data files
set (testdata hemisphere.cfg hemisphere.ply hemisphere_input.cfg sphere.cfg sphere.ply sphere_input.cfg lowz.cfg lowz.ply lowz.random)
foreach(file1 ${testdata})
  configure_file(${file1} ${CMAKE_CURRENT_BINARY_DIR}/${file1} COPYONLY)
endforeach(file1)


# Build tests
set (testlist partition_test area_test RR_test test_rotatePoint mpi_RR_test mpi_area_test test_mangle)

foreach (test1 ${testlist})
add_executable(${test1} ${test1}.cpp)
target_link_libraries(${test1} gtest bosslib)
endforeach(test1)

# Code below is for debugging
#add_executable(test_mangle_debug test_mangle.cpp ../mangle.cpp)
#set_target_properties(test_mangle_debug PROPERTIES COMPILE_FLAGS "-O0 -DDEBUGMODE")
