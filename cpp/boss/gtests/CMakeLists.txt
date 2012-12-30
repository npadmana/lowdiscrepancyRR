
# Test data files
set (testdata hemisphere.cfg hemisphere.ply hemisphere_input.cfg sphere.cfg sphere.ply sphere_input.cfg)
foreach(file1 ${testdata})
  configure_file(${file1} ${CMAKE_CURRENT_BINARY_DIR}/${file1} COPYONLY)
endforeach(file1)


# Build tests
set (testlist partition_test area_test)

foreach (test1 ${testlist})
add_executable(${test1} ${test1}.cpp)
target_link_libraries(${test1} gtest gtest_main bosslib)
endforeach(test1)

