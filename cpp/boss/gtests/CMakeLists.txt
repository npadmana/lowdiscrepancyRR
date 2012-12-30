set (testlist partition_test)

foreach (test1 ${testlist})
add_executable(${test1} ${test1}.cpp)
target_link_libraries(${test1} gtest gtest_main bosslib)
endforeach(test1)

