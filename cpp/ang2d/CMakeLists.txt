add_executable(capintegrate capintegrate.cpp ang2d.cpp)
target_link_libraries(capintegrate npgsl m boost_program_options)

add_executable(deep2 deep2.cpp ang2d.cpp)
target_link_libraries(deep2 npgsl cfitsio m boost_program_options)

add_executable(test_mangle test_mangle.cpp mangle.cpp ang2d.cpp)
add_executable(test_bossmask test_bossmask.cpp BossMask.cpp mangle.cpp ang2d.cpp)
add_executable(boss boss.cpp mangle.cpp BossMask.cpp ang2d.cpp)
target_link_libraries(test_bossmask npgsl m)
target_link_libraries(boss npgsl m boost_program_options)
configure_file(../../data/boss/boss_survey.ply ${CMAKE_CURRENT_BINARY_DIR}/boss_survey.ply COPYONLY)
configure_file(../../data/boss/boss_geometry_2011_06_10_LOWZ-north.ply 
         ${CMAKE_CURRENT_BINARY_DIR}/boss_geometry_2011_06_10_LOWZ-north.ply COPYONLY)