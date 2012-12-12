add_executable(capintegrate capintegrate.cpp)
target_link_libraries(capintegrate npgsl m boost_program_options)

add_executable(deep2 deep2.cpp)
target_link_libraries(deep2 npgsl cfitsio m boost_program_options)

add_executable(test_mangle test_mangle.cpp mangle.cpp)
configure_file(../../data/boss/boss_survey.ply ${CMAKE_CURRENT_BINARY_DIR}/boss_survey.ply COPYONLY)