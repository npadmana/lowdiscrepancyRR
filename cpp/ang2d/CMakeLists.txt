add_executable(capintegrate capintegrate.cpp ang2d.cpp)
target_link_libraries(capintegrate npgsl m boost_program_options)

add_executable(deep2 deep2.cpp ang2d.cpp)
target_link_libraries(deep2 npgsl cfitsio m boost_program_options)

