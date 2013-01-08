
# Build profiling codes
add_executable(mangle_pixelized mangle_pixelized.cpp)
target_link_libraries(mangle_pixelized bosslib)

add_executable(mangle_compare mangle_compare.cpp ../deprecated/mangle.cpp)
target_link_libraries(mangle_compare bosslib)

add_executable(mangle_compare2 mangle_compare2.cpp)
target_link_libraries(mangle_compare2 bosslib)


add_executable(rng_speed rng_speed.cpp)
target_link_libraries(rng_speed bosslib)