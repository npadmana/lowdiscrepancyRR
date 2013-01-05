
# Build profiling codes
add_executable(mangle_pixelized mangle_pixelized.cpp)
target_link_libraries(mangle_pixelized bosslib)

add_executable(mangle_compare mangle_compare.cpp ../deprecated/mangle.cpp)
target_link_libraries(mangle_compare bosslib)
