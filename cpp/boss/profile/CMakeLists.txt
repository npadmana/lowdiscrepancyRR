
# Build profiling codes
set (profilelist mangle_pixelized)

foreach (prof1 ${profilelist})
add_executable(${prof1} ${prof1}.cpp)
target_link_libraries(${prof1} bosslib)
endforeach(prof1)
