add_executable              (relative-zz-filtration     relative-zz-filtration.cpp)
target_link_libraries       (relative-zz-filtration     ${libraries})

add_executable              (relative-lzz     relative-lzz.cpp ../src/cnpy.cpp)
target_link_libraries       (relative-lzz     ${libraries})
