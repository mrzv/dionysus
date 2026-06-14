add_executable              (relative-zz-filtration     relative-zz-filtration.cpp)
target_link_libraries       (relative-zz-filtration     dionysus_examples)

add_executable              (relative-lzz     relative-lzz.cpp ../src/cnpy.cpp)
target_link_libraries       (relative-lzz     dionysus_examples)

add_executable              (extended-persistence     extended-persistence.cpp ../src/cnpy.cpp)
target_link_libraries       (extended-persistence     dionysus_examples)
