add_executable              (dionysus-core-tests core-tests.cpp)
target_link_libraries       (dionysus-core-tests PRIVATE dionysus)

add_test                    (NAME dionysus-core-tests COMMAND dionysus-core-tests)
