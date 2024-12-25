find_package (Threads)
set (libraries ${libraries}  ${CMAKE_THREAD_LIBS_INIT} hera)

add_executable(bottleneck_dist bottleneck_dist.cpp)
target_link_libraries(bottleneck_dist PUBLIC ${libraries})

if(HERA_BUILD_TESTS)
    add_subdirectory(tests)
endif()
