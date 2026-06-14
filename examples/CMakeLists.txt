add_library             (dionysus_examples INTERFACE)
target_link_libraries   (dionysus_examples INTERFACE dionysus)
target_include_directories(dionysus_examples INTERFACE include)

add_subdirectory        (filtration)
add_subdirectory        (rips)
add_subdirectory        (lzz)
add_subdirectory        (omni-field)
