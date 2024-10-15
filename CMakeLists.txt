cmake_minimum_required(VERSION 2.8)

project(Inertial_algorithm)

set(SOURCE_EXE main.cpp ../include/std.h ../include/matrix.h)

set(SOURCE_LIB matrix.cpp ../include/std.h ../include/matrix.h)

add_library(matrix STATIC ${SOURCE_LIB})

add_executable(main ${SOURCE_EXE})

target_link_libraries(main foo)
