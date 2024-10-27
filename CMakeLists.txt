cmake_minimum_required(VERSION 3.11)

project("Inertial_algorithm")

set(SRCS
    mathematics.cpp
    matrix.cpp
    )

set(SOURCE_LIB
    ../include/std.h
    ../include/matrix.h
    ../include/mathematics.h
    )

add_library(lclib STATIC ${SRCS})
target_include_directories(lclib PUBLIC "include")

add_executable(main main.cpp)
target_link_libraries(main lclib)
