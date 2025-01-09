cmake_minimum_required(VERSION 3.11)

project("Inertial_algorithm")

set(SRCS
    mathematics.cpp
    matrix.cpp
    )

set(SOURCE_LIB
    include/std.h
    include/matrix.h
    include/mathematics.h
    )

add_library(lclib STATIC ${SRCS})
target_include_directories(lclib PUBLIC "include")

# Tell CMake to run moc when necessary:
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
# As moc files are generated in the binary dir, tell CMake
# to always look for includes there:
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Widgets finds its own dependencies.
find_package(Qt5Widgets REQUIRED)

add_executable(main main.cpp)

qt5_use_modules(main Core)

target_link_libraries(main lclib)
