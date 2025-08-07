cmake_minimum_required(VERSION 3.11)

project("Inertial_algorithm" LANGUAGES CXX)

set(SRCS
    mathematics.cpp
    matrix.cpp
	quaternions.cpp
    )

set(SOURCE_LIB
    ../include/std.h
    ../include/matrix.h
    ../include/mathematics.h
    /home/nikita_viksne/Документы/CLibraries/alglib-cpp-v4.04.0/src/
    )


add_executable(main main.cpp)

target_link_libraries(main lclib "/home/nikita_viksne/InertialNavigation/libkalman.a" "/home/nikita_viksne/InertialNavigation/libalglib.a")

