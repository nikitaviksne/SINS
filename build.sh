#!/bin/bash

# matrix.cpp
g++ -g -I /home/nikita_viksne/Документы/C_Cpp_progs/InertialNavigation/include -I /home/nikita_viksne/Документы/C_Cpp_progs/InertialNavigation/ -c ./matrix.cpp -O7 -lm -o matrix.o

# mathematics.cpp
g++ -g -I /home/nikita_viksne/Документы/C_Cpp_progs/InertialNavigation/include -I /home/nikita_viksne/Документы/C_Cpp_progs/InertialNavigation/ -c ./mathematics.cpp -O7 -lm -o mathematics.o

# main.cpp
g++ -g -I /home/nikita_viksne/Документы/C_Cpp_progs/InertialNavigation/include -I /home/nikita_viksne/Документы/C_Cpp_progs/InertialNavigation/ -c ./main.cpp -O7 -lm -o main.o

g++ -g main.o matrix.o mathematics.o -O7 -o a.out 
