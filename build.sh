#!/bin/bash

g++ -g -I /home/nikita_viksne/Документы/C_Cpp_progs/InertialNavigation/include -I /home/nikita_viksne/Документы/C_Cpp_progs/InertialNavigation/ -c ./matrix.cpp -O2 -lm -o matrix.o

g++ -g -I /home/nikita_viksne/Документы/C_Cpp_progs/InertialNavigation/include -I /home/nikita_viksne/Документы/C_Cpp_progs/InertialNavigation/ -c ./main.cpp -O2 -lm -o main.o

g++ -g main.o matrix.o -O2 -o a.out 
