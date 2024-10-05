#!/bin/bash

g++ -g -I /home/nikita_viksne/Документы/C_Cpp_progs/InertialNavigation/include -I /home/nikita_viksne/Документы/C_Cpp_progs/InertialNavigation/ -c ./matrix.cpp -lm -o matrix.o

g++ -g -I /home/nikita_viksne/Документы/C_Cpp_progs/InertialNavigation/include -I /home/nikita_viksne/Документы/C_Cpp_progs/InertialNavigation/ -c ./main.cpp -lm -o main.o

g++ -g main.o matrix.o -o a.out 
