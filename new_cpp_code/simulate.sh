#!/bin/bash

# This script is used to simulate the data for the paper
rm *out

# compile the program
g++ -O3 para.cpp kd_tree.cpp particle.cpp test.cpp

./a.out >> data.out
