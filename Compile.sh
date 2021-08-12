#! /bin/bash
g++ FitSpectra.C -o FitSpectra -lgsl -lgslcblas `root-config --libs --cflags`
echo "Compiled"
