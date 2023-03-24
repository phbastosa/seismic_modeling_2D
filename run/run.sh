#!/bin/bash

io="../src/utils/io/io.cpp"

eikonal_model="../src/model/eikonal/eikonal_model.cpp"
scalar_model="../src/model/scalar/scalar_model.cpp"

gaussian_2nd="../src/wavelet/gaussian_2nd/gaussian_2nd.cpp"

eikonal_modeling="../src/modeling/eikonal/eikonal_modeling.cpp"
scalar_modeling="../src/modeling/scalar/scalar_modeling.cpp"

regular_geom="../src/geometry/regular/regular.cpp"
streamer_geom="../src/geometry/streamer/streamer.cpp"

main="../src/modeling_main.cpp"

flags="-fopenmp -lm -O3"

g++ $flags $io $eikonal_model $scalar_model $gaussian_2nd $regular_geom $streamer_geom $eikonal_modeling $scalar_modeling $main -o ../bin/modeling.exe

./../bin/modeling.exe parameters.txt
