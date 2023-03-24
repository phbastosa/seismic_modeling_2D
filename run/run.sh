#!/bin/bash

io="../src/utils/io/io.cpp"

eikonal_model="../src/model/eikonal/eikonal_model.cpp"
scalar_model="../src/model/scalar/scalar_model.cpp"
acoustic_model="../src/model/acoustic/acoustic_model.cpp"

gaussian_1st="../src/wavelet/gaussian_1st/gaussian_1st.cpp"
gaussian_2nd="../src/wavelet/gaussian_2nd/gaussian_2nd.cpp"

eikonal_modeling="../src/modeling/eikonal/eikonal_modeling.cpp"
scalar_modeling="../src/modeling/scalar/scalar_modeling.cpp"
acoustic_modeling="../src/modeling/acoustic/acoustic_modeling.cpp"

regular_geom="../src/geometry/regular/regular.cpp"
streamer_geom="../src/geometry/streamer/streamer.cpp"

main="../src/modeling_main.cpp"

flags="-fopenmp -lm -O3"

g++ $flags $io $eikonal_model $scalar_model $acoustic_model $gaussian_1st $gaussian_2nd $regular_geom $streamer_geom $eikonal_modeling $scalar_modeling $acoustic_modeling $main -o ../bin/modeling.exe

./../bin/modeling.exe parameters.txt
