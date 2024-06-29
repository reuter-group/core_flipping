cmake_minimum_required(VERSION 3.0)
project(wham)

file(GLOB_RECURSE sources *.cu *.h)

find_package(CUDA 9.0)

# add_definitions(-DNOMPI)
cuda_add_executable(wham ${sources})
# target_compile_options(md PUBLIC -DNOMPI -g -O3 -arch=sm_20)
