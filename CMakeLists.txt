cmake_minimum_required(VERSION 3.4)
project(gravitas)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)
include(cmake_modules/FindGC.cmake)

find_package(GSL REQUIRED)
find_package(GC REQUIRED)
find_package(BLAS REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)


include_directories(${BOEHM_GC_INCLUDE_DIR} src/)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu11 -O3 -fopenmp")

set(SOURCE_FILES main.c)
set(LIB_FILES src/libgravitas/particle.h src/libgravitas/particle.c src/libgravitas/color.h src/libgravitas/world.h src/libgravitas/world.c src/libgravitas/color.c)

add_library(gravitas-core SHARED ${LIB_FILES})

add_executable(gravitas ${SOURCE_FILES})
target_link_libraries(gravitas-core gc gsl blas)
target_link_libraries(gravitas gravitas-core GL GLU glut)