//
// Created by nsamson on 3/26/16.
//

#ifndef GRAVITAS_WORLD_H
#define GRAVITAS_WORLD_H

#include "particle.h"
#include <gsl/gsl_rng.h>

typedef struct {
    size_t x_size;
    size_t y_size;
    gsl_vector* pos_x;
    gsl_vector* pos_y;
    gsl_vector* vel_x;
    gsl_vector* vel_y;
    gsl_vector* acc_x;
    gsl_vector* acc_y;
    double* temperatures;
    size_t size;
    double mass;
} world_t;

typedef struct {
    double x;
    double y;
} coordinate_t;

world_t* world_init(gsl_rng* rand, size_t num_particles, size_t gen_radius, double init_speed, double mass);

world_t* world_step(world_t* world);

void print_world(const world_t* world);

color_t* pixel_buffer(const world_t* world);

coordinate_t center_of_mass(const world_t* world);

#endif //GRAVITAS_WORLD_H
