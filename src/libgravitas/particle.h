//
// Created by nsamson on 3/26/16.
//

#ifndef GRAVITAS_PARTICLE_H
#define GRAVITAS_PARTICLE_H

#include <gsl/gsl_vector.h>
#include <inttypes.h>

#include "color.h"

typedef struct {
    gsl_vector* pos;
    gsl_vector* vel;
    gsl_vector* acc;
    color_t color;
} particle_t;

particle_t* particle_init(double rr, double phi, double init_speed);

#endif //GRAVITAS_PARTICLE_H
