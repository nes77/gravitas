//
// Created by nsamson on 3/26/16.
//


#include <gsl/gsl_rng.h>
#include "particle.h"
#include "gc.h"
#include "math.h"
#include "world.h"

static void particle_finalizer(GC_PTR obj, GC_PTR client_data) {
    particle_t* ptr = (particle_t*) obj;
    gsl_vector_free(ptr->acc);
    gsl_vector_free(ptr->pos);
    gsl_vector_free(ptr->vel);
}

particle_t *particle_init(double rr, double phi, double init_speed) {
    particle_t* out = GC_MALLOC_ATOMIC(sizeof(particle_t));
    GC_register_finalizer(out, particle_finalizer, NULL, NULL, NULL);

    out->pos = gsl_vector_alloc(2);
    out->vel = gsl_vector_alloc(2);
    out->acc = gsl_vector_calloc(2);

    gsl_vector_set(out->pos, 0, rr * cos(phi));
    gsl_vector_set(out->pos, 1, rr * sin(phi));
    gsl_vector_set(out->vel, 0, rr * cos(phi + M_PI_2) * init_speed);
    gsl_vector_set(out->vel, 1, rr * sin(phi + M_PI_2) * init_speed);

    return out;
}





