//
// Created by nsamson on 3/26/16.
//
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_double.h>
#include <string.h>
#include "particle.h"
#include "gc.h"
#include "math.h"
#include "world.h"
#include <omp.h>

#define MIN_DIST 5.0
#define MASS 1.00
#define BOUND 50.0
#define GCONST 1.0


double constrain_to_bound(double value, double bound);

static void world_finalizer(GC_PTR obj, GC_PTR __) {
    world_t* world = (world_t*) obj;

    gsl_vector_free(world->acc_x);
    gsl_vector_free(world->acc_y);
    gsl_vector_free(world->pos_x);
    gsl_vector_free(world->pos_y);
    gsl_vector_free(world->vel_x);
    gsl_vector_free(world->vel_y);
}

static world_t* empty_world(size_t size, size_t num_particles) {
    world_t* out = GC_MALLOC(sizeof(world_t));
    GC_register_finalizer(out, world_finalizer, NULL, NULL, NULL);

    out->size = num_particles;
    out->x_size = size;
    out->y_size = size;

    out->acc_x = gsl_vector_calloc(num_particles);
    out->acc_y = gsl_vector_calloc(num_particles);
    out->vel_x = gsl_vector_alloc(num_particles);
    out->vel_y = gsl_vector_alloc(num_particles);
    out->pos_x = gsl_vector_alloc(num_particles);
    out->pos_y = gsl_vector_alloc(num_particles);
    out->temperatures = GC_MALLOC_ATOMIC(num_particles * sizeof(double));
    memset(out->temperatures, 0, num_particles * sizeof(double));

    return out;

}

static void copy_world_data(world_t* dest, const world_t* src) {
    gsl_vector_memcpy(dest->vel_x, src->vel_x);
    gsl_vector_memcpy(dest->pos_y, src->pos_y);
    gsl_vector_memcpy(dest->vel_y, src->vel_y);
    gsl_vector_memcpy(dest->pos_x, src->pos_x);
    gsl_vector_memcpy(dest->acc_y, src->acc_y);
    gsl_vector_memcpy(dest->acc_x, src->acc_x);
    dest->mass = src->mass;
}

world_t *world_init(gsl_rng* rand_gen, size_t num_particles, size_t gen_radius, double init_speed, double mass) {
    world_t* out = GC_MALLOC(sizeof(world_t));
    GC_register_finalizer(out, world_finalizer, NULL, NULL, NULL);
    out->mass = mass;
    out->size = num_particles;
    out->x_size = gen_radius;
    out->y_size = gen_radius;

    out->acc_x = gsl_vector_calloc(num_particles);
    out->acc_y = gsl_vector_calloc(num_particles);
    out->vel_x = gsl_vector_alloc(num_particles);
    out->vel_y = gsl_vector_alloc(num_particles);
    out->pos_x = gsl_vector_alloc(num_particles);
    out->pos_y = gsl_vector_alloc(num_particles);
    out->temperatures = GC_MALLOC_ATOMIC(num_particles * sizeof(double));
    memset(out->temperatures, 0, num_particles * sizeof(double));
    gsl_vector* rr = gsl_vector_alloc(num_particles);
    gsl_vector* phi = gsl_vector_alloc(num_particles);
    gsl_vector* cos_phi = gsl_vector_alloc(num_particles);
    gsl_vector* sin_phi = gsl_vector_alloc(num_particles);
    gsl_vector* cos_2_phi = gsl_vector_alloc(num_particles);
    gsl_vector* sin_2_phi = gsl_vector_alloc(num_particles);

    gsl_vector_set_all(phi, 2*M_PI);

    for (size_t i = 0; i < num_particles; i++) {
        gsl_vector_set(rr, i, (gsl_rng_uniform(rand_gen) * gen_radius));
        gsl_vector_set(phi, i, gsl_vector_get(phi, i) * gsl_rng_uniform(rand_gen));
        gsl_vector_set(cos_phi, i, cos(gsl_vector_get(phi, i)));
        gsl_vector_set(sin_phi, i, sin(gsl_vector_get(phi, i)));
    }

    gsl_vector_add_constant(phi, M_PI_2);

    for (size_t i = 0; i < num_particles; i++) {
        gsl_vector_set(cos_2_phi, i, cos(gsl_vector_get(phi, i)));
        gsl_vector_set(sin_2_phi, i, sin(gsl_vector_get(phi, i)));
    }

    gsl_vector_memcpy(out->pos_x, rr);
    gsl_vector_mul(out->pos_x, cos_phi);

    gsl_vector_memcpy(out->pos_y, rr);
    gsl_vector_mul(out->pos_y, sin_phi);

    gsl_vector_memcpy(out->vel_x, rr);
    gsl_vector_mul(out->vel_x, cos_2_phi);
    gsl_vector_scale(out->vel_x, init_speed);

    gsl_vector_memcpy(out->vel_y, rr);
    gsl_vector_mul(out->vel_y, sin_2_phi);
    gsl_vector_scale(out->vel_y, init_speed);

    gsl_vector_free(rr);
    gsl_vector_free(phi);
    gsl_vector_free(cos_phi);
    gsl_vector_free(cos_2_phi);
    gsl_vector_free(sin_2_phi);
    gsl_vector_free(sin_phi);

    return out;
}

world_t* world_step(world_t *world) {

    world_t* new_data = empty_world(world->x_size, world->size);
    copy_world_data(new_data, world);

    #pragma omp parallel for
    for (size_t i = 0; i < world->size; i++) {
        for (size_t j = 0; j < i; j++) {
            double dx = gsl_vector_get(world->pos_x, i) - gsl_vector_get(world->pos_x, j);
            double dy = gsl_vector_get(world->pos_y, i) - gsl_vector_get(world->pos_y, j);

            dx = constrain_to_bound(dx, BOUND);
            dy = constrain_to_bound(dy, BOUND);

            double dist = sqrt(dx*dx + dy*dy);
            dist = constrain_to_bound(dist, BOUND);
            double inverse_dist = 1.0/dist;

            new_data->temperatures[i] = new_data->temperatures[i] + dist * 0.000001;

            double xd = dx * inverse_dist;
            double yd = dy * inverse_dist;
            double force_r = GCONST * world->mass * inverse_dist * inverse_dist;

            gsl_vector_set(new_data->vel_x, i, gsl_vector_get(new_data->vel_x, i) - force_r * xd);
            gsl_vector_set(new_data->vel_y, i, gsl_vector_get(new_data->vel_y, i) - force_r * yd);
        }

        for (size_t j = i + 1; j < world->size; j++) {
            double dx = gsl_vector_get(world->pos_x, i) - gsl_vector_get(world->pos_x, j);
            double dy = gsl_vector_get(world->pos_y, i) - gsl_vector_get(world->pos_y, j);
            dx = constrain_to_bound(dx, BOUND);
            dy = constrain_to_bound(dy, BOUND);

            double dist = sqrt(dx*dx + dy*dy);
            dist = constrain_to_bound(dist, BOUND);
            double inverse_dist = 1.0/dist;

            new_data->temperatures[i] = new_data->temperatures[i] + dist * 0.000001;

            double xd = dx * inverse_dist;
            double yd = dy * inverse_dist;
            double force_r = GCONST * world->mass * inverse_dist * inverse_dist;

            gsl_vector_set(new_data->vel_x, i, gsl_vector_get(new_data->vel_x, i) - force_r * xd);
            gsl_vector_set(new_data->vel_y, i, gsl_vector_get(new_data->vel_y, i) - force_r * yd);
        }
    }

    gsl_vector_add(new_data->pos_x, new_data->vel_x);
    gsl_vector_add(new_data->pos_y, new_data->vel_y);

    return new_data;
}

void print_world(const world_t* world) {

    for (size_t i = 0; i < world->size; i++) {
        printf("%" PRId64 ": %f, %f\n", i, gsl_vector_get(world->pos_x, i), gsl_vector_get(world->pos_y, i));
    }

}

coordinate_t center_of_mass(const world_t* world) {

    coordinate_t out;

    const double* x_data = gsl_vector_const_ptr(world->pos_x, 0);
    const double* y_data = gsl_vector_const_ptr(world->pos_y, 0);

    out.x = gsl_stats_mean(x_data, 1, world->size);
    out.y = gsl_stats_mean(y_data, 1, world->size);

    return out;
}


double inline constrain_to_bound(double value, double bound) {
    _Bool neg = value < 0.0;

    if (neg && value > -bound) {
        value = -bound;
    } else if (!neg && value < bound) {
        value = bound;
    }

    return value;
}