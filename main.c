#include <stdio.h>
#include <gc.h>
#include <stdlib.h>
#include "libgravitas/world.h"
#include <GL/glut.h>
#include <unistd.h>

world_t* world;

void init_world(double);
void display();

uint64_t counter = 0;

int main(int argc, char** argv) {
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);

    gsl_rng_set(rng, 50UL);

    size_t size = 200;

    world = world_init(rng, 250, size, 0.0, 100.0);

    gsl_rng_free(rng);

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize((int) 1000, (int) 1000);
    glutInitWindowPosition(000, 000);
    glutCreateWindow("Gravitas!");
    init_world((double) 1000);
    glutDisplayFunc(display);
    glutIdleFunc(display);

    glutMainLoop();

}

void init_world(double size) {
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glMatrixMode(GL_PROJECTION);
    gluOrtho2D(-size, size, -size, size);
    glPointSize(3.0f);
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT);

    glBegin(GL_POINTS);

    glPointSize(3.0f);
    for (size_t i = 0; i < world->size; i++) {

        glColor3b((GLbyte) (i % 4 == 3 ? 127 : 63), (GLbyte) (i % 4 == 0 ? 127 : 63), (GLbyte) (i % 2 == 0 ? 127 : 63));
        glVertex2d(
                gsl_vector_get(world->pos_x, i),
                gsl_vector_get(world->pos_y, i)
        );
    }

    coordinate_t com = center_of_mass(world);

    glColor3b(0, 0, 0);

    glPointSize(25.0f);
    glVertex2d(com.x, com.y);

    glEnd();

    world = world_step(world);
    usleep(50000);

    glFlush();
    glutSwapBuffers();
}