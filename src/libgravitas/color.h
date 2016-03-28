//
// Created by nsamson on 3/26/16.
//

#ifndef GRAVITAS_COLOR_H
#define GRAVITAS_COLOR_H

#include <stdint.h>

typedef struct {
    uint8_t red;
    uint8_t green;
    uint8_t blue;
    uint8_t alpha;
} color_t;

color_t convert_temperature(double temp);

#endif //GRAVITAS_COLOR_H
