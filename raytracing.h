#ifndef __RAYTRACING_H
#define __RAYTRACING_H

#include "objects.h"
#include <stdint.h>

typedef struct {
    uint8_t *pixels;
    double* background_color;
    rectangular_node rectangulars;
    sphere_node spheres;
    light_node lights;
    const viewpoint *view;
    int width;
    int height;
    int thread_index;
    int thread_num;
} raypara_struct;

void raytracing(uint8_t *pixels, color background_color,
                rectangular_node rectangulars, sphere_node spheres,
                light_node lights, const viewpoint *view,
                int width, int height);
void ray_setpara(raypara_struct *raypara,
                 uint8_t *pixels, double* background_color,
                 rectangular_node rectangulars, sphere_node spheres,
                 light_node lights, const viewpoint *view,
                 int width, int height, int thread_index, int thread_num);
void ray_thread(void *ptr_raypara);
#endif
