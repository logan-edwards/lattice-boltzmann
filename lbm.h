#include <SDL2/SDL.h>
#include <stdio.h>

#ifndef LBM_HEADER
#define LBM_HEADER

typedef struct vec2 {
	double x;
	double y;
} vec2;

typedef struct gridpoint {
	vec2 coordinate;
	vec2 velocity;
	double density;
	double f[9];
} gridpoint;

typedef struct boundary {
	int condition;
	vec2 coordinate;
	vec2 velocity;
} boundary;

typedef struct color {
	unsigned int r;
	unsigned int g;
	unsigned int b;
	unsigned int alpha;
} color;

/*
TO-DO with BCs:
	Generalized Bounce-Back. For fixed wall, v=0 and this is a no-slip condition
	For wall with v = [vx, vy] this can function as a velocity inlet

	Zero-Gradient Boundary. This is an "open boundary", i.e. a pressure outlet.

*/

void grid_initialize(gridpoint** grid);
void grid_step(gridpoint** grid, gridpoint** swap_grid);
void grid_draw(gridpoint** grid, unsigned int screen_width, unsigned int screen_height, int var_type);
double dotprod2(vec2 u, vec2 v);

#endif