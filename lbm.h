#include <SDL2/SDL.h>

#ifndef LBM_HEADER
#define LBM_HEADER

typedef struct vec2 {
	double x;
	double y;
} vec2;

typedef struct gridpoint {
	vec2 velocity;
	double density;
	double f[9];
} gridpoint;

typedef struct color {
	unsigned int r;
	unsigned int g;
	unsigned int b;
	unsigned int alpha;
} color;

static const double LBM_weight[9] = {4.0/9.0, 
    1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 
    1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0
};

static const vec2 LBM_e[9] = {
    {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0},
    {-1.0, 0.0}, {0.0, -1.0}, {1.0, 1.0},
    {-1.0, 1.0}, {-1.0, -1.0}, {1.0, -1.0}
};

static const double LBM_cs = 1/1.73205080757; // lattice speed of sound = 1/sqrt(3)

void grid_initialize(int N, double rho, gridpoint** grid);
void grid_collision(int N, double tau, gridpoint** grid, gridpoint** swap_grid);
void grid_stream(int N, gridpoint** grid, gridpoint** swap_grid);
void grid_draw(int N, gridpoint** grid, unsigned int screen_width, unsigned int screen_height);

int is_in_domain(int N, int x, int y);
#endif