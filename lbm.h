#include <SDL2/SDL.h>
#include <math.h>

/* TO-DO:
1. Add convergence checks, i.e. store a grid and compare values every 100 steps or so, or calculate max difference between two grids and 
2. Add save and load feature for grids for two reasons: 1. so data can be saved automatically and reloaded (i.e. very 500 iterations) in event of crash, 
	and 2. for future visualizations so data, i.e. velocity, can be read from a structured csv file
3. CUDA speedups (this will take a good deal of work)

plan is to add these before continuing with demos (i.e. channel flow) but I do need to do another set of demos, i.e. to capture an outflow/pressure diff bc
	this can be accomplished with cavity flow first, then ideally (finally) a basic aerofoil or flow around a cylinder or something similar. or perhaps both
	time willing
*/

#ifndef LBM_HEADER
#define LBM_HEADER

typedef struct vec2 {
	double x;
	double y;
} vec2;

typedef struct gridpoint {
	vec2 velocity;
	double density;
	double pressure;
	double f[9];
	double fstar[9];
	double feq[9];
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

static const double LBM_cs = 1.0/1.73205080757; // lattice speed of sound = 1/sqrt(3)

/*
To implement:
grid_duplicate()
grid_save()
max_delta()

also clean up:
plot_grid()
plot_update()

*/

// grid manipulations
void grid_initialize(int Nx, int Ny, double rho, gridpoint** grid);
void grid_collision(int Nx, int Ny, double tau, gridpoint** grid);
void grid_stream(int Nx, int Ny, gridpoint** grid);
void grid_duplicate(int Nx, int Ny, gridpoint** grid_src, gridpoint** grid_target);
//void grid_save(int Nx, int Ny, gridpoint** grid, char* filename);
//void grid_load(gridpoint** grid, char* filename);

// helper functions
double compute_time_constant(double lattice_velocity, int lattice_length, double reynolds_number);
void compute_pressure_field(int Nx, int Ny, gridpoint** grid);
void compute_density_field(int Nx, int Ny, gridpoint** grid);
void compute_velocity_field(int Nx, int Ny, gridpoint** grid);
void compute_equilibrium_field(int Nx, int Ny, gridpoint** grid);
double vec2_magnitude(vec2 u);
int is_in_domain(int Nx, int Ny, int x, int y);

// plotting functions
void plot_grid(int Nx, int Ny, gridpoint** grid, unsigned int screen_width, unsigned int screen_height, char mode);
void plot_update(int Nx, int Ny, gridpoint** grid, SDL_Event event, 
	SDL_Renderer *renderer, SDL_Window *window, int screen_width, 
	unsigned int screen_height, char mode);

#endif