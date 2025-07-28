#include <SDL2/SDL.h>
#include <math.h>

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

void grid_initialize(int Nx, int Ny, double rho, gridpoint** grid);
void grid_collision(int Nx, int Ny, double tau, gridpoint** grid);
void grid_stream(int Nx, int Ny, gridpoint** grid);
void grid_draw(int Nx, int Ny, gridpoint** grid, unsigned int screen_width, unsigned int screen_height, char mode);
void grid_plot_realtime(int Nx, int Ny, gridpoint** grid, SDL_Event event, SDL_Renderer *renderer, SDL_Window *window, int screen_width, unsigned int screen_height, char mode);

int is_in_domain(int Nx, int Ny, int x, int y);
double vec2_magnitude(vec2 u);

double compute_time_constant(double lattice_velocity, int lattice_length, double reynolds_number);
void compute_pressure_field(int Nx, int Ny, gridpoint** grid);
void compute_density_field(int Nx, int Ny, gridpoint** grid);
void compute_velocity_field(int Nx, int Ny, gridpoint** grid);
void compute_equilibrium_field(int Nx, int Ny, gridpoint** grid);


#endif