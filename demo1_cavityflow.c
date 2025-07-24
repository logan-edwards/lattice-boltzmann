#include "lbm.h"

int main(int argc, double** argv) {
    int N = 40;
    double rho = 1;
    int n_timesteps = 50;
    double tau = 0.55;

    gridpoint** grid = malloc(N * sizeof(gridpoint*));
    gridpoint** grid_copy = malloc(N * sizeof(gridpoint*));
    for(int i = 0; i < N; i++) {
        grid[i] = malloc(N * sizeof(gridpoint));
        grid_copy[i] = malloc(N * sizeof(gridpoint));
    }

    grid_initialize(N, rho, grid);
	grid_collision(N, tau, grid, grid_copy);
	grid_stream(N, grid, grid_copy);
    for(int t = 0; t < n_timesteps; t++) {
        grid_collision(N, tau, grid, grid_copy);
        grid_stream(N, grid, grid_copy);
    }
    grid_draw(N, grid, 800, 800);


    for(int i = 0; i < N; i++) {
        free(grid[i]);
        free(grid_copy[i]);
    }
    free(grid);
    free(grid_copy);

    return 0;
}