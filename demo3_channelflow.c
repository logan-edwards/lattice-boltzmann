#include "lbm.h"

/*
Demo 3: Simple free channel flow

This demo involves a velocity inlet on the left wall, a pressure outlet on the
right wall, and bounceback conditions on the top and bottom walls. Length is
chosen such that L = 2*H
*/

void handle_bcs_channel(int Nx, int Ny, double inlet_speed, gridpoint** grid) {
    double rho_approx;
	for(int x = 0; x < Nx; x++) {
        /* Bounceback on top wall */
        grid[x][Ny-1].f[4] = grid[x][Ny-1].f[2];
        grid[x][Ny-1].f[7] = grid[x][Ny-1].f[5];
        grid[x][Ny-1].f[8] = grid[x][Ny-1].f[6];

        /* Bounceback on bottom wall */
        grid[x][0].f[2] = grid[x][0].f[4];
        grid[x][0].f[5] = grid[x][0].f[7];
        grid[x][0].f[6] = grid[x][0].f[8];
    }

    for(int y = 1; y < Ny-1; y++) {
        rho_approx = grid[0][y].f[0] + grid[0][y].f[2] + grid[0][y].f[4] + 2 * (grid[0][y].f[6] + grid[0][y].f[7] + grid[0][y].f[3]);
        /* Inlet on left wall */
        grid[0][y].f[1] = grid[Nx-1][y].f[3] + 2 * rho_approx * LBM_weight[1] * inlet_speed / (LBM_cs * LBM_cs);
        grid[0][y].f[5] = grid[Nx-1][y].f[7] - 2 * rho_approx * LBM_weight[5] * inlet_speed / (LBM_cs * LBM_cs);
        grid[0][y].f[8] = grid[Nx-1][y].f[6] - 2 * rho_approx * LBM_weight[8] * inlet_speed / (LBM_cs * LBM_cs);

        /* Outlet on right wall */
        grid[Nx-1][y].f[3] = grid[Nx-2][y].f[3];
        grid[Nx-1][y].f[6] = grid[Nx-2][y].f[6];
        grid[Nx-1][y].f[7] = grid[Nx-2][y].f[7];
    }
}

int main() {
	/* Physical setup */
    int N = 40;
	int n_timesteps = 500;
    double rho = 1.0;
    double reynolds_number = 200.0;
    double vel = 0.2;
    double tau = compute_time_constant(vel, N, reynolds_number);
    printf("Tau = %f\n", tau);

	/* Memory allocation for grid */
    gridpoint** grid = malloc(2*N * sizeof(gridpoint*));
    for(int i = 0; i < 2*N; i++) {
        grid[i] = malloc(N * sizeof(gridpoint));
	}

	/* Initialize NxN grid with constant density */
	grid_initialize(2*N, N, rho, grid);

	/* Set up renderer to display each timestep */
    SDL_Event event;
    SDL_Renderer *renderer;
    SDL_Window *window;
    SDL_Init(SDL_INIT_VIDEO);
    SDL_CreateWindowAndRenderer(2*800, 800, 0, &window, &renderer);
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
    SDL_RenderClear(renderer);

	/* The main loop: collision step --> streaming step */
	for(int t = 0; t < n_timesteps; t++) {
		printf("Iteration %d...", t);
        compute_density_field(2*N, N, grid);
        compute_velocity_field(2*N, N, grid);
        compute_equilibrium_field(2*N, N, grid);

        grid_collision(2*N, N, tau, grid);
        handle_bcs_channel(2*N, N, vel, grid);
        grid_stream(2*N, N, grid);

        plot_update(2*N, N, grid, event, renderer, window, 2*800, 800, 'v');

		printf("done.\n");
    }

	SDL_RenderPresent(renderer);
    while (1) {
        if (SDL_PollEvent(&event) && event.type == SDL_QUIT) {
            break;
        }
    }
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

	/* Plot velocity and density after simulation ends */
	plot_grid(N, N, grid, 2*800, 800, 'v');
	plot_grid(N, N, grid, 2*800, 800, 'd');

    return 0;
}