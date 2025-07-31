#include "lbm.h"

/*
Demo 2: Cavity Flow
Example of laminar lid-driven cavity flow in a square cavity

This validation requires a constant velocity boundary condition as well as
bounceback conditions on cavity walls. 
*/

void handle_bcs_cavity(int N, double lid_speed, gridpoint** grid, double tau) {
	double rho_approx;
	for(int x = 0; x < N; x++) {
		/* Bounceback on bottom wall */
        grid[x][0].f[2] = grid[x][0].f[4];
        grid[x][0].f[5] = grid[x][0].f[7];
        grid[x][0].f[6] = grid[x][0].f[8];
    }

    for(int y = 0; y < N-1; y++) {
		/* Bounceback on left wall */
        grid[0][y].f[1] = grid[0][y].f[3];
        grid[0][y].f[5] = grid[0][y].f[7];
        grid[0][y].f[8] = grid[0][y].f[6];

		/* Bounceback on right wall */
        grid[N-1][y].f[3] = grid[N-1][y].f[1];
        grid[N-1][y].f[7] = grid[N-1][y].f[5];
        grid[N-1][y].f[6] = grid[N-1][y].f[8];
    }

	/* Velocity condition on lid (top wall) */
	for(int x = 1; x < N-1; x++) {
		rho_approx = grid[x][N-1].f[0] + grid[x][N-1].f[1] + grid[x][N-1].f[3] 
			+ 2 * (grid[x][N-1].f[2] + grid[x][N-1].f[5] + grid[x][N-1].f[6]);

		grid[x][N-1].f[4] = grid[x][N-1].f[2];
		grid[x][N-1].f[7] = grid[x][N-1].f[5] 
			- 2.0 * rho_approx * lid_speed * LBM_weight[7] / (LBM_cs * LBM_cs);
		grid[x][N-1].f[8] = grid[x][N-1].f[6] 
			+ 2.0 * rho_approx * lid_speed * LBM_weight[8] / (LBM_cs * LBM_cs);
	}
}

int main() {
	/* Physical setup */
    int N = 80;
	int n_timesteps = 200000;
    double rho = 1.0;
    double reynolds_number = 200.0;
    double vel = 0.1;
    double tau = compute_time_constant(vel, 2*N, reynolds_number);
    printf("Tau = %f\n", tau);

	/* Memory allocation for grid */
    gridpoint** grid = malloc(2*N * sizeof(gridpoint*));
    for(int i = 0; i < 2*N; i++) {
        grid[i] = malloc(N * sizeof(gridpoint));
	}

	/* Initialize NxN grid with constant density */
	grid_initialize(N, N, rho, grid);

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
        compute_density_field(N, N, grid);
        compute_velocity_field(N, N, grid);
        compute_equilibrium_field(N, N, grid);

        grid_collision(N, N, tau, grid);
        handle_bcs_cavity(N, vel, grid, tau);
        grid_stream(N, N, grid);

        plot_update(N, N, grid, event, renderer, window, 800, 800, 'v');

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
	plot_grid(N, N, grid, 800, 800, 'v');
	plot_grid(N, N, grid, 800, 800, 'd');

    return 0;
}