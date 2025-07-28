#include "lbm.h"

void handle_bcs_cavity(int N, double lid_speed, gridpoint** grid) {
	double rho_approx;
	// zhou-he bounceback
	for(int x = 0; x < N; x++) {
        grid[x][0].f[2] = grid[x][0].f[4];
        grid[x][0].f[5] = grid[x][0].f[7];
        grid[x][0].f[6] = grid[x][0].f[8];
    
    	//grid[x][N-1].f[4] = grid[x][N-1].f[2];
        //grid[x][N-1].f[7] = grid[x][N-1].f[5];
        //grid[x][N-1].f[8] = grid[x][N-1].f[6];
    }

    for(int y = 0; y < N-1; y++) {
        grid[0][y].f[1] = grid[0][y].f[3];
        grid[0][y].f[5] = grid[0][y].f[7];
        grid[0][y].f[8] = grid[0][y].f[6];

        grid[N-1][y].f[3] = grid[N-1][y].f[1];
        grid[N-1][y].f[7] = grid[N-1][y].f[5];
        grid[N-1][y].f[6] = grid[N-1][y].f[8];
    }

	// velocity condition
	for(int x = 0; x < N; x++) {
		rho_approx = grid[x][N-1].f[N-1] + grid[x][N-1].f[1] + grid[x][N-1].f[3] + 2 * (grid[x][N-1].f[4] + grid[x][N-1].f[7] + grid[x][N-1].f[8]);
		grid[x][N-1].f[4] = grid[x][N-1].f[2];
		grid[x][N-1].f[7] = grid[x][N-1].f[5] + 0.5 * (grid[x][N-1].f[1] - grid[x][N-1].f[3]) + (1.0/6.0) * rho_approx * lid_speed;
		grid[x][N-1].f[8] = grid[x][N-1].f[6] + 0.5 * (grid[x][N-1].f[3] - grid[x][N-1].f[1]) - (1.0/6.0) * rho_approx * lid_speed;
	}
}

int main() {
    int N = 80;
    double rho = 1.0;
    double reynolds_number = 200.0;
    int n_timesteps = 4000;
    double vel = 0.1;
    double tau = compute_time_constant(vel, N, reynolds_number);
    printf("Tau = %f\n", tau);

    gridpoint** grid = malloc(N * sizeof(gridpoint*));
    for(int i = 0; i < N; i++) {
        grid[i] = malloc(N * sizeof(gridpoint));
    }

    grid_initialize(N, N, rho, grid);

    SDL_Event event;
    SDL_Renderer *renderer;
    SDL_Window *window;
    SDL_Init(SDL_INIT_VIDEO);
    SDL_CreateWindowAndRenderer(800, 800, 0, &window, &renderer);
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
    SDL_RenderClear(renderer);

    for(int t = 0; t < n_timesteps; t++) {
        printf("Timestep %d... ", t);
        compute_density_field(N, N, grid);
        compute_velocity_field(N, N, grid);
        compute_equilibrium_field(N, N, grid);

        grid_collision(N, N, tau, grid);
        handle_bcs_cavity(N, vel, grid);
        grid_stream(N, N, grid);

        printf("done.\n");
        grid_plot_realtime(N, N, grid, event, renderer, window, 800, 800, 'v');
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

    for(int i = 0; i < N; i++) {
        free(grid[i]);
    }
    free(grid);

    return 0;
}