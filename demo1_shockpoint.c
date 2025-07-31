#include "lbm.h"

/*
Demo 1: Shock Point
Example of shock induced by initial point of relative high pressure

This validation is chosen for two reasons:
    (1) setup is simple, involving only no-slip (bounceback) wall boundary
        conditions and a single initial condition on density.
    (2) quick verification by evaluating symmetry. Solution is expected to be
        symmetric ring as shockwave travels through domain.
*/

void handle_bcs_bounceback(int N, gridpoint** grid) {
    for(int x = 0; x < N; x++) {
        /* Bounceback condition on domain bottom */
        grid[x][0].f[2] = grid[x][0].f[4];
        grid[x][0].f[5] = grid[x][0].f[7];
        grid[x][0].f[6] = grid[x][0].f[8];
    
        /* Bounceback condition on domain top */
        grid[x][N-1].f[4] = grid[x][N-1].f[2];
        grid[x][N-1].f[7] = grid[x][N-1].f[5];
        grid[x][N-1].f[8] = grid[x][N-1].f[6];
    }


    for(int y = 0; y < N; y++) {
        /* Bounceback condition on left wall */
        grid[0][y].f[1] = grid[0][y].f[3];
        grid[0][y].f[5] = grid[0][y].f[7];
        grid[0][y].f[8] = grid[0][y].f[6];

        /* Bounceback condition on right wall */
        grid[N-1][y].f[3] = grid[N-1][y].f[1];
        grid[N-1][y].f[7] = grid[N-1][y].f[5];
        grid[N-1][y].f[6] = grid[N-1][y].f[8];
    }
}

int main(int argc, double** argv) {
    /* Physical setup */
    int n_timesteps = 100;
    int N = 200;
    double rho = 1.0;
    double reynolds_number = 2000.0;
    double vel = 0.1;
    double tau = compute_time_constant(vel, N, reynolds_number);
    printf("Tau = %f\n", tau);

    /* Memory allocation for grid */
    gridpoint** grid = malloc(N * sizeof(gridpoint*));
    for(int i = 0; i < N; i++) {
        grid[i] = malloc(N * sizeof(gridpoint));
    }

    /* Initialize NxN grid with constant density */
    grid_initialize(N, N, rho, grid);

    /* Set initial condition as 5x density at center point */
    for(int i = 0; i < 9; i++) grid[N/2][N/2].f[i] = LBM_weight[i] * 5.0;

    /* Set up renderer to display each timestep */
    SDL_Event event;
    SDL_Renderer *renderer;
    SDL_Window *window;
    SDL_Init(SDL_INIT_VIDEO);
    SDL_CreateWindowAndRenderer(800, 800, 0, &window, &renderer);
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
    SDL_RenderClear(renderer);

    /* The main loop: collision step --> streaming step */
    for(int t = 0; t < n_timesteps; t++) {
        printf("Timestep %d... ", t);
        compute_density_field(N, N, grid);
        compute_velocity_field(N, N, grid);
        compute_equilibrium_field(N, N, grid);

        grid_collision(N, N, tau, grid);
        handle_bcs_bounceback(N, grid);
        grid_stream(N, N, grid);

        printf("done.\n");
        plot_update(N, N, grid, event, renderer, window, 800, 800, 'v');
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

    /* Plot velocity and density after simulation ends */
    plot_grid(N, N, grid, 800, 800, 'v');
	plot_grid(N, N, grid, 800, 800, 'd');

    return 0;
}