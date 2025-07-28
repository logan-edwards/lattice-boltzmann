#include "lbm.h"

void handle_bcs_bounceback(int N, gridpoint** grid) {
    for(int x = 0; x < N; x++) {
        grid[x][0].f[2] = grid[x][0].f[4];
        grid[x][0].f[5] = grid[x][0].f[7];
        grid[x][0].f[6] = grid[x][0].f[8];
    
        grid[x][N-1].f[4] = grid[x][N-1].f[2];
        grid[x][N-1].f[7] = grid[x][N-1].f[5];
        grid[x][N-1].f[8] = grid[x][N-1].f[6];
    }

    for(int y = 0; y < N; y++) {
        grid[0][y].f[1] = grid[0][y].f[3];
        grid[0][y].f[5] = grid[0][y].f[7];
        grid[0][y].f[8] = grid[0][y].f[6];

        grid[N-1][y].f[3] = grid[N-1][y].f[1];
        grid[N-1][y].f[7] = grid[N-1][y].f[5];
        grid[N-1][y].f[6] = grid[N-1][y].f[8];
    }
}

int main(int argc, double** argv) {
    int N = 64;
    double rho = 1.0;
    double reynolds_number = 200.0;
    int n_timesteps = 50;
    double vel = 0.1;
    double tau = compute_time_constant(vel, N, reynolds_number);    // 200 is reynolds number for this flow so we get ~laminar
    printf("Tau = %f\n", tau);

    gridpoint** grid = malloc(N * sizeof(gridpoint*));
    for(int i = 0; i < N; i++) {
        grid[i] = malloc(N * sizeof(gridpoint));
    }

    grid_initialize(N, N, rho, grid);

    for(int i = 0; i < 9; i++) grid[N/2][N/2].f[i] = LBM_weight[i] * 5.0;

    for(int t = 0; t < n_timesteps; t++) {
        printf("Timestep %d... ", t);
        compute_density_field(N, N, grid);
        compute_velocity_field(N, N, grid);
        compute_equilibrium_field(N, N, grid);


        grid_collision(N, N, tau, grid);
        handle_bcs_bounceback(N, grid);
        grid_stream(N, N, grid);
        printf("done.\n");
    }

    grid_draw(N, N, grid, 800, 800, 'd');
    grid_draw(N, N, grid, 800, 800, 'v');

    for(int i = 0; i < N; i++) {
        free(grid[i]);
    }
    free(grid);

    return 0;
}