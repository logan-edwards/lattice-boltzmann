#include "lbm.h"

void handle_bcs_cavityflow(int N, gridpoint** grid, double lid_speed) {
    float rho_approx;

    // zhou-he bounce-back on bottom wall:
    for(int x = 0; x < N; x++) {
        grid[x][0].f[2] = grid[x][0].f[4];
        grid[x][0].f[5] = grid[x][0].f[7];
        grid[x][0].f[6] = grid[x][0].f[8];
    }
    // zhou-he bounce-back on side walls:
    for(int y = 0; y < N; y++) {
        grid[0][y].f[1] = grid[0][y].f[3];
        grid[0][y].f[5] = grid[0][y].f[7];
        grid[0][y].f[8] = grid[0][y].f[6];

        grid[N-1][y].f[3] = grid[N-1][y].f[1];
        grid[N-1][y].f[7] = grid[N-1][y].f[5];
        grid[N-1][y].f[6] = grid[N-1][y].f[8];
    }

    // known velocity condition at top wall (lid):
    for(int x = 0; x < N; x++) {
        rho_approx = (1.0 / (1.0-lid_speed)) * (grid[x][N-1].f[0] + grid[x][N-1].f[1] + grid[x][N-1].f[3] + 2*(grid[x][N-1].f[2] + grid[x][N-1].f[5] + grid[x][N-1].f[6]));
        grid[x][N-1].f[4] = grid[x][N-1].f[2];
        //grid[x][N-1].f[7] = grid[x][N-1].f[5] + 0.5 * (grid[x][N-1].f[3] - grid[x][N-1].f[1]) - (1.0/6.0)*rho_approx*lid_speed;
        //grid[x][N-1].f[8] = grid[x][N-1].f[6] + 0.5 * (grid[x][N-1].f[1] - grid[x][N-1].f[3]) - (1.0/6.0)*rho_approx*lid_speed;
        grid[x][N-1].f[7] = grid[x][N-1].f[5] - 2 * rho_approx * lid_speed * LBM_weight[5] / (LBM_cs * LBM_cs);
        grid[x][N-1].f[8] = grid[x][N-1].f[6] + 2 * rho_approx * lid_speed * LBM_weight[6] / (LBM_cs * LBM_cs);

        /* SOMETHING IS FUNKY WITH THESE BCs. Neither of these are working (commented set is from AI, working set is from a paper lbm3) */
    }
}

int main(int argc, double** argv) {
    int N = 40;
    double rho = 1;
    int n_timesteps = 500;
    double vel = 0.1;
    double tau = compute_time_constant(vel, N, 200.0);    // 200 is reynolds number for this flow so we get ~laminar
    printf("Tau = %f\n", tau);

    gridpoint** grid = malloc(N * sizeof(gridpoint*));
    gridpoint** grid_copy = malloc(N * sizeof(gridpoint*));
    for(int i = 0; i < N; i++) {
        grid[i] = malloc(N * sizeof(gridpoint));
        grid_copy[i] = malloc(N * sizeof(gridpoint));
    }
    gridpoint** swap;

    grid_initialize(N, rho, grid);

    for(int t = 0; t < n_timesteps; t++) {
        grid_collision(N, tau, grid, grid_copy);
        handle_bcs_cavityflow(N, grid_copy, vel);
        grid_stream(N, grid, grid_copy);

        swap = grid;
        grid = grid_copy;
        grid_copy = swap;
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