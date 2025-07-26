#include "lbm.h"

void handle_bcs_cavityflow(int N, gridpoint** grid, double lid_speed) {

    float rho_approx;

    // known velocity condition at top wall (lid):
    for(int x = 0; x < N; x++) {
        rho_approx = grid[x][N-1].f[0] + grid[x][N-1].f[1] + grid[x][N-1].f[3] + 2 * (grid[x][N-1].f[2] + grid[x][N-1].f[5] + grid[x][N-1].f[6]);
        grid[x][N-1].f[4] = grid[x][N-1].f[2];
        grid[x][N-1].f[7] = 0.5 * (grid[x][N-1].f[1] - grid[x][N-1].f[3] + 2 * grid[x][N-1].f[5] - rho_approx * lid_speed);
        grid[x][N-1].f[8] = -0.5 * (grid[x][N-1].f[1] - grid[x][N-1].f[3] - 2 * grid[x][N-1].f[6] - rho_approx * lid_speed); 
    }

    // zhou-he bounce-back on bottom wall:
    for(int x = 0; x < N; x++) {
        grid[x][0].f[2] = grid[x][0].f[4];
        grid[x][0].f[5] = grid[x][0].f[7];
        grid[x][0].f[6] = grid[x][0].f[8];
    }

    // zhou-he bounce-back on side walls excluding topmost node (velocity boundary):
    for(int y = 0; y < N-1; y++) {
        grid[0][y].f[1] = grid[0][y].f[3];
        grid[0][y].f[5] = grid[0][y].f[7];
        grid[0][y].f[8] = grid[0][y].f[6];

        grid[N-1][y].f[3] = grid[N-1][y].f[1];
        grid[N-1][y].f[7] = grid[N-1][y].f[5];
        grid[N-1][y].f[6] = grid[N-1][y].f[8];
    }
}

void handle_bcs_bounceback(int N, gridpoint** grid) {

    float rho_approx;

    // zhou-he bounce-back on top and bottom wall:
    for(int x = 0; x < N; x++) {
        grid[x][0].fstar[2] = grid[x][0].fstar[4];
        grid[x][0].fstar[5] = grid[x][0].fstar[7];
        grid[x][0].fstar[6] = grid[x][0].fstar[8];

        grid[x][N-1].fstar[4] = grid[x][N-1].fstar[2];
        grid[x][N-1].fstar[7] = grid[x][N-1].fstar[5];
        grid[x][N-1].fstar[8] = grid[x][N-1].fstar[6];
    }

    // zhou-he bounce-back on side walls excluding topmost node (velocity boundary):
    for(int y = 0; y < N-1; y++) {
        grid[0][y].fstar[1] = grid[0][y].fstar[3];
        grid[0][y].fstar[5] = grid[0][y].fstar[7];
        grid[0][y].fstar[8] = grid[0][y].fstar[6];

        grid[N-1][y].fstar[3] = grid[N-1][y].fstar[1];
        grid[N-1][y].fstar[7] = grid[N-1][y].fstar[5];
        grid[N-1][y].fstar[6] = grid[N-1][y].fstar[8];
    }
}

int main(int argc, double** argv) {
    int N = 40;
    double rho = 1.0;
    double reynolds_number = 200.0;
    int n_timesteps = 6;
    double vel = 0.1;
    double tau = compute_time_constant(vel, N, reynolds_number);    // 200 is reynolds number for this flow so we get ~laminar
    printf("Tau = %f\n", tau);

    gridpoint** grid = malloc(N * sizeof(gridpoint*));
    for(int i = 0; i < N; i++) {
        grid[i] = malloc(N * sizeof(gridpoint));
    }

    grid_initialize(N, N, rho, grid);

    for(int i = 0; i < 9; i++) {
        grid[N/2][N/2].f[i] = LBM_weight[i] * 2.0;
    }

    for(int t = 0; t < n_timesteps; t++) {
        grid_collision(N, N, tau, grid);
//        handle_bcs_cavityflow(N, grid, vel); // only one input N because we assume square domain for cavity flow
        handle_bcs_bounceback(N, grid);
       for(int x = 0; x < N; x++) for(int y = 0; y < N; y++) if (grid[x][y].density < 0) printf("t=%d\tNegative density at (%d,%d): rho = %f\n", t, x, y, grid[x][y].density);
        grid_stream(N, N, grid);
    }
    grid_draw(N, N, grid, 800, 800, 'd');
    grid_draw(N, N, grid, 800, 800, 'v');

    for(int i = 0; i < N; i++) {
        free(grid[i]);
    }
    free(grid);

    return 0;
}

/* NOTE TO SELF ON CURRENT ISSUES:

When running w/ current boundary conditions, the following happens:
    density is initialized properly, but the density at the velocity boundary approaches infinity while the remainder of the domain
    either (a) remains fixed at 1 or decreases (i cannot tell yet haven't debugged this)
        - this is caused by f_eq increasing far above 1 at the lid for each i = 1,...9 (only at the lid though, nowhere else)
    velocity remains 0 throughout the entire domain except at the boundary, where velocity is constant and properly applied

the question is then twofold: (1) why is density bunching up at the top? and (2) why is momentum not being transferred away from the lid?

*/