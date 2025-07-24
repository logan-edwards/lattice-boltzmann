#include "lbm.h"
#include <SDL2/SDL.h>

void grid_initialize(int N, double rho, gridpoint** grid) {
    for(int x = 0; x < N; x++) {
        for(int y = 0; y < N; y++) {
            for(int i = 0; i < 9; i++) grid[x][y].f[i] = rho * LBM_weight[i];
        }
    }
}

void grid_collision(int N, double tau, gridpoint** grid, gridpoint** swap_grid) {
    double lattice_density;
    vec2 lattice_momentum;
    double f_eq[9];
    double e_dot_u, u_dot_u;

    for(int x = 0; x < N; x++) {
        for(int y = 0; y < N; y++) {
            /* Computation of density and velocity */
            lattice_density = 0;
            lattice_momentum.x = 0;
            lattice_momentum.y = 0;
            for(int i = 0; i < 9; i++) lattice_density += grid[x][y].f[i];
            grid[x][y].density = lattice_density;
            for(int i = 0; i < 9; i++) {
                lattice_momentum.x += grid[x][y].f[i] * LBM_e[i].x;
                lattice_momentum.y += grid[x][y].f[i] * LBM_e[i].y;
            }
            grid[x][y].velocity.x = lattice_momentum.x / lattice_density;
            grid[x][y].velocity.y = lattice_momentum.y / lattice_density;

            /* Computation of f_eq */
            for(int i = 0; i < 9; i++) {
                e_dot_u = LBM_e[i].x * grid[x][y].velocity.x + LBM_e[i].y * grid[x][y].velocity.y;
                u_dot_u = grid[x][y].velocity.x * grid[x][y].velocity.x + grid[x][y].velocity.y * grid[x][y].velocity.y;
                
                f_eq[i] = LBM_weight[i]*grid[x][y].density * (
                1 +
                3 * e_dot_u / (LBM_cs * LBM_cs) 
                + 9 * e_dot_u * e_dot_u / (2 * LBM_cs * LBM_cs * LBM_cs * LBM_cs) 
                - 3 * u_dot_u / (2 * LBM_cs * LBM_cs)
                );

                /* Collision step */
                swap_grid[x][y].f[i] = grid[x][y].f[i] + (f_eq[i] - grid[x][y].f[i]) / tau;
            }
        }
    }
}

void grid_stream(int N, gridpoint** grid, gridpoint** swap_grid) {
    for(int x = 0; x < N; x++) {
        for(int y = 0; y < N; y++) {
            for(int i = 0; i < 9; i++) {
                if(is_in_domain(N, x+LBM_e[i].x, y+LBM_e[i].y) == 1) {
                    grid[x + (int)LBM_e[i].x][y + (int)LBM_e[i].y].f[i] = swap_grid[x][y].f[i];
                }
            }
        }
    }
}

int is_in_domain(int N, int x, int y) {
    if(x > (N-1) || x < 0) return 0;
    if(y > (N-1) || y < 0) return 0;
    
    return 1;
}

void grid_draw(int N, gridpoint** grid, unsigned int screen_width, unsigned int screen_height) {
    double max_val, min_val;
    int x, y; // x,y value on original grid

    max_val = grid[0][0].density;
    min_val = max_val;
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            if(grid[i][j].density > max_val) max_val = grid[i][j].density;
            if(grid[i][j].density < min_val) min_val = grid[i][j].density;
        }
    }
    printf("Minimum density = %f\n Maximum density = %f", min_val, max_val);
    // initialize color grid:
    color** color_grid = malloc(screen_width * sizeof(color*));
    for(int i = 0; i < screen_width; i++) {
        color_grid[i] = malloc(screen_height * sizeof(color));
    }
    for(int i = 0; i < screen_width; i++) {
        for(int j = 0; j < screen_height; j++) {
            x = (double)i * N / screen_width;
            y = -1.0 * j * N / screen_height + N;
            if(grid[x][y].density >= 0) {
                color_grid[i][j].r = 0;
                color_grid[i][j].g = 0;
                color_grid[i][j].b = grid[x][y].density * 255 / (max_val - 0);  // interpolate such that 0 density is black. if desired, replace 0 with min_val
                color_grid[i][j].alpha = 255;
            }
            else {
                color_grid[i][j].r = 255;
                color_grid[i][j].g = 0;
                color_grid[i][j].b = 0;
                color_grid[i][j].alpha = 255;                                   // sets negative densities to red to detect errors
            }
        }
    }
    SDL_Event event;
    SDL_Renderer *renderer;
    SDL_Window *window;
    SDL_Init(SDL_INIT_VIDEO);
    SDL_CreateWindowAndRenderer(screen_width, screen_height, 0, &window, &renderer);
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
    SDL_RenderClear(renderer);
    for(int i = 0; i < screen_width; i++) {
        for(int j = 0; j < screen_height; j++) {
            SDL_SetRenderDrawColor(renderer, 0, 0, color_grid[i][j].b, 0);
            SDL_RenderDrawPoint(renderer, i, j);
        }
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
}

double compute_time_constant(double lattice_velocity, int lattice_length, double reynolds_number) {
    double lattice_viscosity, tau;
    lattice_viscosity = lattice_velocity * lattice_length / reynolds_number;
    tau = 3.0 * lattice_viscosity + 0.5;
    return(tau);
}