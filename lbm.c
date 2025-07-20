#include "lbm.h"
#include <stdio.h>

/*
The general operational steps:
    (0) Convert the values (i.e. initial conditions) into lattice units
    (I) start with a set of f_i (initial) values in each cell
    (II) algorithm:
        1. compute rho = sum f_i, u = (1/rho) * sum(v_i f_i)    where v_i is a VECTOR
        2. determine feq_i using velocity and density
        3. perform collision step:      f*_i = f_i - (1/tau) * [f_i - feq_i]
        4. handle boundary conditions
        5. perform streaming step:      f*_i = f_i      i.e. f_i <--- f*_i since computations are finished
    (III) [i imagine]: convert the normalized (lattice unit) values into SI units
        - some trickery with the Reynold's number is needed here
*/

const double weight[9] = {4.0/9.0, 
    1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 
    1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0
};
const double e_x[9] = {0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0};
const double e_y[9] = {0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0};

const vec2 direction[9] = {
    {0.0,0.0}, {1.0, 0.0}, {0.0, 1.0},
    {-1.0, 0.0}, {0.0, -1.0}, {1.0, 1.0},
    {-1.0, 1.0}, {-1.0, -1.0}, {1.0, -1.0}
};

const double c_s = 1/1.73205080757; // lattice speed of sound = 1/sqrt(3)

double tau; // note: this is ideally calculated in terms of the initial (physical, SI) values, i.e. viscosity
double dx;
double dy;
unsigned int Nx;
unsigned int Ny;

double dotprod2(vec2 u, vec2 v) {
    double sum;
    sum = u.x + v.x + u.y + v.y;
    return sum;
}

void grid_initialize(gridpoint** grid) {
    // perchance do step 0 up here as well
    for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            grid[i][j].coordinate.x = i * dx;
            grid[i][j].coordinate.y = j * dy;
            for(int k = 0; k < 9; k++) {
                grid[i][j].f[k] = weight[k]; // In reality, this is mathematically f_i = w_i*rho_i for an initial rho. For unit density, f_i = w_i
            }
        }
    }
}

void grid_step(gridpoint** grid, gridpoint** swap_grid) {
    float f_eq[9];
    for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            /* Computation of density at each node, using rho = sum_i f_i */
            grid[i][j].density = 0;
            for(int k = 0; k < 9; k++) {
                grid[i][j].density = grid[i][j].density + grid[i][j].f[k];
            }

            /* Computation of velocity at each node, using the vector sum
            u = (1/rho) * sum_i e_i f_i */
            grid[i][j].velocity.x = 0;
            grid[i][j].velocity.y = 0;
            for(int k = 0; k < 9; k++) {
                grid[i][j].velocity.x = grid[i][j].velocity.x + e_x[k]*grid[i][j].f[k];
                grid[i][j].velocity.y = grid[i][j].velocity.y + e_y[k]*grid[i][j].f[k];
            }
            grid[i][j].velocity.x /= grid[i][j].density;
            grid[i][j].velocity.y /= grid[i][j].density;

            /* computation of f_eq using taylor expansion (the long equation) */
            for(int k = 0; k < 9; k++) {
                f_eq[k] = weight[k] * grid[i][j].density * (
                    1 + 
                    dotprod2(direction[k], grid[i][j].velocity)/(c_s*c_s)
                    + dotprod2(direction[k], grid[i][j].velocity)*dotprod2(direction[k], grid[i][j].velocity)/(2*c_s*c_s*c_s*c_s)
                    - dotprod2(grid[i][j].velocity, grid[i][j].velocity)/(2*c_s*c_s)
                );
            }

            /* Collision step */
            for(int k = 0; k < 9; k++) {
                swap_grid[i][j].f[k] = grid[i][j].f[k] - (1/tau) * (grid[i][j].f[k] - f_eq[k]);
            }

            /* 
            
            HANDLE BOUNDARY CONDITIONS HERE
            
            */

            /* Streaming Step */
            for(int k = 0; k < 9; k++) {
                grid[i][j].f[k] = swap_grid[i][j].f[k];
            }
        }
    }
}

void grid_draw(gridpoint** grid, unsigned int screen_width, unsigned int screen_height, int var_type) {
    /* var_type index:
        1 = LB Density
        2 = LB Velocity
        3 = LB Pressure (using R-B visualizer, R if <0 and B if >0)
        4 = SI Density
        5 = SI Velocity
        6 = SI Pressure (using R-B visualizer, R if <0 and B if >0)
    */
    double max_val, min_val;
    int x, y; // x,y value on original grid

    // right now, we just do this with density:
    max_val = grid[0][0].density;
    min_val = max_val;
    for(int i = 0; i < Nx; i++) {
        for(int j = 0; j < Ny; j++) {
            if(grid[i][j].density > max_val) max_val = grid[i][j].density;
            if(grid[i][j].density < min_val) min_val = grid[i][j].density;
        }
    }

    // initialize color grid:
    color** color_grid = malloc(screen_width * sizeof(color*));
    for(int i = 0; i < screen_width; i++) {
        color_grid[i] = malloc(screen_height * sizeof(color));
    }
    for(int i = 0; i < screen_width; i++) {
        for(int j = 0; j < screen_height; j++) {
            x = (double)i * Nx / screen_width;
            y = -1.0 * j * Ny / screen_height + Ny;
            color_grid[i][j].r = 0;
            color_grid[i][j].g = 0;
            color_grid[i][j].b = grid[x][y].density * 255 / (max_val - min_val); // interpolation works, but grid[x][y].density is not returning proper values 
            color_grid[i][j].alpha = 255;
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

/* TO-DO: read a config file from main containing:
    - Time of Simulation
    - Viscosity
    - Length
    - Height
    - Top Wall BC
    - Right Wall BC
    - Left Wall BC
    - Bottom Wall BC

*/
int main(int argc, double** argv) {
    double viscocity;
    double time;
    double length;
    double height;
    // things to be read by config file:
    Nx = 40;
    Ny = 40;
    viscocity = 1;
    time = 100;
    length = 1;
    height = 1;

    gridpoint** grid = malloc(Nx * sizeof(gridpoint*));
    gridpoint** grid_copy = malloc(Nx * sizeof(gridpoint*));
    for(int i = 0; i < Nx; i++) {
        grid[i] = malloc(Ny * sizeof(gridpoint));
        grid_copy[i] = malloc(Ny * sizeof(gridpoint));
    }

    grid_initialize(grid);
    grid_step(grid, grid_copy);
    grid_draw(grid, 1000, 800, 1);

    for(int i = 0; i < Nx; i++) {
        free(grid[i]);
        free(grid_copy[i]);
    }
    free(grid);
    free(grid_copy);
}