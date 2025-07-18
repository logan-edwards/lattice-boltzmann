#include "lbm.h"
#include <SDL2/SDL.h>
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

void grid_draw(gridpoint** grid, unsigned int variable_type, char color_style) {

    /* Outline of style variables:
        color_style:
            'r' = red-variable theme
            'g' = green-variable theme
            'b' = blue-variable theme (default for other variables)
            'h' = high-resolution theme (need to figure this out, will vary over several colors)
        variable_type:
            0 = LB-Density (default for variable_type > 5)
            1 = LB-Velocity
            2 = LB-Pressure
            3 = SI-Density
            4 = SI-Velocity
            5 = SI-Pressure
    */
    double max, min;
    unsigned int rgb_max, rgb_min;
    unsigned int screen_width, screen_height;

    // temp stuff for testing:
    screen_width = 800;
    screen_height = 800;

    SDL_Event event;
    SDL_Renderer *renderer;
    SDL_Window *window;

    SDL_Init(SDL_INIT_VIDEO);
    SDL_CreateWindowAndRenderer(screen_width, screen_height, 0, &window, &renderer);
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
    SDL_RenderClear(renderer);
    for(int i = 0; i < screen_width; i++) {
        for(int j = 0; j < screen_height; j++) {
            SDL_SetRenderDrawColor(renderer, 0, 0, i/4, 255);
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


void main() {
    Nx = 5;
    Ny = 5;

    gridpoint** grid = malloc(Nx * sizeof(gridpoint*));
    gridpoint** grid_copy = malloc(Nx * sizeof(gridpoint*));
    for(int i = 0; i < Nx; i++) {
        grid[i] = malloc(Ny * sizeof(gridpoint));
        grid_copy[i] = malloc(Ny * sizeof(gridpoint));
    }

    grid_initialize(grid);
    grid_step(grid, grid_copy);
    grid_draw(grid, 0, 'b');
}

/*

// SDL example main function:

int main(void) {
    SDL_Event event;
    SDL_Renderer *renderer;
    SDL_Window *window;
    int i;

    int windowwidth = 800;

    SDL_Init(SDL_INIT_VIDEO);
    SDL_CreateWindowAndRenderer(windowwidth, windowwidth, 0, &window, &renderer);
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0);
    SDL_RenderClear(renderer);
    SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
    for (i = 0; i < windowwidth; ++i)
        SDL_RenderDrawPoint(renderer, i, i);
    SDL_RenderPresent(renderer);
    while (1) {
        if (SDL_PollEvent(&event) && event.type == SDL_QUIT)
            break;
    }
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return EXIT_SUCCESS;
}
*/