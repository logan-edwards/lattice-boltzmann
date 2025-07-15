//#include "lbm.h"
#include <SDL2/SDL.h>
#include "lbm.h"

/*
pseudocode:

- setup grid w/ initial conditions
- for each cell:
	- apply bcs if applicable
	- then, compute density
	- then, compute velocity
	- reload visualizer
*/

const float weight[9] = {4.0/9.0, 
    1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 
    1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0
};

void copy_grid(unsigned int nx, unsigned int ny, gridpoint** target, gridpoint** source) {
    for(int i = 0; i < nx; i++) {
        for(int j = 0; j < ny; j++) {
            target[i][j] = source[i][j];
        }
    }
}

void initialize_grid(unsigned int nx, unsigned int ny, double xmax, double ymax, gridpoint** grid) {
    double dx, dy;
    dx = xmax / (nx + 1);
    dy = ymax / (ny + 1);
    for(int i = 0; i < nx; i++) {
        for(int j = 0; j < ny; j++) {
            grid[i][j].coordinate.x = i*dx;
            grid[i][j].coordinate.y = j*dy;
            for(int k = 0; k < 9; k++) {
                grid[i][j].weight[k] = weight[k];
            }
            grid[i][j].velocity.x = 0;
            grid[i][j].velocity.y = 0;
            grid[i][j].density = 1.0;
            grid[i][j].boundary_condition = 0;
        }
    }
}

void step_grid(unsigned int nx, unsigned int ny, gridpoint** grid1, gridpoint** grid2) {
    double density;
    vec2 velocity;
    for(int i = 0; i < nx; i++) {
        for(int j = 0; j < ny; j++) {
            for(int k = 0; k < )
        }
    }
}

void main() {
    unsigned int nx, ny;
    float xmax, float ymax;
    nx = 5;
    ny = 5;

    gridpoint** grid1 = malloc(nx * sizeof(gridpoint*));
    gridpoint** grid2 = malloc(nx * sizeof(gridpoint*));
    for(int i = 0; i < nx; i++) {
        grid1[i] = malloc(ny * sizeof(gridpoint));
        grid2[i] = malloc(ny * sizeof(gridpoint));
    }
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