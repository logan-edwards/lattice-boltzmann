#include "lbm.h"
#include <SDL2/SDL.h>

void grid_initialize(int Nx, int Ny, double rho, gridpoint** grid) {
    for(int x = 0; x < Nx; x++) {
        for(int y = 0; y < Ny; y++) {
            for(int i = 0; i < 9; i++) grid[x][y].f[i] = rho * LBM_weight[i];
        }
    }
}

void grid_collision(int Nx, int Ny, double tau, gridpoint** grid, gridpoint** swap_grid) {
    double lattice_density;
    vec2 lattice_momentum;
    double f_eq[9];
    double e_dot_u, u_dot_u;

    for(int x = 0; x < Nx; x++) {
        for(int y = 0; y < Ny; y++) {
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

                if(f_eq[i] > 1) printf("f_eq[%d] = %f\n at (%d,%d)", i, f_eq[i], x, y); // THIS IS THE PROBLEM POINT!!! top becomes infinitely dense

                /* Collision step */
                swap_grid[x][y].f[i] = grid[x][y].f[i] + (f_eq[i] - grid[x][y].f[i]) / tau;
            }
        }
        // TEST for velocity ~= 0.1:
        if(abs(grid[x][Ny-1].velocity.x - 0.1) > 1e-8) printf("Bad BC at x=%d, velocity is (%f,%f)", x, grid[x][Ny-1].velocity.x,grid[x][Ny-1].velocity.y);
    }
}

void grid_stream(int Nx, int Ny, gridpoint** grid, gridpoint** swap_grid) {
    for(int x = 1; x < Nx-1; x++) {
        for(int y = 1; y < Ny-1; y++) {
            for(int i = 0; i < 9; i++) {
                if(is_in_domain(Nx, Ny, x + LBM_e[i].x, y + LBM_e[i].y) == 1) {
                    grid[x + (int)LBM_e[i].x][y + (int)LBM_e[i].y].f[i] = swap_grid[x][y].f[i];
                }
            }
        }
    }
}

int is_in_domain(int Nx, int Ny, int x, int y) {
    if(x >= (Nx-1) || x <= 0) return 0;
    if(y >= (Ny-1) || y <= 0) return 0;
    
    return 1;
}

void grid_draw(int Nx, int Ny, gridpoint** grid, unsigned int screen_width, unsigned int screen_height, char mode) {
    
    int x, y; // x,y value on original grid
    // initialize color grid:
    color** color_grid = malloc(screen_width * sizeof(color*));
    for(int i = 0; i < screen_width; i++) {
        color_grid[i] = malloc(screen_height * sizeof(color));
    }
    if(mode == 'd') {
        double max_val, min_val;

        max_val = grid[0][0].density;
        min_val = max_val;
        for(int i = 0; i < Nx; i++) {
            for(int j = 0; j < Ny; j++) {
                if(grid[i][j].density > max_val) max_val = grid[i][j].density;
                if(grid[i][j].density < min_val) min_val = grid[i][j].density;
            }
        }
        printf("Minimum density = %f\n Maximum density = %f", min_val, max_val);
        for(int i = 0; i < screen_width; i++) {
            for(int j = 0; j < screen_height; j++) {
                x = (double)i * Nx / screen_width;
                y = -1.0 * j * Ny / screen_height + Ny;
                if(grid[x][y].density >= 0) {
                    color_grid[i][j].r = 0;
                    //color_grid[i][j].r = grid[x][y].velocity.x * 255.0 / (maxvel_x - minvel_x);
                    color_grid[i][j].g = 0;
                    color_grid[i][j].b = grid[x][y].density * 255 / (max_val - 0);  // interpolate such that 0 density is black. if desired, replace 0 with min_val
                    //color_grid[i][j].b = vec2_magnitude(grid[x][y].velocity) * 255.0 / (maxmag-minmag);
                    //color_grid[i][j].b = grid[x][y].velocity.y * 255.0 / (maxvel_y - minvel_y);
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
}
        else if (mode == 'v') {
            double maxvel_x, maxvel_y, minvel_x, minvel_y;
            maxvel_x = grid[0][0].velocity.x;
            maxvel_y = grid[0][0].velocity.y;
            minvel_x = maxvel_x;
            minvel_y = maxvel_y;
            for(int i = 0; i < Nx; i++) {
                for(int j = 0; j < Ny; j++) {
                    if(grid[i][j].velocity.x > maxvel_x) maxvel_x = grid[i][j].velocity.x;
                    if(grid[i][j].velocity.x < minvel_x) minvel_x = grid[i][j].velocity.x;
                    if(grid[i][j].velocity.y > maxvel_y) maxvel_y = grid[i][j].velocity.y;
                    if(grid[i][j].velocity.y < minvel_y) minvel_y = grid[i][j].velocity.y;
                }
            }
            printf("Max velocity components: %f in x, %f in y\n", maxvel_x, maxvel_y);
            printf("Min veloicty components: %f in x, %f in y\n", minvel_x, minvel_y);

            double maxmag, minmag;
            maxmag = vec2_magnitude(grid[0][0].velocity);
            for(int i = 0; i < Nx; i++) {
                for(int j = 0; j < Ny; j++) {
                    if(vec2_magnitude(grid[i][j].velocity) > maxmag) maxmag = vec2_magnitude(grid[i][j].velocity);
                    if(vec2_magnitude(grid[i][j].velocity) < minmag) minmag = vec2_magnitude(grid[i][j].velocity);
                }
            }
            printf("Max magnitude: %f\nMin magnitude: %f", maxmag, minmag);

            for(int i = 0; i < screen_width; i++) {
                for(int j = 0; j < screen_height; j++) {
                    x = (double)i * Nx / screen_width;
                    y = -1.0 * j * Ny / screen_height + Ny;
                    //if(grid[x][y].density >= 0) {
                        color_grid[i][j].r = 0;
                        //color_grid[i][j].r = grid[x][y].velocity.x * 255.0 / (maxvel_x - minvel_x);
                        color_grid[i][j].g = 0;
                        //color_grid[i][j].b = grid[x][y].density * 255 / (maxmag - minmag);  // interpolate such that 0 density is black. if desired, replace 0 with min_val
                        color_grid[i][j].b = vec2_magnitude(grid[x][y].velocity) * 255.0 / (maxmag-minmag);
                        //color_grid[i][j].b = grid[x][y].velocity.y * 255.0 / (maxvel_y - minvel_y);
                        color_grid[i][j].alpha = 255;
                    //}
                    /*else {
                        color_grid[i][j].r = 255;
                        color_grid[i][j].g = 0;
                        color_grid[i][j].b = 0;
                        color_grid[i][j].alpha = 255;                                   // sets negative densities to red to detect errors
                    }*/
                }
            }
}
        else if (mode == 'p') {
            printf("Pressure condition not yet implemented\n");
}
            else printf("Unknown mode \"%d\"\n", mode);
    
    

    // this is just for troubleshooting:
    
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

double vec2_magnitude(vec2 u) {
    return(sqrt(u.x * u.x + u.y * u.y));
}