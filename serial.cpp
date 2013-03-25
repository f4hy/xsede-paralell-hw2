#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <vector>

#define TIMERS 0


//
//  benchmarking program
//
int main(int argc, char **argv)
{
    int navg, nabsavg = 0;
    double davg, dmin, absmin = 1.0, absavg = 0.0;

    if(find_option(argc, argv, "-h") >= 0) {
        printf("Options:\n");
        printf("-h to see this help\n");
        printf("-n <int> to set the number of particles\n");
        printf("-o <filename> to specify the output file name\n");
        printf("-s <filename> to specify a summary file name\n");
        printf("-no turns off all correctness checks and particle output\n");
        return 0;
    }

    int n = read_int(argc, argv, "-n", 1000);

    char *savename = read_string(argc, argv, "-o", NULL);
    char *sumname = read_string(argc, argv, "-s", NULL);

    FILE *fsave = savename ? fopen(savename, "w") : NULL;
    FILE *fsum = sumname ? fopen(sumname, "a") : NULL;

    particle_t *particles = (particle_t*) malloc(n * sizeof(particle_t));
    set_size(n);
    init_particles(n, particles);

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer();
#if TIMERS==1
    double force_time = 0.0;
    double move_time = 0.0;
    double running_time;
#endif

    for(int step = 0; step < NSTEPS; step++) {
        navg = 0;
        davg = 0.0;
        dmin = 1.0;
        //
        //  compute forces
        //
#if TIMERS==1
        running_time = read_timer();
#endif
        std::vector< std::vector<particle_t*> > blocks;
        std::vector< std::vector<particle_t*> > blocks_buffered;
        
        double size =  sqrt(0.0005 * n);
        double buffer = 0.01;

        int subdiv = 6;
        for(int sx = 0; sx<subdiv; sx++){
            for(int sy = 0; sy<subdiv; sy++){
                double left = sx*(size/subdiv);
                double right = (sx+1)*(size/subdiv);
                double bot = sy*(size/subdiv);
                double top = (sy+1)*(size/subdiv);
                std::vector<particle_t*> block;
                std::vector<particle_t*> block_buffered;
                
                for(int i = 0; i < n; i++) {
                    double x = particles[i].x;
                    double y = particles[i].y;

                    if(left<x && x<=right && bot<y && y<=top){
                        block.push_back(particles+i);
                    }
                    if(left-buffer<x && x<=right+buffer && bot-buffer<y && y<=top+buffer){
                        block_buffered.push_back(particles+i);
                    }
                }
                blocks.push_back(block);
                blocks_buffered.push_back(block_buffered);
            }
        }
        std::vector< std::vector<particle_t*> >::iterator b_b = blocks_buffered.begin();
        for(std::vector< std::vector<particle_t*> >::iterator block = blocks.begin(); block<blocks.end(); block++){
            for(std::vector<particle_t*>::iterator i = block->begin(); i<block->end(); i++){
                (*i)->ax = (*i)->ay = 0;
                for(std::vector<particle_t*>::iterator j = b_b->begin(); j<b_b->end(); j++){
                    apply_force(**i,**j, &dmin, &davg, &navg);
                }
            }
            b_b++;
        }
#if TIMERS==1
        force_time += read_timer() - running_time;
        running_time = read_timer();
#endif
//
//  move particles
//
        for(int i = 0; i < n; i++) {
            move(particles[i]);
        }
#if TIMERS==1
        move_time += read_timer() - running_time;
#endif

        if(find_option(argc, argv, "-no") == -1) {
            //
            // Computing statistical data
            //
            if(navg) {
                absavg +=  davg / navg;
                nabsavg++;
            }
            if(dmin < absmin) {
                absmin = dmin;
            }

            //
            //  save if necessary
            //
            if(fsave && (step % SAVEFREQ) == 0) {
                save(fsave, n, particles);
            }
        }
    }
    simulation_time = read_timer() - simulation_time;

    printf("n = %d, simulation time = %g seconds", n, simulation_time);

    if(find_option(argc, argv, "-no") == -1) {
        if(nabsavg) {
            absavg /= nabsavg;
        }
        //
        //  -the minimum distance absmin between 2 particles during the run of the simulation
        //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
        //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
        //
        //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
        //
        printf(", absmin = %lf, absavg = %lf", absmin, absavg);
        if(absmin < 0.4) {
            printf("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
        }
        if(absavg < 0.8) {
            printf("\nThe average distance is below 0.8 meaning that most particles are not interacting");
        }
    }
    printf("\n");
#if TIMERS==1
    printf(" force time = %g seconds", force_time);
    printf(" move time = %g seconds\n", move_time);
#endif

//
// Printing summary data
//
    if(fsum) {
        fprintf(fsum, "%d %g\n", n, simulation_time);
    }

//
// Clearing space
//
    if(fsum) {
        fclose(fsum);
    }
    free(particles);
    if(fsave) {
        fclose(fsave);
    }

    return 0;
}
