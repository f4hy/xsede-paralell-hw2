
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <vector>

#define TIMERS 1


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
        std::vector<particle_t> TOP_LEFT;
        std::vector<particle_t> TOP_RIGHT;

        std::vector<particle_t> BOT_LEFT;
        std::vector<particle_t> BOT_RIGHT;

        double size =  sqrt(0.0005 * n);
        
        for(int i = 0; i < n; i++) {
            double x = particles[i].x;
            double y = particles[i].y;
            if(x<size/2){        // left
                if(y<size/2) {
                    BOT_LEFT.push_back(particles[i]);
                }
                else{
                    TOP_LEFT.push_back(particles[i]);
                }
            }
            else{
                if(y<size/2) {
                    BOT_RIGHT.push_back(particles[i]);
                }
                else{
                    TOP_RIGHT.push_back(particles[i]);
                }
            }
        }
        
        for(std::vector<particle_t>::iterator i = TOP_LEFT.begin(); i<TOP_LEFT.end(); i++){
            i->ax = i->ay = 0;
            for(std::vector<particle_t>::iterator j= TOP_LEFT.begin(); j<TOP_LEFT.end(); j++){
                apply_force(*i,*j, &dmin, &davg, &navg);
            }            
        }
        for(std::vector<particle_t>::iterator i= TOP_RIGHT.begin(); i<TOP_RIGHT.end(); i++){
            i->ax = i->ay = 0;
            for(std::vector<particle_t>::iterator j= TOP_RIGHT.begin(); j<TOP_RIGHT.end(); j++){
                apply_force(*i,*j, &dmin, &davg, &navg);
            }            
        }
        for(std::vector<particle_t>::iterator i= BOT_LEFT.begin(); i<BOT_LEFT.end(); i++){
            i->ax = i->ay = 0;
            for(std::vector<particle_t>::iterator j= BOT_LEFT.begin(); j<BOT_LEFT.end(); j++){  
                apply_force(*i,*j, &dmin, &davg, &navg);
            }            
        }
        for(std::vector<particle_t>::iterator i= BOT_RIGHT.begin(); i<BOT_RIGHT.end(); i++){ 
            i->ax = i->ay = 0;
            for(std::vector<particle_t>::iterator j= BOT_RIGHT.begin(); j<BOT_RIGHT.end(); j++){
                apply_force(*i,*j, &dmin, &davg, &navg);
            }            
        }

        // for(int i = 0; i < n; i++) {
        //     particles[i].ax = particles[i].ay = 0;
        //     for(int j = 0; j < n; j++) {
        //         apply_force(particles[i], particles[j], &dmin, &davg, &navg);
        //     }
        // }
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
