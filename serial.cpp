#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <vector>

#define TIMERS 0

#include <stdint.h>
static inline uint32_t log2(const uint32_t x) {
  uint32_t y;
  asm ( "\tbsr %1, %0\n"
        : "=r"(y)
        : "r" (x)
        );
  return y;
}
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

    const uint32_t n = read_int(argc, argv, "-n", 1000);

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

    double size =  sqrt(0.0005 * n);
    double interaction_length = 0.01;


    int blocksize = (int) ceil(size/interaction_length);
    double block_width = size/blocksize;

    printf("blocksize :%d\n", blocksize);

    particle_t*** blocks = (particle_t***) malloc(blocksize*blocksize * sizeof(particle_t**));
    for(int b=0; b<blocksize*blocksize; b++){
        blocks[b] = (particle_t**)malloc(n*sizeof(particle_t*));
    }

    
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
        // std::vector< std::vector<particle_t*> > blocks;
        // std::vector< std::vector<particle_t*> > blocks_buffered;

        
        // particle_t* blocks = (particle_t*) malloc(blocksize*blocksize*n * sizeof(particle_t));
        // particle_t* blocks = new particle_t[blocksize*blocksize][n];
            
        // particle_t* blocks[blocksize*blocksize][n];
        int number_in_block[blocksize*blocksize];

        // printf("initalize\n");
        for(int b=0; b<blocksize*blocksize; b++){
            number_in_block[b] = 0; // starts with no particles in any box;
        }
        // printf("setup\n");

        for(int p = 0; p < n; p++) {
            double x = particles[p].x;
            double y = particles[p].y;

            int x_index = (int)floor(x/block_width);
            int y_index = (int)floor(y/block_width);
            // assert (x_index >=0 && x_index <blocksize );
            // assert (y_index >=0 && y_index <blocksize );
            int particle_index = number_in_block[x_index + y_index*blocksize]++;
            // assert (particle_index < n);            
            blocks[x_index + y_index*blocksize][particle_index] = particles+p;
        }
        // printf("applying force\n");
        // for(int p=0; p<n; p++){
        //     printf("%d    : x %f, y %f, vx %f, vy %f, ax %f, ay %f \n",p,particles[p].x,particles[p].y,
        //            particles[p].vx,particles[p].vy,
        //            particles[p].ax,particles[p].ay);
        // }
        // int count = 0;
        // for(int i=0; i<blocksize; i++){
        //     for(int j=0; j<blocksize; j++){
        //         count += number_in_block[i + j*blocksize];
        //         printf("blocksize %d thisblock[%d][%d] %d \n",blocksize,i,j,number_in_block[i + j*blocksize]);
        //     }
        // }
        // assert(count == n);
        // if(count != n){
        //     printf("ERROR");
        //     exit(1);
        // }
        bool mustquit = false;
        for(int i=0; i<blocksize; i++){
            for(int j=0; j<blocksize; j++){
                for(int p=0; p<number_in_block[i + j*blocksize]; p++ ){
                    blocks[i + j*blocksize][p]->ax = blocks[i + j*blocksize][p]->ay = 0;
                    // printf("%d %d %d %d\n",step, i,j,p);
                    // printf("%d,%d,%d: x %f, y %f, vx %f, vy %f, ax %f, ay %f \n",i,j,p,(blocks[i + j*blocksize][p])->x,(blocks[i + j*blocksize][p])->y,
                    //        (blocks[i + j*blocksize][p])->vx,(blocks[i + j*blocksize][p])->vy,
                    //        (blocks[i + j*blocksize][p])->ax,(blocks[i + j*blocksize][p])->ay);
                    // (blocks[i + j*blocksize][p])->vx = 0;
                    // particle_t interacting_p = blocks[i + j*blocksize][p];

                    
                    // double preax = (blocks[i + j*blocksize][p])->ax;
                    // double preay = (blocks[i + j*blocksize][p])->ay;
                    
                    // for(int foo=0; foo<n; foo++){
                    //     apply_force(*(blocks[i + j*blocksize][p]), particles[foo], &dmin, &davg, &navg);

                    // }
                    // double finalax = (blocks[i + j*blocksize][p])->ax;
                    // double finalay = (blocks[i + j*blocksize][p])->ay;
                    // blocks[i + j*blocksize][p]->ax = blocks[i + j*blocksize][p]->ay = 0;
                    //interact with this block
                    for(int xoffset=-2; xoffset<3; xoffset++){
                        int xblockindex = i+xoffset;
                        if(xblockindex<0 || xblockindex >= blocksize ){
                            continue;
                        }
                        for(int yoffset=-2; yoffset<3; yoffset++){
                            int yblockindex = j+yoffset;
                            if(yblockindex<0 || yblockindex >= blocksize ){
                                continue;
                            }
                            // printf("%d %d  ",xblockindex, yblockindex);
                            for(int num=0; num<number_in_block[xblockindex + yblockindex*blocksize]; num++ ){
                                apply_force(*(blocks[i + j*blocksize][p]), *(blocks[xblockindex + yblockindex*blocksize][num]), &dmin, &davg, &navg);

                                }
                        }
                    }
                    // printf(" \n");
                    double postax = (blocks[i + j*blocksize][p])->ax;
                    double postay = (blocks[i + j*blocksize][p])->ay;

                    
                    // if(postax != finalax){
                    //     printf("%d %d %d %d: %f\n",step,i,j,p,postax);
                    //     printf("%d         : %f\n",step,finalax);
                    //     printf("%d %d,%d,%d: x %f, y %f, vx %f, vy %f, ax %f, ay %f \n",step,i,j,p,(blocks[i + j*blocksize][p])->x,(blocks[i + j*blocksize][p])->y,
                    //            (blocks[i + j*blocksize][p])->vx,(blocks[i + j*blocksize][p])->vy,
                    //            (blocks[i + j*blocksize][p])->ax,(blocks[i + j*blocksize][p])->ay);
                    //     mustquit = 1;
                    // }
                    // assert(postax == finalax);
                    // assert(postax == finalax);
                    // for(int n=0; n<number_in_block[i + j*blocksize]; n++ ){

                    //     apply_force(*(blocks[i + j*blocksize][p]), *(blocks[i + j*blocksize][n]), &dmin, &davg, &navg);
                    // }
                    // printf("%d,%d,%d: x %f, y %f, vx %f, vy %f, ax %f, ay %f \n",i,j,p,(blocks[i + j*blocksize][p])->x,(blocks[i + j*blocksize][p])->y,
                    //        (blocks[i + j*blocksize][p])->vx,(blocks[i + j*blocksize][p])->vy,
                    //        (blocks[i + j*blocksize][p])->ax,(blocks[i + j*blocksize][p])->ay);
                    // // interact with i-1 neighbors
                    // if(i!=0){
                    //     for(int n=0; n<number_in_block[i-1 + j*blocksize]; n++ ){
                    //         apply_force(*(blocks[i + j*blocksize][p]), *(blocks[i-1 + j*blocksize][n]), &dmin, &davg, &navg);
                    //     }
                    // }
                    // printf("%d,%d,%d: x %f, y %f, vx %f, vy %f, ax %f, ay %f \n",i,j,p,(blocks[i + j*blocksize][p])->x,(blocks[i + j*blocksize][p])->y,
                    //        (blocks[i + j*blocksize][p])->vx,(blocks[i + j*blocksize][p])->vy,
                    //        (blocks[i + j*blocksize][p])->ax,(blocks[i + j*blocksize][p])->ay);
                    // // interact with j-1 neighbors
                    // if(j!=0){
                    //     for(int n=0; n<number_in_block[i + (j-1)*blocksize]; n++ ){
                    //         apply_force(*(blocks[i + j*blocksize][p]), *(blocks[i + (j-1)*blocksize][n]), &dmin, &davg, &navg);
                    //     }
                    // }
                    // printf("%d,%d,%d: x %f, y %f, vx %f, vy %f, ax %f, ay %f \n",i,j,p,(blocks[i + j*blocksize][p])->x,(blocks[i + j*blocksize][p])->y,
                    //        (blocks[i + j*blocksize][p])->vx,(blocks[i + j*blocksize][p])->vy,
                    //        (blocks[i + j*blocksize][p])->ax,(blocks[i + j*blocksize][p])->ay);
                    // // interact with i+1 neighbors
                    // if(i!=blocksize-1){
                    //     for(int n=0; n<number_in_block[(i+1) + j*blocksize]; n++ ){
                    //         apply_force(*(blocks[i + j*blocksize][p]), *(blocks[(i+1) + j*blocksize][n]), &dmin, &davg, &navg);
                    //     }
                    // }
                    // printf("%d,%d,%d: x %f, y %f, vx %f, vy %f, ax %f, ay %f \n",i,j,p,(blocks[i + j*blocksize][p])->x,(blocks[i + j*blocksize][p])->y,
                    //        (blocks[i + j*blocksize][p])->vx,(blocks[i + j*blocksize][p])->vy,
                    //        (blocks[i + j*blocksize][p])->ax,(blocks[i + j*blocksize][p])->ay);
                    // // interact with j+1 neighbors
                    // if(j!=blocksize-1){
                    //     for(int n=0; n<number_in_block[i + (j+1)*blocksize]; n++ ){
                    //         apply_force(*(blocks[i + j*blocksize][p]), *(blocks[i + (j+1)*blocksize][n]), &dmin, &davg, &navg);
                    //     }
                    // }
                    // printf("%d,%d,%d: x %f, y %f, vx %f, vy %f, ax %f, ay %f \n",i,j,p,(blocks[i + j*blocksize][p])->x,(blocks[i + j*blocksize][p])->y,
                    //        (blocks[i + j*blocksize][p])->vx,(blocks[i + j*blocksize][p])->vy,
                    //        (blocks[i + j*blocksize][p])->ax,(blocks[i + j*blocksize][p])->ay);

                }
            }
        }
        if(mustquit){
            for(int i=0; i<blocksize; i++){
                for(int j=0; j<blocksize; j++){
                    for(int p=0; p<number_in_block[i + j*blocksize]; p++ ){
                        printf("%d,%d,%d: x %f, y %f, vx %f, vy %f, ax %f, ay %f \n",i,j,p,(blocks[i + j*blocksize][p])->x,(blocks[i + j*blocksize][p])->y,
                               (blocks[i + j*blocksize][p])->vx,(blocks[i + j*blocksize][p])->vy,
                               (blocks[i + j*blocksize][p])->ax,(blocks[i + j*blocksize][p])->ay);
                    }
                }
            }
            for(int p=0; p<n; p++){
                printf("%d    : x %f, y %f, vx %f, vy %f, ax %f, ay %f \n",p,particles[p].x,particles[p].y,
                       particles[p].vx,particles[p].vy,
                       particles[p].ax,particles[p].ay);
            }
        exit(1);
        }

        // int subdiv = max(4,log2(n));
        // for(int sx = 0; sx<subdiv; sx++){
        //     for(int sy = 0; sy<subdiv; sy++){
        //         double left = sx*(size/subdiv);
        //         double right = (sx+1)*(size/subdiv);
        //         double bot = sy*(size/subdiv);
        //         double top = (sy+1)*(size/subdiv);
        //         std::vector<particle_t*> block;
        //         std::vector<particle_t*> block_buffered;

        //         for(int i = 0; i < n; i++) {
        //             double x = particles[i].x;
        //             double y = particles[i].y;

        //             if(left<x && x<=right && bot<y && y<=top){
        //                 block.push_back(particles+i);
        //             }
        //             if(left-buffer<x && x<=right+buffer && bot-buffer<y && y<=top+buffer){
        //                 block_buffered.push_back(particles+i);
        //             }
        //         }
        //         blocks.push_back(block);
        //         blocks_buffered.push_back(block_buffered);
        //     }
        // }
        // std::vector< std::vector<particle_t*> >::iterator b_b = blocks_buffered.begin();
        // for(std::vector< std::vector<particle_t*> >::iterator block = blocks.begin(); block<blocks.end(); block++){
        //     for(std::vector<particle_t*>::iterator i = block->begin(); i<block->end(); i++){
        //         (*i)->ax = (*i)->ay = 0;
        //         for(std::vector<particle_t*>::iterator j = b_b->begin(); j<b_b->end(); j++){
        //             apply_force(**i,**j, &dmin, &davg, &navg);
        //         }
        //     }
        //     b_b++;
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
