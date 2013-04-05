#ifdef RELEASE
#define DBOUT( x )
#else
#define DBOUT( x )  if(rank == 0) x
#endif


#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include "common.h"
#define MAXPARTILCESPERBOX 8



//
//  benchmarking program
//
int main(int argc, char **argv)
{
    int navg, nabsavg = 0;
    double dmin, absmin = 1.0, davg, absavg = 0.0;
    double rdavg, rdmin;
    int rnavg;

    //
    //  process command line parameters
    //
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

    //
    //  set up MPI
    //
    int n_proc, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //
    //  allocate generic resources
    //
    // FILE *fsave = savename && rank == 0 ? fopen(savename, "w") : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen(sumname, "a") : NULL;
    FILE *fsave = savename ? fopen(savename, "a") : NULL;


    particle_t *particles = (particle_t*) malloc(n * sizeof(particle_t));
    double size =  sqrt(0.0005 * n);
    double interaction_length = 0.01;

    int blocksize = (int) ceil(size/interaction_length);
    DBOUT(printf("blocksize = %d \n" ,blocksize));

    double block_width = size/blocksize;

    particle_t*** blocks = (particle_t***) malloc(blocksize*blocksize * sizeof(particle_t**));
    for(int b=0; b<blocksize*blocksize; b++){
        blocks[b] = (particle_t**)malloc(MAXPARTILCESPERBOX*sizeof(particle_t*));
    }
    int number_in_block[blocksize*blocksize];
    for(int b=0; b<blocksize*blocksize; b++){
        number_in_block[b] = 0; // starts with no particles in any box;
    }

    
    DBOUT(printf("setup particle mpi types\n"));
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous(6, MPI_DOUBLE, &PARTICLE);
    MPI_Type_commit(&PARTICLE);

    DBOUT(printf("setup block mpi types\n"));
    MPI_Datatype BLOCK;
    MPI_Type_contiguous(6*MAXPARTILCESPERBOX, MPI_DOUBLE, &BLOCK);
    MPI_Type_commit(&BLOCK);

    //
    //  set up the data partitioning across processors
    //
    DBOUT(printf("setup partition offsets \n"));
    int blocks_per_proc = (blocksize + n_proc - 1) / n_proc;
    int *partition_offsets = (int*) malloc((n_proc + 1) * sizeof(int));
    for(int i = 0; i < n_proc + 1; i++) {
        DBOUT(printf("%d \n" ,i));
        partition_offsets[i] = min(i * blocks_per_proc, blocksize);
        DBOUT(printf("partition_ofsets[%d] = %d \n" ,i, partition_offsets[i]));
    }

    DBOUT(printf("setup partition sizes \n"));
    int *partition_sizes = (int*) malloc(n_proc * sizeof(int));
    for(int i = 0; i < n_proc; i++) {
        DBOUT(printf("%d \n" ,i));
        partition_sizes[i] = partition_offsets[i + 1] - partition_offsets[i];
        DBOUT(printf("partition_sizes[%d] = %d \n" ,i, partition_sizes[i]));
    }

    
    //allocate storage for local partition
    
    // DBOUT(printf("setup local blocks\n"));
    // int num_local_blocks = partition_sizes[rank];

    // particle_t *local_blocks = (particle_t*) malloc(num_local_blocks*MAXPARTILCESPERBOX * sizeof(particle_t));
    
    
    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //

    set_size(n);
    if(rank == 0) {
        DBOUT(printf("init particles\n"));
        init_particles(n, particles);


        DBOUT(printf("set blocks the first time\n"));
        for(int p = 0; p < n; p++) {
            double x = particles[p].x;
            double y = particles[p].y;

            int x_index = (int)floor(x/block_width);
            int y_index = (int)floor(y/block_width);
            int particle_index = number_in_block[x_index + y_index*blocksize]++;
            blocks[x_index + y_index*blocksize][particle_index] = particles+p;
        }

    }
    particle_t *contig_particles = (particle_t*) malloc(n * sizeof(particle_t));

    int particles_on_proc[n_proc];
    int particles_on_proc_offsets[n_proc+1];
    particles_on_proc_offsets[0] = 0;
    DBOUT(printf("set contig blocks\n"));
    if(rank == 0) {
        int contig_index = 0;
        for (int r=0; r<n_proc; r++){
            particles_on_proc[r] = 0;
            for(int i=partition_offsets[r]; i<partition_offsets[r+1]; i++){
                for(int j=0; j<blocksize; j++){
                    for(int p=0; p<number_in_block[i + j*blocksize]; p++ ){
                        // for(int p=0; p<number_in_block[i + j*blocksize]; p++ ){
                        contig_particles[contig_index] = *(blocks[i + j*blocksize][p]);
                        contig_index++;
                        particles_on_proc[r]++;
                    }
                }
                
            }
            particles_on_proc_offsets[r+1] = particles_on_proc_offsets[r]+particles_on_proc[r];
        }
        assert(contig_index == n);
        int total_particles = 0;
        for (int r=0; r<n_proc; r++){
            total_particles += particles_on_proc[r];
            DBOUT(printf("partciles on proc[%d] = %d\n",r,particles_on_proc[r]));
        }
        assert(total_particles == n);
    }
    MPI_Bcast(particles_on_proc, n_proc, MPI_INT, 0, MPI_COMM_WORLD);

    int myparticle_count = particles_on_proc[rank];
    particle_t *myparticles = (particle_t*) malloc((myparticle_count)* sizeof(particle_t) );
    
    printf("scatterv\n");
    // MPI_Scatterv(contig_block_particles, partition_sizes, partition_offsets, BLOCK,
    //              local_blocks, num_local_blocks, PARTICLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(contig_particles, particles_on_proc, particles_on_proc_offsets, PARTICLE, myparticles, myparticle_count, PARTICLE, 0, MPI_COMM_WORLD);

    printf("scatteredv\n");

    // std::cout << rank << " " << myparticle_count << std::endl;
    // for(int p=0; p<myparticle_count; p++){
    //     std::cout << rank << " "<<myparticles[p].x << " " << myparticles[p].y<< std::endl;
    // }
    // exit(1);
    // printf("inital scatter\n");
    // MPI_Scatter(particles, n, PARTICLE, particles, n, PARTICLE, 0, MPI_COMM_WORLD);
    // printf("inital scattered\n");

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer();
    for(int step = 0; step < NSTEPS; step++) {
        navg = 0;
        dmin = 1.0;
        davg = 0.0;
        //
        //  collect all global data locally (not good idea to do)
        //
        // MPI_Allgatherv(local, nlocal, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD);

        for(int b=0; b<blocksize*blocksize; b++){
            number_in_block[b] = 0; // starts with no particles in any box;
        }
        for(int p = 0; p < myparticle_count; p++) {
            double x = myparticles[p].x;
            double y = myparticles[p].y;
                
            int x_index = (int)floor(x/block_width);
            int y_index = (int)floor(y/block_width);
            int particle_index = number_in_block[x_index + y_index*blocksize]++;
            blocks[x_index + y_index*blocksize][particle_index] = myparticles+p;
        }
        

        //
        //  save current step if necessary (slightly different semantics than in other codes)
        //
        if(find_option(argc, argv, "-no") == -1)
            // if(fsave && (step % SAVEFREQ) == 0) {
            if( (step % SAVEFREQ) == 0) {
                // MPI_Scatterv(contig_particles, particles_on_proc, particles_on_proc_offsets, PARTICLE, myparticles, myparticle_count, PARTICLE, 0, MPI_COMM_WORLD);
                
                // MPI_Gatherv(myparticles, myparticle_count, PARTICLE, particles, particles_on_proc, particles_on_proc_offsets, PARTICLE, 0, MPI_COMM_WORLD);
                // save(fsave, n, particles);
                for(int r=0; r<n_proc; r++){
                    if(rank ==r)
                        save(fsave, myparticle_count, myparticles);
                    MPI_Barrier(MPI_COMM_WORLD );
                }
                
            }

        //
        //  compute all forces
        //
        for(int i=0; i<blocksize; i++){
            for(int j=0; j<blocksize; j++){
                for(int p=0; p<number_in_block[i + j*blocksize]; p++ ){
                    blocks[i + j*blocksize][p]->ax = blocks[i + j*blocksize][p]->ay = 0;
                    //interact with blocks and neighbors
                    for(int xoffset=-1; xoffset<2; xoffset++){
                        int xblockindex = i+xoffset;
                        // dont go out of bounds
                        if(xblockindex<0 || xblockindex >= blocksize ){
                            continue;
                        }
                        for(int yoffset=-1; yoffset<2; yoffset++){
                            int yblockindex = j+yoffset;
                            // dont go out of bounds
                            if(yblockindex<0 || yblockindex >= blocksize ){
                                continue;
                            }
                            for(int num=0; num<number_in_block[xblockindex + yblockindex*blocksize]; num++ ){
                                // apply the force
                                apply_force(*(blocks[i + j*blocksize][p]), // this particle
                                            *(blocks[xblockindex + yblockindex*blocksize][num]), // its neighbor
                                            &dmin, &davg, &navg);
                            }
                        }
                    }
                }
            }
        }
        // for(int i = 0; i < nlocal; i++) {
        //     local[i].ax = local[i].ay = 0;
        //     for(int j = 0; j < n; j++) {
        //         apply_force(local[i], particles[j], &dmin, &davg, &navg);
        //     }
        // }

        if(find_option(argc, argv, "-no") == -1) {

            MPI_Reduce(&davg, &rdavg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&navg, &rnavg, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&dmin, &rdmin, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);


            if(rank == 0) {
                //
                // Computing statistical data
                //
                if(rnavg) {
                    absavg +=  rdavg / rnavg;
                    nabsavg++;
                }
                if(rdmin < absmin) {
                    absmin = rdmin;
                }
            }
        }

        //
        //  move particles
        //x
        for(int i = 0; i < myparticle_count; i++) {
            move(myparticles[i]);
        }
        // for(int i = 0; i < nlocal; i++) {
        //     move(local[i]);
        // }
    }
    simulation_time = read_timer() - simulation_time;

    if(rank == 0) {
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

        //
        // Printing summary data
        //
        if(fsum) {
            fprintf(fsum, "%d %d %g\n", n, n_proc, simulation_time);
        }
    }

    //
    //  release resources
    //
    if(fsum) {
        fclose(fsum);
    }
    // free(partition_offsets);
    // free(partition_sizes);
    // free(local);
    free(particles);
    if(fsave) {
        fclose(fsave);
    }

    MPI_Finalize();

    return 0;
}
