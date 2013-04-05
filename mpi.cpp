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


int belongs_on(int num_procs, double ymin[], double ymax[], particle_t &p){
    for(int r=0; r<num_procs; r++){
        if(p.y <= ymax[r] && p.y >= ymin[r]){
            return r;
        }
    }
    //should never happen, should be in one of them!!
    assert(false);
    return -1;
}

void recompute_offsets(int num_procs, int sizes[], int *offsets){
    offsets[0] = 0;
    for (int r=0; r<num_procs; r++){
        offsets[r+1] = offsets[r]+sizes[r];
    }
    return;
}

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
    FILE *fsave = savename && rank == 0 ? fopen(savename, "w") : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen(sumname, "a") : NULL;
    // FILE *fsave = savename ? fopen(savename, "a") : NULL;
    int saveoutput = 0;
    if(rank==0 && fsave){
        saveoutput = 1;
    }
    MPI_Bcast(&saveoutput, 1, MPI_INT, 0, MPI_COMM_WORLD);
    DBOUT(printf("saveoutput = %d \n" ,saveoutput));
    
    
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

    particle_t **particles_to_send = (particle_t**) malloc(n_proc * sizeof(particle_t*));
    for(int r = 0; r < n_proc; r++){
        particles_to_send[r] = (particle_t*) malloc(n * sizeof(particle_t));
    }

    particle_t **particles_to_recv = (particle_t**) malloc(n_proc * sizeof(particle_t*));
    for(int r = 0; r < n_proc; r++){
        particles_to_recv[r] = (particle_t*) malloc(n * sizeof(particle_t));
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

    double ymax[n_proc];
    double ymin[n_proc];
    // ymin[0] = 0;
    for(int i = 0; i < n_proc; i++) {
        ymin[i] = i*(size/n_proc);
        ymax[i] = (i+1)*(size/n_proc);
    }

    if(rank==0)
        printf("size %f\n",size);
    for(int i = 0; i < n_proc; i++) {
        printf("%d min %f, max %f\n",i, ymin[i], ymax[i]);
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
            for(int j=partition_offsets[r]; j<partition_offsets[r+1]; j++){
                for(int i=0; i<blocksize; i++){
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
    // particle_t *myparticles = (particle_t*) malloc((myparticle_count)* sizeof(particle_t) );
    particle_t *myparticles = (particle_t*) malloc((n)* sizeof(particle_t) );
    
    printf("scatterv\n");
    // MPI_Scatterv(contig_block_particles, partition_sizes, partition_offsets, BLOCK,
    //              local_blocks, num_local_blocks, PARTICLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(contig_particles, particles_on_proc, particles_on_proc_offsets, PARTICLE, myparticles, myparticle_count, PARTICLE, 0, MPI_COMM_WORLD);

    printf("scatteredv\n");

    // for(int i=0; i<myparticle_count; i++){
    //     if(rank != belongs_on(n_proc, ymin, ymax, myparticles[i])){
    //         printf("rank %d,y %f, myymin %f, ymax %f , i %d/%d\n", rank, myparticles[i].y, ymin[rank], ymax[rank], i, myparticle_count);
    //     }
    //     assert(belongs_on(n_proc, ymin, ymax, myparticles[i]) == rank);
    // }

    
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
        // printf("step: %d\n", step);
        navg = 0;
        dmin = 1.0;
        davg = 0.0;

        /*****************************************************************************************************
         *
         * Begin particle communication
         *
         *****************************************************************************************************/ 
        int number_to_send_to[n_proc];
        int number_to_recv_from[n_proc];
        for(int r=0; r<n_proc; r++){
            number_to_send_to[r] = 0;
            number_to_recv_from[r] = 0;
        }
        

        for(int p = myparticle_count-1; p>=0  ; p--) {
            int this_belongs_on = belongs_on(n_proc, ymin, ymax, myparticles[p]);
            if(this_belongs_on != rank){
                // printf("Oh no particle doesnt belong!\n");
                // printf("We are %d, this belongs on %d\n",rank, this_belongs_on);
                particles_to_send[this_belongs_on][number_to_send_to[this_belongs_on]] = myparticles[p];
                number_to_send_to[this_belongs_on]++;
                // copy the end of the array to this position, should
                // work even if this one is the end;
                myparticle_count--;
                myparticles[p] = myparticles[myparticle_count];
            }
            // assert(this_belongs_on == rank);
            // assert(y <= ymax[rank] && y >= ymin[rank]);
        }
        // printf("big communication\n");
        for(int r=0; r<n_proc; r++){
            if(r != rank){
                // printf("rank: %d sending to %d \n", rank, r);
                MPI_Send(number_to_send_to+r, 1, MPI_INT, r, 1, MPI_COMM_WORLD);
            }
            else{
                for(int recv=0; recv<n_proc; recv++){
                    if(recv == rank){
                        continue;
                    }
                    // printf("rank: %d recving from %d \n", rank, recv);
                    MPI_Recv(number_to_recv_from+recv, 1, MPI_INT, recv, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }

        }
        

    
        for(int r=0; r<n_proc; r++){
            if(r != rank){
                if(number_to_send_to[r] > 0){
                    // printf("step %d rank: %d sending %d particles to %d \n", step, rank, number_to_send_to[r], r);
                    MPI_Send(particles_to_send[r], number_to_send_to[r], PARTICLE, r, 2, MPI_COMM_WORLD);
                }
            }
            else{
                for(int recv=0; recv<n_proc; recv++){
                    if(recv == rank){
                        continue;
                    }
                    if(number_to_recv_from[recv] > 0){
                        // printf("step %d rank: %d receiving %d particles from %d \n", step, rank, number_to_recv_from[recv], recv);
                        MPI_Recv(particles_to_recv[recv], number_to_recv_from[recv], PARTICLE, recv, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
            }
        }


        for(int r=0; r<n_proc; r++){
            for(int i=0; i<number_to_recv_from[r]; i++){
                double y = particles_to_recv[r][i].y;
                assert(y <= ymax[rank] && y >= ymin[rank]);
                myparticles[myparticle_count] = particles_to_recv[r][i];
                myparticle_count++;
                
            }
        }
        for(int i=0; i<myparticle_count; i++){
            // printf("my rank %d, thisbelongs on %d\n",rank,belongs_on(n_proc, ymin, ymax, myparticles[i]));
            if(rank != belongs_on(n_proc, ymin, ymax, myparticles[i])){
                printf("rank %d,y %f, myymin %f, ymax %f , i %d/%d\n", rank, myparticles[i].y, ymin[rank], ymax[rank], i, myparticle_count);
            }
            assert(belongs_on(n_proc, ymin, ymax, myparticles[i]) == rank);
        }

        int total_particle_count = 0;
        MPI_Reduce(&myparticle_count, &total_particle_count, 1, MPI_INT, 
                   MPI_SUM, 0, MPI_COMM_WORLD);
        if(rank == 0){
            assert(total_particle_count == n);
            // printf("total %d\n", total_particle_count);
        }
        /*****************************************************************************************************
         *
         * End particle communication
         *
         *****************************************************************************************************/



        for(int b=0; b<blocksize*blocksize; b++){
            number_in_block[b] = 0; // starts with no particles in any box;
        }
        for(int p = 0; p < myparticle_count; p++) {
            double x = myparticles[p].x;
            double y = myparticles[p].y;
            // assert(y <= ymax[rank] && y >= ymin[rank]);
            int x_index = (int)floor(x/block_width);
            int y_index = (int)floor(y/block_width);
            int particle_index = number_in_block[x_index + y_index*blocksize]++;
            blocks[x_index + y_index*blocksize][particle_index] = myparticles+p;
        }

        // for(int j=0; j<blocksize; j++){
        //     printf("\nrank, j: %d, %d numinblock:",rank,j);
        //     for(int i=0; i<blocksize; i++){
        //         printf(" %d", number_in_block[i + j*blocksize]);
        //     }
        // }
        //
        //  save current step if necessary (slightly different semantics than in other codes)
        //
        if(find_option(argc, argv, "-no") == -1){
            // if(fsave && (step % SAVEFREQ) == 0) {
            if( saveoutput && ( step % SAVEFREQ) == 0) {
                MPI_Gather(&myparticle_count, 1, MPI_INT, 
                           particles_on_proc, 1, MPI_INT, 
                           0, MPI_COMM_WORLD);
                recompute_offsets(n_proc, particles_on_proc, particles_on_proc_offsets);
                MPI_Gatherv(myparticles, myparticle_count, PARTICLE, particles, particles_on_proc, particles_on_proc_offsets, PARTICLE, 0, MPI_COMM_WORLD);
                if(rank == 0){
                    save(fsave, n, particles);
                    // save(fsave, myparticle_count, myparticles);
                }
                // for(int r=0; r<n_proc; r++){
                //     if(rank ==r){
                //         static bool first = true;
                //         if(first && rank==0) {
                //             fprintf(fsave, "%d %g\n", n, size);
                //             first = false;
                //         }
                //         for(int i = 0; i < myparticle_count; i++) {
                //             fprintf(fsave, "%g %g\n", myparticles[i].x, myparticles[i].y);
                //         }
                //         // save(fsave, myparticle_count, myparticles);
                //     }
                //     MPI_Barrier(MPI_COMM_WORLD );
                // }
                
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
