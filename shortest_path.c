#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "structs.h"
#define ROOT 0

void setup_grid(GRID_INFO_TYPE* grid) {
    int old_rank;
    int dimensions[2];
    int periods[2];
    int coordinates[2];
    int varying_coords[2];

    // Global grid information

    MPI_Comm_size(MPI_COMM_WORLD, &(grid->p));
    MPI_Comm_rank(MPI_COMM_WORLD, &old_rank);
    grid->q = (int)sqrt((double)grid->p);
    dimensions[0] = dimensions[1] = grid->q;
    periods[0] = periods[1] = 1;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, periods, 1, &(grid->comm));
    MPI_Comm_rank(grid->comm, &(grid->my_rank));
    MPI_Cart_coords(grid->comm, grid->my_rank, 2, coordinates);
    grid->my_row = coordinates[0];
    grid->my_col = coordinates[1];

    // Set up the row communicators

    varying_coords[0] = 0;
    varying_coords[1] = 1;
    MPI_Cart_sub(grid->comm, varying_coords, &(grid->row_comm));

    // Set up the column communicators

    varying_coords[0] = 1;
    varying_coords[1] = 0;
    MPI_Cart_sub(grid->comm, varying_coords, &(grid->col_comm));
}

// implementation of the fox algorithm
void fox() {}

/**
 * @brief Check if the Fox algorithm can be applied
 *
 * @param nprocess number of processes
 * @param nodes number of nodes
 * @return int 1 if true 0 otherwise
 */
int check_fox(int nprocess, int nodes) {
    int q = sqrt(nprocess);
    if (nodes % q != 0) {
        printf("Algorithm cannot be applied\n");
        return 0;
    }
    return 1;
}
/**
 * @brief Allocates memory for a matrix
 *
 * @param m matrix
 * @param n number of nodes
 */
int** allocate_memory(int** m, int n) {
    m = (int**)malloc(sizeof(int) * n * n);
    for (int i = 0; i < n; i++) {
        m[i] = (int*)malloc(sizeof(int) * n);
    }

    return m;
}
void recieve_input(int** m, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            scanf("%d", &m[i][j]);
        }
    }
}

void print_matrix(int** m, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%d ", m[i][j]);
        }
        printf("\n");
    }
}

int main(int argc, char* argv[]) {
    int my_rank, nprocess, nodes, **matrix;
    double finish, start;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocess);

    if (my_rank == ROOT) {
        scanf("%d", &nodes);
        if (check_fox(nprocess, nodes) == 0) {
            MPI_Finalize();
            exit(1);
        } else {
            matrix = allocate_memory(matrix, nodes);
            recieve_input(matrix, nodes);
            print_matrix(matrix, nodes);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    /* CODIGO PARA RESOLVER O PROBLEMA*/

    MPI_Barrier(MPI_COMM_WORLD);
    finish = MPI_Wtime();

    if (my_rank == ROOT) {
        printf("Execution time: %lf\n", finish - start);
    }

    MPI_Finalize();
    return 0;
}