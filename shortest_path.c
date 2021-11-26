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
void matrix_multiply(int n, int** m1, int** m2, int** m3) {
    int i, j, k;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                if ((m1[i][k] + m2[k][j]) < m3[i][j])
                    m3[i][j] = m1[i][k] + m2[k][j];
            }
        }
    }
}

void set_to_zero(int** m, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            m[i][j] = 0;
        }
    }
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

// implementation of the fox algorithm
void fox(int n, GRID_INFO_TYPE* grid, int** local_A, int** local_B, int** local_C) {
    int step, bcast_root, n_bar, source, dest, tag = 43, **temp_A;
    MPI_Status status;

    n_bar = n / grid->q;
    /*CRIAR FUNCAO*/
    set_to_zero(local_C, n);

    source = (grid->my_row + 1) % grid->q;
    dest = (grid->my_row + grid->q - 1) % grid->q;

    temp_A = allocate_memory(temp_A, n_bar);

    for (step = 0; step < grid->q; step++) {
        bcast_root = (grid->my_row + step) % grid->q;
        if (bcast_root == grid->my_col) {
            MPI_Bcast(local_A, 1, MPI_INT, bcast_root, grid->row_comm);
            matrix_multiply(n, local_A, local_B, local_C);
        } else {
            MPI_Bcast(temp_A, 1, MPI_INT, bcast_root, grid->row_comm);
            matrix_multiply(n, temp_A, local_B, local_C);
        }
        MPI_Send(local_B, 1, MPI_INT, dest, tag, grid->col_comm);
        MPI_Recv(local_B, 1, MPI_INT, source, tag, grid->col_comm, &status);
    }
}

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
    int my_rank, nprocess, nodes, **matrix, valid;
    double finish, start;
    GRID_INFO_TYPE grid;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocess);

    if (my_rank == ROOT) {
        scanf("%d", &nodes);
        matrix = allocate_memory(matrix, nodes);
        recieve_input(matrix, nodes);
        print_matrix(matrix, nodes);
        valid = check_fox(nprocess, nodes);
    }

    MPI_Bcast(&valid, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    if (valid == 0) {
        MPI_Finalize();
        return 0;
    }

    setup_grid(&grid);

    //share the number of nodes with everyone
    MPI_Bcast(&nodes, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

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