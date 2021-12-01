#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "structs.h"
#define ROOT 0
#define INF 99999
/**
 * @brief Allocates memory for a matrix
 *
 * @param m matrix
 * @param n number of nodes
 */
void allocate_memory(int*** m, int n) {
    int* p = (int*)malloc(n * n * sizeof(int*));
    (*m) = (int**)malloc(n * sizeof(int*));

    for (int i = 0; i < n; i++) {
        (*m)[i] = &(p[i * n]);
    }
}
void recieve_input(int** m, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            scanf("%d", &m[i][j]);
            // if m[i][j]==0 it means that there is no path
            // i == j its the distance between a node and itself
            if (i != j && (m[i][j] == 0)) m[i][j] = INF;
        }
    }
}

void print_matrix(int** m, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // remove the INF
            if (m[i][j] == INF) {
                m[i][j] = 0;
            }
            printf("%d ", m[i][j]);
        }
        printf("\n");
    }
}
void submatrix(int** m, GRID_INFO_TYPE* grid, int** matrix_aux, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int aux_i = grid->my_row * n + i;
            int aux_j = grid->my_col * n + j;
            matrix_aux[i][j] = m[aux_i][aux_j];
        }
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
                if ((m1[i][k] + m2[k][j]) < m3[i][j]) {
                    m3[i][j] = m1[i][k] + m2[k][j];
                }
            }
        }
    }
}

void set_to_inf(int** m, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            m[i][j] = INF;
        }
    }
}
int** transfer_matrix(int** m, int n) {
    int** aux;
    allocate_memory(&aux, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aux[i][j] = m[i][j];
        }
    }
    return aux;
}

// implementation of the fox algorithm
void fox(int n, GRID_INFO_TYPE* grid, int** local_A, int** local_B,
         int** local_C) {
    int step, bcast_root, n_bar, source, dest, tag = 43, **temp_A;
    MPI_Status status;

    set_to_inf(local_C, n);
    source = (grid->my_row + 1) % grid->q;
    dest = (grid->my_row + grid->q - 1) % grid->q;

    allocate_memory(&temp_A, n);

    for (step = 0; step < grid->q; step++) {
        bcast_root = (grid->my_row + step) % grid->q;
        if (bcast_root == grid->my_col) {
            MPI_Bcast(&local_A[0][0], n * n, MPI_INT, bcast_root,
                      grid->row_comm);
            matrix_multiply(n, local_A, local_B, local_C);
        } else {
            MPI_Bcast(&temp_A[0][0], n * n, MPI_INT, bcast_root,
                      grid->row_comm);
            matrix_multiply(n, temp_A, local_B, local_C);
        }
        MPI_Send(&local_B[0][0], n * n, MPI_INT, dest, tag, grid->col_comm);
        MPI_Recv(&local_B[0][0], n * n, MPI_INT, source, tag, grid->col_comm,
                 &status);
    }
}

int main(int argc, char* argv[]) {
    int my_rank, nprocess, nodes, **matrix, valid, half_size;
    double finish, start;
    GRID_INFO_TYPE grid;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocess);

    if (my_rank == ROOT) {
        scanf("%d", &nodes);
        valid = check_fox(nprocess, nodes);
    }

    MPI_Bcast(&valid, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    // Close the program if the input is not valid
    if (valid == 0) {
        MPI_Finalize();
        return 0;
    }

    // Prepare to recieve the input
    MPI_Bcast(&nodes, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    allocate_memory(&matrix, nodes);

    setup_grid(&grid);
    half_size = nodes / grid.q;

    // Root reads the input
    if (my_rank == ROOT) {
        recieve_input(matrix, nodes);
        // print_matrix(matrix, nodes);
    }
    // MPI_Bcast(&(matrix[0][0]), nodes * nodes, MPI_INT, ROOT, MPI_COMM_WORLD);
    int **m1, **res;
    allocate_memory(&m1, half_size);
    allocate_memory(&res, half_size);
    MPI_Scatter(&matrix[0][0], half_size * half_size, MPI_INT, &m1[0][0],
                half_size * half_size, MPI_INT, ROOT, MPI_COMM_WORLD);

    // submatrix(matrix, &grid, m1, half_size);

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    /* CODIGO PARA RESOLVER O PROBLEMA*/

    for (int i = 1; i < nodes - 1; i *= 2) {
        fox(half_size, &grid, m1, m1, res);
        m1 = transfer_matrix(res, half_size);
    }

    // printf("Print do %d\n", my_rank);
    print_matrix(res, half_size);
    MPI_Barrier(MPI_COMM_WORLD);
    finish = MPI_Wtime();

    if (my_rank == ROOT) {
        // printf("Execution time: %lf\n", finish - start);
    }
    MPI_Finalize();
    return 0;
}