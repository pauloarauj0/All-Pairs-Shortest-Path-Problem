//#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#define ROOT 0

int main() {
    int my_rank, nprocess, nodes;
    // MPI_Init(&argc, &argv);
    // MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    // MPI_Comm_size(MPI_COMM_WORLD, &nprocess);

    scanf("%d", &nodes);
    int vertices[nodes][nodes];

    for (int i = 0; i < nodes; i++) {
        for (int j = 0; j < nodes; j++) {
            scanf("%d ", &vertices[i][j]);
        }
    }

    printf("\n");
    for (int i = 0; i < nodes; i++) {
        for (int j = 0; j < nodes; j++) {
            printf("%d ", vertices[i][j]);
        }
        printf("\n");
    }

    //MPI_Finalize();
    return 0;
}