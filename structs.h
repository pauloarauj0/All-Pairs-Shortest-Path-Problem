// Struct of a grid based on the book "A User's guide to MPI" by Peter S.Pacheco
typedef struct {
    int p;              // Number of processes
    MPI_Comm comm;      // Communicator for the entire grid
    MPI_Comm row_comm;  // Communicator for my row
    MPI_Comm col_comm;  // Communicator for my col
    int q;              // Order of grid
    int my_row;
    int my_col;
    int my_rank;
} GRID_INFO_TYPE;