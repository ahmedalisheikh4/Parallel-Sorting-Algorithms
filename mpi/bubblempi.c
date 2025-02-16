#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void bubble_sort(int *data, int size) {
    int temp;
    for (int i = 0; i < size - 1; i++) {
        for (int j = 0; j < size - i - 1; j++) {
            if (data[j] > data[j + 1]) {
                temp = data[j];
                data[j] = data[j + 1];
                data[j + 1] = temp;
            }
        }
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int *data = NULL;
    int data_size = (100 * 1024 * 1024)/(sizeof(int));  // 100mb
    if (rank == 0) {
        // Generate or input data to sort
        data = (int *)malloc(data_size * sizeof(int));
        // Replace with your data generation or input mechanism
        for (int i = 0; i < data_size; i++) {
            data[i] = rand() % 100; // Generate random data
        }
        //printf("Original data: ");
        //for (int i = 0; i < data_size; i++) {
         //   printf("%d ", data[i]);
       // }
        //printf("\n");
    }

    int *local_data = (int *)malloc(data_size / size * sizeof(int));
    double t1 = MPI_Wtime();
    MPI_Scatter(data, data_size / size, MPI_INT, local_data, data_size / size, MPI_INT, 0, MPI_COMM_WORLD);

    // Perform sorting algorithm (Bubble Sort)
    bubble_sort(local_data, data_size / size);

    int *sorted_data = NULL;
    if (rank == 0) {
        sorted_data = (int *)malloc(data_size * sizeof(int));
    }

    MPI_Gather(local_data, data_size / size, MPI_INT, sorted_data, data_size / size, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // Perform final sorting (Bubble Sort) on gathered data
        bubble_sort(sorted_data, data_size);
	double t2 = MPI_Wtime();
        printf("\nSorted data (Bubble Sort): ");
        for (int i = 0; i < data_size; i++) {
            printf("%d ", sorted_data[i]);
        }
        printf("\n");
        printf("Time taken: %f\n",t2-t1);
        free(data);
        free(sorted_data);
    }

    free(local_data);
    MPI_Finalize();
    return 0;
}

