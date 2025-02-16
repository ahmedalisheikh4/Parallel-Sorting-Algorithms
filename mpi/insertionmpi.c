#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void insertion_sort(int *arr, int size) {
    for (int i = 1; i < size; i++) {
        int key = arr[i];
        int j = i - 1;

        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}

void parallel_insertion_sort(int *local_data, int local_size) {
    insertion_sort(local_data, local_size);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    int data_size = (1 * 1024 * 1024)/(sizeof(int));  // 100mb 
    int *data = NULL;

    if (rank == 0) {
        data = (int *)malloc(data_size * sizeof(int));
        for (int i = 0; i < data_size; i++) {
            data[i] = rand() % 100; // Generate random data
        }
    }

    int local_size = data_size / num_procs;
    int *local_data = (int *)malloc(local_size * sizeof(int));
double start_time = MPI_Wtime(); // Start timing
    MPI_Scatter(data, local_size, MPI_INT, local_data, local_size, MPI_INT, 0, MPI_COMM_WORLD);



    parallel_insertion_sort(local_data, local_size);

    

    MPI_Gather(local_data, local_size, MPI_INT, data, local_size, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
      
        insertion_sort(data, data_size); // Final sorting after gathering
        double end_time = MPI_Wtime(); // End timing

        printf("Sorted data (Insertion Sort): ");
        for (int i = 0; i < data_size; i++) {
            printf("%d ", data[i]);
        }
        printf("\n");
        free(data);

        printf("Sorting time: %f seconds\n", end_time - start_time);
    }

    free(local_data);
    MPI_Finalize();
    return 0;
}

