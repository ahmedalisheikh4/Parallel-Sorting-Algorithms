#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int compare(const void *a, const void *b) {
    return (*(int *)a - *(int *)b);
}

void SequentialSort(int *arr, int size) {
    qsort(arr, size, sizeof(int), compare);
}

void BitonicMerge(int *arr, int low, int cnt, int dir) {
    if (cnt > 1) {
        int k = cnt / 2;
        for (int i = low; i < low + k; i++) {
            if ((arr[i] > arr[i + k]) == dir) {
                int temp = arr[i];
                arr[i] = arr[i + k];
                arr[i + k] = temp;
            }
        }
        BitonicMerge(arr, low, k, dir);
        BitonicMerge(arr, low + k, k, dir);
    }
}

void BitonicSort(int *arr, int low, int cnt, int dir) {
    if (cnt > 1) {
        int k = cnt / 2;
        BitonicSort(arr, low, k, 1);
        BitonicSort(arr, low + k, k, 0);
        BitonicMerge(arr, low, cnt, dir);
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int data_size = (100 * 1024 * 1024)/(sizeof(int));  // 100mb
    int local_data_size = data_size / size;

    // Allocate memory for local data
    int *local_data = (int *)malloc(local_data_size * sizeof(int));

    // Generate random data for each process
    srand(rank); // Seed using rank for different random values
    for (int i = 0; i < local_data_size; i++) {
        local_data[i] = rand() % 100; // Generate random data
    }

    MPI_Barrier(MPI_COMM_WORLD); // Synchronize before starting timing
    double start_time = MPI_Wtime();

    // Perform local sorting using sequential quicksort
    SequentialSort(local_data, local_data_size);

    // Perform bitonic sorting
    BitonicSort(local_data, 0, local_data_size, 1); // 1 for ascending order

    MPI_Barrier(MPI_COMM_WORLD); // Synchronize after timing

    // Gather all local sorted data into one array on rank 0
    int *sorted_data = NULL;
    if (rank == 0) {
        sorted_data = (int *)malloc(data_size * sizeof(int));
    }

    MPI_Gather(local_data, local_data_size, MPI_INT, sorted_data, local_data_size, MPI_INT, 0, MPI_COMM_WORLD);

    // Merge all gathered sorted arrays on rank 0
    if (rank == 0) {
        SequentialSort(sorted_data, data_size); // Sort the merged array
        double end_time = MPI_Wtime();
        printf("Sorted data (Bitonic Sort): ");
        for (int i = 0; i < data_size; i++) {
            printf("%d ", sorted_data[i]);
        }
        printf("\n");

        printf("Time taken for sorting: %f seconds\n", end_time - start_time);
        free(sorted_data);
    }

    free(local_data);
    MPI_Finalize();
    return 0;
}

