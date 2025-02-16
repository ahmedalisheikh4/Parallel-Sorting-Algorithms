#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void selection_sort(int *arr, int size) {
    for (int i = 0; i < size - 1; i++) {
        int min_index = i;
        for (int j = i + 1; j < size; j++) {
            if (arr[j] < arr[min_index]) {
                min_index = j;
            }
        }
        if (min_index != i) {
            int temp = arr[i];
            arr[i] = arr[min_index];
            arr[min_index] = temp;
        }
    }
}

void parallel_selection_sort(int *local_data, int local_size) {
    for (int i = 0; i < local_size - 1; i++) {
        int min_index = i;
        for (int j = i + 1; j < local_size; j++) {
            if (local_data[j] < local_data[min_index]) {
                min_index = j;
            }
        }
        if (min_index != i) {
            int temp = local_data[i];
            local_data[i] = local_data[min_index];
            local_data[min_index] = temp;
        }
    }
}

void merge_sorted_arrays(int *data, int *local_data, int local_size, int data_size, int rank) {
    int *temp = (int *)malloc(data_size * sizeof(int));
    MPI_Gather(local_data, local_size, MPI_INT, temp, local_size, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        selection_sort(temp, data_size);
        for (int i = 0; i < data_size; i++) {
            data[i] = temp[i];
        }
    }

    free(temp);
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
        //printf("Original data: ");
        //for (int i = 0; i < data_size; i++) {
          //  printf("%d ", data[i]);
        //}
        printf("\n");
    }

    MPI_Barrier(MPI_COMM_WORLD); // Synchronize before sorting

    int local_size = data_size / num_procs;
    int *local_data = (int *)malloc(local_size * sizeof(int));

    double start_time = MPI_Wtime(); // Start timing
    MPI_Scatter(data, local_size, MPI_INT, local_data, local_size, MPI_INT, 0, MPI_COMM_WORLD);

    parallel_selection_sort(local_data, local_size);
    merge_sorted_arrays(data, local_data, local_size, data_size, rank); // Pass 'rank'
    double end_time = MPI_Wtime(); // End timing

    if (rank == 0) {
        printf("Sorted data (Selection Sort): ");
        for (int i = 0; i < data_size; i++) {
            printf("%d ", data[i]);
        }
        printf("\n");
        free(data);
    }

    free(local_data);

    MPI_Barrier(MPI_COMM_WORLD); // Synchronize after sorting

    if (rank == 0) {
        printf("Time taken for sorting: %f seconds\n", end_time - start_time);
    }

    MPI_Finalize();
    return 0;
}

