#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void swap(int *a, int *b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

int partition(int *arr, int low, int high) {
    int pivot = arr[high];
    int i = low - 1;

    for (int j = low; j < high; j++) {
        if (arr[j] <= pivot) {
            i++;
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return i + 1;
}

void quick_sort(int *arr, int low, int high) {
    if (low < high) {
        int pi = partition(arr, low, high);

        quick_sort(arr, low, pi - 1);
        quick_sort(arr, pi + 1, high);
    }
}

void parallel_quick_sort(int *local_data, int local_size, int rank, int num_procs) {
    quick_sort(local_data, 0, local_size - 1);

    int *recv_data = (int *)malloc(local_size * num_procs * sizeof(int));

    MPI_Allgather(local_data, local_size, MPI_INT, recv_data, local_size, MPI_INT, MPI_COMM_WORLD);

    int *sorted_data = (int *)malloc(local_size * num_procs * sizeof(int));

    for (int i = 0; i < local_size * num_procs; i++) {
        sorted_data[i] = recv_data[i];
    }

    quick_sort(sorted_data, 0, local_size * num_procs - 1);

    for (int i = 0; i < local_size; i++) {
        local_data[i] = sorted_data[rank * local_size + i];
    }

    free(recv_data);
    free(sorted_data);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    double start_time, end_time;

    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    int data_size = (100 * 1024 * 1024)/(sizeof(int));  // 100mb
    int *data = NULL;

    if (rank == 0) {
        data = (int *)malloc(data_size * sizeof(int));
        for (int i = 0; i < data_size; i++) {
            data[i] = rand() % 100; // Generate random data
        }
        //printf("Original data: ");
        //for (int i = 0; i < data_size; i++) {
            //printf("%d ", data[i]);
        //}
        printf("\n");
    }

    int local_size = data_size / num_procs;
    int *local_data = (int *)malloc(local_size * sizeof(int));
start_time = MPI_Wtime(); // Record starting time
    MPI_Scatter(data, local_size, MPI_INT, local_data, local_size, MPI_INT, 0, MPI_COMM_WORLD);

    

    parallel_quick_sort(local_data, local_size, rank, num_procs);

    MPI_Gather(local_data, local_size, MPI_INT, data, local_size, MPI_INT, 0, MPI_COMM_WORLD);
end_time = MPI_Wtime(); // Record ending time
    if (rank == 0) {
        printf("Sorted data (Quick Sort): ");
        for (int i = 0; i < data_size; i++) {
            printf("%d ", data[i]);
        }
        printf("\n");

        // Calculate and print sorting time
        printf("Sorting time: %f seconds\n", end_time - start_time);

        free(data);
    }

    free(local_data);
    MPI_Finalize();
    return 0;
}
