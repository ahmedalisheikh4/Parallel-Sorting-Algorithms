#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void heapify(int *arr, int size, int i) {
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;

    if (left < size && arr[left] > arr[largest]) {
        largest = left;
    }

    if (right < size && arr[right] > arr[largest]) {
        largest = right;
    }

    if (largest != i) {
        int temp = arr[i];
        arr[i] = arr[largest];
        arr[largest] = temp;

        heapify(arr, size, largest);
    }
}

void heap_sort(int *arr, int size) {
    for (int i = size / 2 - 1; i >= 0; i--) {
        heapify(arr, size, i);
    }

    for (int i = size - 1; i > 0; i--) {
        int temp = arr[0];
        arr[0] = arr[i];
        arr[i] = temp;

        heapify(arr, i, 0);
    }
}

void parallel_heap_sort(int *local_data, int local_size) {
    heap_sort(local_data, local_size);
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
         //   printf("%d ", data[i]);
        //}
        printf("\n");
    }

    int local_size = data_size / num_procs;
    int remainder = data_size % num_procs;

    int *local_data = (int *)malloc((local_size + (rank == num_procs - 1 ? remainder : 0)) * sizeof(int));
    start_time = MPI_Wtime(); // Start time measurement	
    MPI_Scatter(data, local_size, MPI_INT, local_data, local_size, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        MPI_Scatter(&data[data_size - remainder], remainder, MPI_INT, &local_data[local_size], remainder, MPI_INT, 0, MPI_COMM_WORLD);
    }


    parallel_heap_sort(local_data, local_size + (rank == num_procs - 1 ? remainder : 0));


    int *sorted_data = (int *)malloc(data_size * sizeof(int));
    MPI_Gather(local_data, local_size + (rank == num_procs - 1 ? remainder : 0), MPI_INT, sorted_data, local_size + (rank == num_procs - 1 ? remainder : 0), MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        heap_sort(sorted_data, data_size); // Final sorting of merged data
        end_time = MPI_Wtime(); // End time measurement
        printf("Sorted data (Heap Sort): ");
        for (int i = 0; i < data_size; i++) {
            printf("%d ", sorted_data[i]);
        }
        printf("\n");

        printf("Time taken for parallel heap sort: %f seconds\n", end_time - start_time);

        free(data);
        free(sorted_data);
    }

    free(local_data);
    MPI_Finalize();
    return 0;
}

