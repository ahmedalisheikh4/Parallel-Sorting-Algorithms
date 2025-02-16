#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int getNextGap(int gap) {
    gap = (gap * 10) / 13;
    if (gap < 1) {
        return 1;
    }
    return gap;
}

void comb_sort(int *arr, int size) {
    int gap = size;
    int swapped = 1;

    while (gap != 1 || swapped == 1) {
        gap = getNextGap(gap);
        swapped = 0;

        for (int i = 0; i < size - gap; i++) {
            if (arr[i] > arr[i + gap]) {
                int temp = arr[i];
                arr[i] = arr[i + gap];
                arr[i + gap] = temp;
                swapped = 1;
            }
        }
    }
}

void parallel_comb_sort(int *local_data, int local_size, int rank, int num_procs, int data_size) {
    // Sort local data segment
    comb_sort(local_data, local_size);

    int *sorted_data = NULL;

    if (rank == 0) {
        // Allocate memory for sorted data in the root process
        sorted_data = (int *)malloc(data_size * sizeof(int));
    }

    // Gather all local_data in sorted_data in rank 0
    MPI_Gather(local_data, local_size, MPI_INT, sorted_data, local_size, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // Sort the gathered data in the root process
        comb_sort(sorted_data, data_size);

        // Distribute sorted_data to all processes
        int *send_counts = (int *)malloc(num_procs * sizeof(int));
        int *displs = (int *)malloc(num_procs * sizeof(int));

        int remainder = data_size % num_procs;
        int quotient = data_size / num_procs;
        int offset = 0;

        for (int i = 0; i < num_procs; i++) {
            send_counts[i] = (i < remainder) ? quotient + 1 : quotient;
            displs[i] = offset;
            offset += send_counts[i];
        }

        // Scatterv the sorted data to all processes
        MPI_Scatterv(sorted_data, send_counts, displs, MPI_INT, local_data, local_size, MPI_INT, 0, MPI_COMM_WORLD);

        // Free allocated memory
        free(send_counts);
        free(displs);
        free(sorted_data);
    }
}


int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    int data_size = (1 * 1024  * 1024)/(sizeof(int));  // 100mb
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

    // Start the timer
    double start_time = MPI_Wtime();

    int local_size = data_size / num_procs;
    int *local_data = (int *)malloc(local_size * sizeof(int));

    MPI_Scatter(data, local_size, MPI_INT, local_data, local_size, MPI_INT, 0, MPI_COMM_WORLD);

    parallel_comb_sort(local_data, local_size, rank, num_procs, data_size);

    // Gather all local_data into a single array on rank 0
    MPI_Gather(local_data, local_size, MPI_INT, data, local_size, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // Sort the fully gathered data on the root process
        comb_sort(data, data_size);

        // Stop the timer
        double end_time = MPI_Wtime();

        // Calculate sorting time
        double sorting_time = end_time - start_time;

        // The root process prints the fully sorted data and sorting time
        printf("Sorted data (Comb Sort): ");
        for (int i = 0; i < data_size; i++) {
            printf("%d ", data[i]);
        }
        printf("\nSorting Time: %f seconds\n", sorting_time);
        free(data);
    }

    free(local_data);
    MPI_Finalize();
    return 0;
}
