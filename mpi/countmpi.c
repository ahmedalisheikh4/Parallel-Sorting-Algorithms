#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void count_sort(int *arr, int size, int range) {
    int count[range + 1];
    int *output = (int *)malloc(size * sizeof(int));

    for (int i = 0; i <= range; i++) {
        count[i] = 0;
    }

    for (int i = 0; i < size; i++) {
        count[arr[i]]++;
    }

    for (int i = 1; i <= range; i++) {
        count[i] += count[i - 1];
    }

    for (int i = size - 1; i >= 0; i--) {
        output[count[arr[i]] - 1] = arr[i];
        count[arr[i]]--;
    }

    for (int i = 0; i < size; i++) {
        arr[i] = output[i];
    }

    free(output);
}

void parallel_count_sort(int *local_data, int local_size, int range) {
    count_sort(local_data, local_size, range);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
   
    int data_size = (100 * 1024 * 1024)/(sizeof(int));  // 100mb
    int range = 100;    // Modify this according to the range of integers in your data

    int *data = NULL;

    if (rank == 0) {
        data = (int *)malloc(data_size * sizeof(int));
        for (int i = 0; i < data_size; i++) {
            data[i] = rand() % range; // Generate random data within the specified range
        }
        //printf("Original data: ");
        //for (int i = 0; i < data_size; i++) {
          //  printf("%d ", data[i]);
       // }
        printf("\n");
    }

    int local_size = data_size / num_procs;
    int *local_data = (int *)malloc(local_size * sizeof(int));
    double start_time = MPI_Wtime();
    MPI_Scatter(data, local_size, MPI_INT, local_data, local_size, MPI_INT, 0, MPI_COMM_WORLD);

   // double start_time = MPI_Wtime(); // Start timing

    parallel_count_sort(local_data, local_size, range);

   

    int *sorted_data = NULL;
    if (rank == 0) {
        sorted_data = (int *)malloc(data_size * sizeof(int));
    }
 //double end_time = MPI_Wtime(); // End timing
    MPI_Gather(local_data, local_size, MPI_INT, sorted_data, local_size, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {

        count_sort(sorted_data, data_size, range); // Final sorting after gathering
        double end_time = MPI_Wtime(); // End timing
 
        printf("Sorted data (Count Sort): ");
        for (int i = 0; i < data_size; i++) {
            printf("%d ", sorted_data[i]);
        }
        printf("\n");
        free(sorted_data);
        free(data);

        printf("Sorting time: %f seconds\n", end_time - start_time);
    }

    free(local_data);
    MPI_Finalize();
    return 0;
}

