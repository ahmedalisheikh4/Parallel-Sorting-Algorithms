#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void merge(int *arr, int l, int m, int r) {
    int n1 = m - l + 1;
    int n2 = r - m;

    int L[n1], R[n2];

    for (int i = 0; i < n1; i++) {
        L[i] = arr[l + i];
    }
    for (int j = 0; j < n2; j++) {
        R[j] = arr[m + 1 + j];
    }

    int i = 0, j = 0, k = l;

    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            i++;
        } else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }

    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }

    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }
}

void merge_sort(int *arr, int l, int r) {
    if (l < r) {
        int m = l + (r - l) / 2;

        merge_sort(arr, l, m);
        merge_sort(arr, m + 1, r);

        merge(arr, l, m, r);
    }
}

void parallel_merge_sort(int *local_data, int local_size) {
    merge_sort(local_data, 0, local_size - 1);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    double start_time, end_time;

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
 start_time = MPI_Wtime(); // Record starting time
    MPI_Scatter(data, local_size, MPI_INT, local_data, local_size, MPI_INT, 0, MPI_COMM_WORLD);

   

    parallel_merge_sort(local_data, local_size);

    

    MPI_Gather(local_data, local_size, MPI_INT, data, local_size, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        for (int step = 1; step < num_procs; step *= 2) {
            for (int i = 0; i < num_procs; i += 2 * step) {
                int left = i * local_size;
                int mid = (i + step) * local_size;
                int right = (i + 2 * step) * local_size - 1;
                if (right >= data_size) {
                    right = data_size - 1;
                }
                merge(data, left, mid - 1, right);
            }
        }
end_time = MPI_Wtime(); // Record ending time
        printf("\nSorted data (Merge Sort): ");
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

