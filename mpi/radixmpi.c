#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

// Function to get the maximum value in arr[]
long long int getMax(long long int arr[], int n) {
    long long int max = arr[0];
    for (int i = 1; i < n; i++) {
        if (arr[i] > max) {
            max = arr[i];
        }
    }
    return max;
}

// A function to do counting sort of arr[] according to the digit represented by exp
void countSort(long long int arr[], int n, int exp) {
    long long int *output = (long long int *)malloc(n * sizeof(long long int));
    int count[10] = {0};

    for (int i = 0; i < n; i++) {
        count[(arr[i] / exp) % 10]++;
    }

    for (int i = 1; i < 10; i++) {
        count[i] += count[i - 1];
    }

    for (int i = n - 1; i >= 0; i--) {
        output[count[(arr[i] / exp) % 10] - 1] = arr[i];
        count[(arr[i] / exp) % 10]--;
    }

    for (int i = 0; i < n; i++) {
        arr[i] = output[i];
    }

    free(output);
}

// Radix sort function
void radixSort(long long int arr[], int n) {
    long long int max = getMax(arr, n);

    for (int exp = 1; max / exp > 0; exp *= 10) {
        countSort(arr, n, exp);
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n = (100 * 1024 * 1024) / sizeof(long long int);  // 100mb
    long long int *arr = (long long int *)malloc(n * sizeof(long long int));

    srand(rank + 1);
    for (int i = 0; i < n; i++) {
        arr[i] = rand() % 1000;
    }

    double start_time, end_time;

    if (rank == 0) {
        printf("\n\n");
    }

    int local_n = n / size;
    long long int *local_arr = (long long int *)malloc(local_n * sizeof(long long int));

    start_time = MPI_Wtime();
    MPI_Scatter(arr, local_n, MPI_LONG_LONG, local_arr, local_n, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    radixSort(local_arr, local_n);

    MPI_Gather(local_arr, local_n, MPI_LONG_LONG, arr, local_n, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        radixSort(arr, n);
        end_time = MPI_Wtime();
        printf("Sorted array: ");
        for (int i = 0; i < n; i++) {
            printf("%lld ", arr[i]);
        }
        printf("\n");

        printf("Time taken for sorting: %f seconds\n", end_time - start_time);
    }

    free(arr);
    free(local_arr);
    MPI_Finalize();
    return 0;
}

