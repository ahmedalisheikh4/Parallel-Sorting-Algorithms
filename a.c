#include<stdio.h>
#include<omp.h>
#include<stdlib.h>
#include<time.h>
void bubbleSort(long long int arr[], int n) { //serial
    for (int i = 0; i < n - 1; i++) {
    
        for (int j = 0; j < n - i - 1; j++) {
        
            if (arr[j] > arr[j + 1]) {
                // Swap arr[j] and arr[j + 1]
                long long int temp = arr[j];
                arr[j] = arr[j + 1];
                arr[j + 1] = temp;
            }
        }
    }
}

void parallelBubbleSort(long long int arr[], int n) {
    int swapped;
    #pragma omp task
    for (int i = 0; i < n - 1; i++) {
        swapped = 0;
        #pragma omp parallel for shared(arr,n) reduction(||:swapped)
        for (int j = 0; j < n - i - 1; j++) {
            if (arr[j] > arr[j + 1]) {
                // Swap arr[j] and arr[j + 1]
                long long int temp = arr[j];
                arr[j] = arr[j + 1];
                arr[j + 1] = temp;
                swapped = 1;
            }
        }
        // If no two elements were swapped in the inner loop, the array is sorted
        if (!swapped) {
            break;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void combSort(long long int arr[], int n) {
    int gap = n;
    double shrink = 1.3;
    int swapped = 1;

    while (gap > 1 || swapped) {
        if (gap > 1) {
            // Shrink the gap
            gap = (int)(gap / shrink);
            if (gap < 1)
                gap = 1;
        }

        swapped = 0;

        // Compare elements with a gap
        for (int i = 0; i < n - gap; i++) {
            int j = i + gap;
            if (arr[i] > arr[j]) {
                // Swap elements
                long long int temp = arr[i];
                arr[i] = arr[j];
                arr[j] = temp;
                swapped = 1;
            }
        }
    }
}

void parallelCombSort(long long int arr[], int n) {
    int gap = n;
    double shrink = 1.3;
    int swapped = 1;

    while (gap > 1 || swapped) {
        if (gap > 1) {
            // Shrink the gap
            gap = (int)(gap / shrink);
            if (gap < 1)
                gap = 1;
        }

        swapped = 0;

        #pragma omp parallel for shared(arr, gap, n) reduction(|| : swapped)
        // Compare elements with a gap
        for (int i = 0; i < n - gap; i++) {
            int j = i + gap;
            if (arr[i] > arr[j]) {
                // Swap elements
                long long int temp = arr[i];
                arr[i] = arr[j];
                arr[j] = temp;
                swapped = 1;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void countingSort(long long int arr[], int n, int range) {
    int output[n];
    int count[range];

    for (int i = 0; i < range; i++) {
        count[i] = 0;
    }

    for (int i = 0; i < n; i++) {
        count[arr[i]]++;
    }

    int output_index = 0;
    for (int i = 0; i < range; i++) {
        while (count[i] > 0) {
            output[output_index] = i;
            output_index++;
            count[i]--;
        }
    }

    for (int i = 0; i < n; i++) {
        arr[i] = output[i];
    }
}

void parallelCountingSort(long long int arr[], int n, int range) {
    int *output = (int *)malloc(n * sizeof(int)); // Allocate a private output array for each thread
    int *count = (int *)malloc(range * sizeof(int));

    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < range; i++) {
            count[i] = 0;
        }

        #pragma omp for
        for (int i = 0; i < n; i++) {
            #pragma omp atomic
            count[arr[i]]++;
        }

        #pragma omp single
        {
            int output_index = 0;
            for (int i = 0; i < range; i++) {
                while (count[i] > 0) {
                    output[output_index] = i;
                    output_index++;
                    count[i]--;
                }
            }
        }

        #pragma omp for
        for (int i = 0; i < n; i++) {
            arr[i] = output[i];
        }
    }

    free(output);
    free(count);
}
 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void heapify(long long int arr[], int n, int i) {
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;

    if (left < n && arr[left] > arr[largest]) {
        largest = left;
    }

    if (right < n && arr[right] > arr[largest]) {
        largest = right;
    }

    if (largest != i) {
        int temp = arr[i];
        arr[i] = arr[largest];
        arr[largest] = temp;
        heapify(arr, n, largest);
    }
}

// Function to perform Heap Sort
void heapSort(long long int arr[], int n) {
    // Build a max heap
    for (int i = n / 2 - 1; i >= 0; i--) {
    	#pragma omp wait
    	{
        heapify(arr, n, i);
        }
    }

    // Extract elements one by one from the heap
    for (int i = n - 1; i > 0; i--) {
    #pragma omp wait
        int temp = arr[0];
        arr[0] = arr[i];
        arr[i] = temp;
        heapify(arr, i, 0);
    }
}

// Parallelized heapify function
void parallelHeapify(long long int arr[], int n, int i) {
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;

    if (left < n && arr[left] > arr[largest]) {
        largest = left;
    }

    if (right < n && arr[right] > arr[largest]) {
        largest = right;
    }

    if (largest != i) {
        int temp = arr[i];
        arr[i] = arr[largest];
        arr[largest] = temp;
        heapify(arr, n, largest);
    }
}

// Parallelized Heap Sort
void parallelHeapSort(long long int arr[], int n) {
    #pragma omp parallel
    {
        #pragma omp single
        for (int i = n / 2 - 1; i >= 0; i--) {
            parallelHeapify(arr, n, i);
        }
        
        #pragma omp single
        for (int i = n - 1; i > 0; i--) {
            int temp = arr[0];
            arr[0] = arr[i];
            arr[i] = temp;
            parallelHeapify(arr, i, 0);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void insertionSort(long long int arr[], int n) {
    for (int i = 1; i < n; i++) {
        long long int key = arr[i];
        int j = i - 1;

        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j--;
        }
        arr[j + 1] = key;
    }
}

void parallelInsertionSort(long long int arr[], int n) {
    int i, j;

    #pragma omp parallel for shared(arr, n) private(i, j) num_threads(8)
    for (i = 1; i < n; i++) {
        long long int key = arr[i];

        // Variable to signal cancellation
        int cancel = 0;

        #pragma omp cancellation point for
        for (j = i - 1; j >= 0; j--) {
            //#pragma omp cancel for
            if (arr[j] > key) {
                arr[j + 1] = arr[j];
            } else {
                // Set the cancel flag to break out of the loop
                #pragma omp atomic write
                cancel = 1;
                break;
            }
        }

        // Check if cancellation was requested
        #pragma omp cancellation point for
        if (cancel) {
        
        }

        arr[j + 1] = key;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////
void merge(long long int arr[], int left, int middle, int right) {
    int n1 = middle - left + 1;
    int n2 = right - middle;

    // Create temporary arrays
    int L[n1], R[n2];

    // Copy data to temp arrays L[] and R[]
    for (int i = 0; i < n1; i++) {
        L[i] = arr[left + i];
    }
    for (int i = 0; i < n2; i++) {
        R[i] = arr[middle + 1 + i];
    }

    // Merge the temporary arrays back into arr[l..r]
    int i = 0, j = 0, k = left;

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

// Recursive merge sort
void mergeSort(long long int arr[], int left, int right) {
    if (left < right) {
        int middle = left + (right - left) / 2;

        // Sort first and second halves
        mergeSort(arr, left, middle);
        mergeSort(arr, middle + 1, right);

        // Merge the sorted halves
        merge(arr, left, middle, right);
    }
}

void parallelMerge(long long int arr[], int left, int middle, int right) {
    int n1 = middle - left + 1;
    int n2 = right - middle;

    long long int L[n1] ;
    long long int R[n2] ;

    for (int i = 0; i < n1; i++) {
        L[i] = arr[left + i];
    }

    for (int i = 0; i < n2; i++) {
        R[i] = arr[middle + 1 + i];
    }

    int i = 0, j = 0, k = left;

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

void parallelMergeSort(long long int arr[], int left, int right) {
    if (left < right) {
        int middle = left + (right - left) / 2;

        #pragma omp task
            parallelMergeSort(arr, left, middle);

        #pragma omp task
            parallelMergeSort(arr, middle + 1, right);
        #pragma omp taskwait
        parallelMerge(arr, left, middle, right);
        }
    }

///////////////////////////////////////////////////////////////////////////////////////////////////////////
void swap(long long int* a,long long  int* b)
{
	long long int t = *a;
	*a = *b;
	*b = t;
}
int partition (long long int arr[], int low, int high)
{
	int pivot = arr[high]; // pivot
	int i = (low - 1); // Index of smaller element
	for (int j = low; j <= high- 1; j++)
	{
		// If current element is smaller than or
		// equal to pivot
		if (arr[j] <= pivot)
		{
			i++; // increment index of smaller element
			swap(&arr[i], &arr[j]);
		}
	}
	swap(&arr[i + 1], &arr[high]);
	return (i + 1);
}

void quickSort(long long int arr[], int low, int high)
{
	if (low < high)
	{
		/* pi is partitioning index, arr[p] is now
		at right place */
		int pi = partition(arr, low, high);
			quickSort(arr,low, pi - 1);
			quickSort(arr, pi + 1, high);
	}
}

int parallelpartition (long long int arr[], int low, int high)
{
	int pivot = arr[high]; // pivot
	int i = (low - 1); // Index of smaller element

	for (int j = low; j <= high- 1; j++)
	{
		if (arr[j] <= pivot)
		{
			i++; 
			swap(&arr[i], &arr[j]);
		}
	}
	swap(&arr[i + 1], &arr[high]);
	return (i + 1);

}

void parallelquickSort(long long int arr[], int left, int right) {
    if (left < right) {
        int pivot = arr[left];
        int i = left, j = right;

        while (i < j) {
            while (arr[j] > pivot)
                j--;
            while (i < j && arr[i] <= pivot)
                i++;
            if (i < j) {
                int temp = arr[i];
                arr[i] = arr[j];
                arr[j] = temp;
            }
        }

        arr[left] = arr[i];
        arr[i] = pivot;

        #pragma omp tasknowait
        {
        	parallelquickSort(arr, left, i - 1);
        }
        
        #pragma omp tasknowait
        {
        	parallelquickSort(arr, i + 1, right);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void selectionSort(long long int arr[], int n) {
    int i, j, min_index;
    for (i = 0; i < n - 1; i++) {
        min_index = i;
        for (j = i + 1; j < n; j++) {
            if (arr[j] < arr[min_index]) {
                min_index = j;
            }
        }
        // Swap the found minimum element with the element at index i
        long long int temp = arr[i];
        arr[i] = arr[min_index];
        arr[min_index] = temp;
    }
}

void parallelSelectionSort(long long int arr[], int n) {
    int i, j, min_index;
    #pragma omp parallel for private(i, j, min_index) shared(arr, n)
    for (i = 0; i < n - 1; i++) {
        min_index = i;
        for (j = i + 1; j < n; j++) {
            if (arr[j] < arr[min_index]) {
                min_index = j;
            }
        }
        // Swap the found minimum element with the element at index i
        long long int temp = arr[i];
        arr[i] = arr[min_index];
        arr[min_index] = temp;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
long long int getMax(long long int arr[], int n) {
    long long int max = arr[0];
    for (int i = 1; i < n; i++) {
        if (arr[i] > max) {
            max = arr[i];
        }
    }
    return max;
}

// Serial counting sort
void serialCountingSortradix(long long int arr[], int n, int exp, long long int output[]) {
    const int range = 10; // Radix is 10
    int count[range];

    // Initialize count array
    for (int i = 0; i < range; i++) {
        count[i] = 0;
    }

    // Count occurrences of each digit
    for (int i = 0; i < n; i++) {
        count[(arr[i] / exp) % range]++;
    }

    // Adjust count array to represent the position of elements
    for (int i = 1; i < range; i++) {
        count[i] += count[i - 1];
    }

    // Place elements in their correct positions
    for (int i = n - 1; i >= 0; i--) {
        output[count[(arr[i] / exp) % range] - 1] = arr[i];
        count[(arr[i] / exp) % range]--;
    }

    // Copy the sorted elements back to the original array
    for (int i = 0; i < n; i++) {
        arr[i] = output[i];
    }
}

// Parallel counting sort for a specific digit (exp)
void parallelCountingSortradix(long long int arr[], int n, int exp, long long int output[]) {
    const int range = 10;
    int count[10] = {0};  // Initialize count array to zeros

    #pragma omp parallel for num_threads(2)
    for (int i = 0; i < n; i++) {
        #pragma omp atomic
        count[(arr[i] / exp) % range]++;
    }

    for (int i = 1; i < range; i++) {
        count[i] += count[i - 1];
    }

    
    for (int i = n - 1; i >= 0; i--) {
        output[count[(arr[i] / exp) % range] - 1] = arr[i];
        #pragma omp atomic
        count[(arr[i] / exp) % range]--;
    }

    #pragma omp parallel for num_threads(2)
    for (int i = 0; i < n; i++) {
        arr[i] = output[i];
    }
}


// Serial Radix Sort
void serialRadixSort(long long int arr[], int n) {
    long long int max = getMax(arr, n);

    for (int exp = 1; max / exp > 0; exp *= 10) {
        long long int* output = (long long int*)malloc(n * sizeof(long long int));
        serialCountingSortradix(arr, n, exp, output);
        free(output);
    }
}
// Parallel Radix Sort
void parallelRadixSort(long long int arr[], int n) {
    long long int max = getMax(arr, n);
    long long int* output = (long long int*)malloc(n * sizeof(long long int));

    for (int exp = 1; max / exp > 0; exp *= 10) {
        parallelCountingSortradix(arr, n, exp, output);
    }

    free(output);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void bitonicMerge(long long int arr[], int low, int count, int dir) {
    if (count > 1) {
        int k = count / 2;
        for (int i = low; i < low + k; i++) {
            if ((arr[i] > arr[i + k]) == dir) {
                // Swap arr[i] and arr[i+k]
                long long int temp = arr[i];
                arr[i] = arr[i + k];
                arr[i + k] = temp;
            }
        }
        bitonicMerge(arr, low, k, dir);
        bitonicMerge(arr, low + k, k, dir);
    }
}

void bitonicSort(long long int arr[], int low, int count, int dir) {
    if (count > 1) {
        int k = count / 2;

        // Sort in ascending order
        bitonicSort(arr, low, k, 1);

        // Sort in descending order
        bitonicSort(arr, low + k, k, 0);

        // Merge the entire sequence
        bitonicMerge(arr, low, count, dir);
    }
}



void bitonicSortParallel(long long int arr[], int low, int count, int dir) {
    if (count > 1) {
        int k = count / 2;

        #pragma omp task
        bitonicSortParallel(arr, low, k, 1);

        #pragma omp task
        bitonicSortParallel(arr, low + k, k, 0);

        #pragma omp tasknowait

        bitonicMerge(arr, low, count, dir);
    }
}

void convert(long long int *original,long long int *copy,int SIZE){
	for(int i =0 ; i < SIZE; i++){
		copy[i] = original[i];
	}	
}

int main()
{
double t1_bubbles,t2_bubbles,t1_bubblep,t2_bubblep,
t1_combs,t2_combs,t1_combp,t2_combp,
t1_counts,t2_counts,t1_countp,t2_countp,
t1_heaps,t2_heaps,t1_heapp,t2_heapp,
t1_insertions,t2_insertions,t1_insertionp,t2_insertionp,
t1_merges,t2_merges,t1_mergep,t2_mergep,
t1_quicks,t2_quicks,t1_quickp,t2_quickp,
t1_selections,t2_selections,t1_selectionp,t2_selectionp,
t1_radixs,t2_radixs,t1_radixp,t2_radixp,
t1_bitonics,t2_bitonics,t1_bitonicp,t2_bitonicp;
    int n;
    //omp_set_nested(1);
    int threads;
    printf("Enter Number of threads:" );
    scanf("%d",&threads);
    omp_set_num_threads(threads);
    omp_set_dynamic(0);
    printf("Enter the memory size in MB: ");
    scanf("%d", &n);
    long long int size = (long long int)n * 1024 * 1024 / sizeof(long long int);
    long long int SIZE = size;
 
    FILE *file = fopen("random_numbers.txt", "w");
    if (file == NULL) {
        printf("Error opening file.\n");
        return -1;
    }

    srand(time(NULL)); // Seed the random number generator

    for (int i = 0; i < SIZE; ++i) {
        fprintf(file, "%lld\n", (long long)(rand()%10000));
    }

    fclose(file);

    // Reading integers from the file into a long long int array
    long long int *initial = (long long int*)(malloc(sizeof (long long int) * SIZE ));

    file = fopen("random_numbers.txt", "r");
    if (file == NULL) {
        printf("Error opening file.\n");
        return -1;
    }

    for (int i = 0; i < SIZE; ++i) {
        if (fscanf(file, "%lld", &initial[i]) != 1) {
            printf("Error reading from file.\n");
            return -1;
        }
    }
    fclose(file);

    // Printing the first few elements of the array as a check
    //printf("First few elements of the array:\n");
    //for (int i = 0; i < 10; ++i) {
      //  printf("%lld\n", numbers[i]);
    //}
    
    //store original array of numbers
    long long int *original = (long long int*)(malloc(sizeof (long long int) * SIZE ));;
	for(int i =0 ; i < SIZE; i++){
		original[i] = initial[i];
	}
	long long int *numbers = (long long int*)(malloc(sizeof (long long int) * SIZE ));;
	convert(original,numbers,SIZE);
/*
 t1_bubbles = omp_get_wtime();                    	
bubbleSort(numbers,SIZE);
 t2_bubbles = omp_get_wtime();

if(threads == 1){
	convert(original,numbers,SIZE);
	 t1_bubblep = omp_get_wtime();                     
	bubbleSort(numbers,SIZE);
	 t2_bubblep = omp_get_wtime();
}
else{
convert(original,numbers,SIZE);
	 t1_bubblep = omp_get_wtime();                     
	parallelBubbleSort(numbers,SIZE);
	 t2_bubblep = omp_get_wtime();
}
*/
convert(original,numbers,SIZE);
 t1_combs = omp_get_wtime();
combSort(numbers,SIZE);
 t2_combs = omp_get_wtime();

if(threads == 1){
convert(original,numbers,SIZE);
 t1_combp = omp_get_wtime();
combSort(numbers,SIZE);
 t2_combp = omp_get_wtime();
}
else{
convert(original,numbers,SIZE);
 t1_combp = omp_get_wtime();
parallelCombSort(numbers, SIZE);
 t2_combp = omp_get_wtime();
}

convert(original,numbers,SIZE);
 t1_heaps = omp_get_wtime();
heapSort(numbers,SIZE);
 t2_heaps = omp_get_wtime();

if(threads > 1){
convert(original,numbers,SIZE);
 t1_heapp = omp_get_wtime();
parallelHeapSort(numbers, SIZE);
 t2_heapp = omp_get_wtime();
}

else{
 convert(original,numbers,SIZE);
 t1_heapp = omp_get_wtime();
heapSort(numbers, SIZE);
 t2_heapp = omp_get_wtime();
}
/*
convert(original,numbers,SIZE);
 t1_insertions = omp_get_wtime();
insertionSort(numbers, SIZE);
 t2_insertions = omp_get_wtime();

if(threads > 1){
convert(original,numbers,SIZE);
 t1_insertionp = omp_get_wtime();
parallelInsertionSort(numbers, SIZE);
 t2_insertionp = omp_get_wtime();
}
else{
convert(original,numbers,SIZE);
 t1_insertionp = omp_get_wtime();
insertionSort(numbers, SIZE);
 t2_insertionp = omp_get_wtime();
}
*/
convert(original,numbers,SIZE);
 t1_merges = omp_get_wtime();
mergeSort(numbers, 0, SIZE - 1);
 t2_merges = omp_get_wtime();

if(threads > 1){
convert(original,numbers,SIZE);
 t1_mergep = omp_get_wtime();
parallelMergeSort(numbers, 0, SIZE - 1);
 t2_mergep = omp_get_wtime();
}
else{
convert(original,numbers,SIZE);
 t1_mergep  = omp_get_wtime();
mergeSort(numbers, 0, SIZE - 1);
 t2_mergep = omp_get_wtime();
}
/*
convert(original,numbers,SIZE);
 t1_selections = omp_get_wtime();
selectionSort(numbers, SIZE);
 t2_selections = omp_get_wtime();

if(threads > 1){
convert(original,numbers,SIZE);
 t1_selectionp = omp_get_wtime();
parallelSelectionSort(numbers, SIZE);
 t2_selectionp = omp_get_wtime();
}
else{
convert(original,numbers,SIZE);
 t1_selectionp = omp_get_wtime();
selectionSort(numbers, SIZE);
 t2_selectionp = omp_get_wtime();
}
*/
convert(original,numbers,SIZE);
 t1_radixs = omp_get_wtime();
serialRadixSort(numbers, SIZE);
 t2_radixs = omp_get_wtime();

if(threads > 1){
convert(original,numbers,SIZE);
 t1_radixp = omp_get_wtime();
parallelRadixSort(numbers, SIZE);
 t2_radixp = omp_get_wtime();
}
else{
convert(original,numbers,SIZE);
 t1_radixp = omp_get_wtime();
serialRadixSort(numbers, SIZE);
 t2_radixp = omp_get_wtime();
}

convert(original,numbers,SIZE);
 t1_bitonics = omp_get_wtime();
bitonicSort(numbers, 0, SIZE,1);
 t2_bitonics = omp_get_wtime();

if(threads > 1){
convert(original,numbers,SIZE);
 t1_bitonicp = omp_get_wtime();
bitonicSortParallel(numbers, 0, SIZE,1);
 t2_bitonicp = omp_get_wtime();
}
else{
convert(original,numbers,SIZE);
 t1_bitonicp = omp_get_wtime();
bitonicSort(numbers, 0, SIZE,1);
 t2_bitonicp = omp_get_wtime();
}

convert(original,numbers,SIZE);
 t1_quicks = omp_get_wtime();
quickSort(numbers, 0, SIZE-1);
 t2_quicks = omp_get_wtime();

if(threads > 1){
convert(original,numbers,SIZE);
 t1_quickp = omp_get_wtime();
parallelquickSort(numbers, 0, SIZE-1);
 t2_quickp = omp_get_wtime();
}
else{
convert(original,numbers,SIZE);
 t1_quickp = omp_get_wtime();
quickSort(numbers, 0, SIZE-1);
 t2_quickp = omp_get_wtime();
}

convert(original,numbers,SIZE);
 t1_counts = omp_get_wtime();
countingSort(numbers,SIZE, 10000);
 t2_counts = omp_get_wtime();

if(threads > 1){
convert(original,numbers,SIZE);
 t1_countp = omp_get_wtime();
parallelCountingSort(numbers,SIZE, 100000);
 t2_countp = omp_get_wtime();
}
else{
 convert(original,numbers,SIZE);
 t1_countp = omp_get_wtime();
countingSort(numbers,SIZE, 100000);
 t2_countp = omp_get_wtime();

}

for(int i =0  ; i< SIZE; i++){
	printf("%lld\n",numbers[i]);
}

// Calculate speedup
double speedupBubble = (t2_bubbles - t1_bubbles)/(t2_bubblep - t1_bubblep);
double speedupComb = (t2_combs - t1_combs)/(t2_combp - t1_combp);
double speedupCounting = (t2_counts - t1_counts)/(t2_countp - t1_countp);
double speedupHeap = (t2_heaps - t1_heaps)/(t2_heapp - t1_heapp);
double speedupInsertion = (t2_insertions - t1_insertions)/(t2_insertionp - t1_insertionp);
double speedupMerge = (t2_merges - t1_merges)/(t2_mergep - t1_mergep);
double speedupQuick = (t2_quicks - t1_quicks)/(t2_quickp - t1_quickp);
double speedupSelection = (t2_selections - t1_selections)/(t2_selectionp - t1_selectionp);
double speedupRadix = (t2_radixs - t1_radixs)/(t2_radixp - t1_radixp);
double speedupBitonic = (t2_bitonics - t1_bitonics)/(t2_bitonicp - t1_bitonicp);

printf("Sorting Algorithms\tSerial\t\tParallel\tSpeedUp\n");
printf("Bubble Sort\t\t%fs\t%fs\t%fs\n",t2_bubbles - t1_bubbles,t2_bubblep - t1_bubblep,speedupBubble);
printf("Comb Sort\t\t%fs\t%fs\t%fs\n",t2_combs - t1_combs,t2_combp - t1_combp,speedupComb);
printf("Count Sort\t\t%fs\t%fs\t%fs\n",t2_counts - t1_counts,t2_countp - t1_countp,speedupCounting);
printf("Heap Sort\t\t%fs\t%fs\t%fs\n",t2_heaps - t1_heaps,t2_heapp - t1_heapp,speedupHeap);
printf("Insertion Sort\t\t%fs\t%fs\t%fs\n",t2_insertions - t1_insertions,t2_insertionp - t1_insertionp,speedupInsertion);
printf("Merge Sort\t\t%fs\t%fs\t%fs\n",t2_merges - t1_merges,t2_mergep - t1_mergep,speedupMerge);
printf("Quick Sort\t\t%fs\t%fs\t%fs\n",t2_quicks - t1_quicks,t2_quickp - t1_quickp,speedupQuick);
printf("Selection Sort\t\t%fs\t%fs\t%fs\n",t2_selections - t1_selections,t2_selectionp - t1_selectionp,speedupSelection);
printf("Radix Sort\t\t%fs\t%fs\t%fs\n",t2_radixs - t1_radixs,t2_radixp - t1_radixp,speedupRadix);
printf("Bitonic Sort\t\t%fs\t%fs\t%fs\n",t2_bitonics - t1_bitonics,t2_bitonicp - t1_bitonicp,speedupBitonic);




free(initial);
free(original);
free(numbers);

    return 0;
}

