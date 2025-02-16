# 🚀 Parallel Sorting Algorithms - Performance Analysis

## 📌 Project Overview
This project conducts a **comparative performance analysis** of **10 different sorting algorithms** using three execution models:  
1️⃣ **Serial Execution** – Baseline performance without parallelism.  
2️⃣ **OpenMP Execution** – Multi-threaded execution in a **shared memory system**.  
3️⃣ **MPI Execution** – Parallel execution in a **distributed memory system**.  

By evaluating execution time, scalability, and efficiency, this project explores the **impact of parallel computing** on sorting large datasets.

---

## 🛠️ Technologies Used
- **C, C++** 🖥️ (Core implementation)
- **POSIX Threads (pthread)** 🔀 (Multi-threading)
- **OpenMP** 🏗️ (Shared-memory parallelism)
- **MPI (Message Passing Interface)** 📡 (Distributed-memory parallelism)
- **Ubuntu (GCC)** 🐧 (Development Environment)

---

## 📊 Sorting Algorithms Implemented
| **Algorithm**      | **Type**             | **Parallelization** |
|--------------------|---------------------|---------------------|
| Bubble Sort       | Comparison-based    | Serial, OpenMP, MPI |
| Selection Sort    | Comparison-based    | Serial, OpenMP, MPI |
| Insertion Sort    | Comparison-based    | Serial, OpenMP, MPI |
| Merge Sort       | Divide & Conquer     | Serial, OpenMP, MPI |
| Quick Sort       | Divide & Conquer     | Serial, OpenMP, MPI |
| Bitonic Sort     | Parallel Sort        | Serial, OpenMP, MPI |
| Comb Sort        | Optimized Bubble     | Serial, OpenMP, MPI |
| Counting Sort    | Non-comparison Sort  | Serial, OpenMP, MPI |
| Heap Sort        | Priority Queue Sort  | Serial, OpenMP, MPI |
| Radix Sort       | Non-comparison Sort  | Serial, OpenMP, MPI |

---

## 📥 Installation & Setup
### **1️⃣ Clone the Repository**
```sh
git clone https://github.com/yourusername/Parallel-Sorting-Algorithms.git
cd Parallel-Sorting-Algorithms
```
### **2️⃣ Compile & Run**
🔹 Serial Execution
```sh
gcc bubble_sort.c -o bubble_sort
./bubble_sort
```
🔹 OpenMP Execution
```sh
gcc -fopenmp quick_sort.c -o quick_sort
./quick_sort
```

🔹 MPI Execution
```sh
mpicc -o mpi_merge_sort mpi_merge_sort.c
mpirun -np 4 ./mpi_merge_sort
```
## 🔬 Research Contribution
### **This project aims to:** 
-✅ Identify the most efficient sorting techniques for large-scale datasets.
-✅ Compare shared-memory vs. distributed-memory performance.
-✅ Optimize sorting algorithms using parallel computing strategies.
-✅ Provide open-source implementations for further research.
