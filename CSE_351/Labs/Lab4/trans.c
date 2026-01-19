/*
 * CSE 351 Lab 4 (Caches and Cache-Friendly Code)
 * Part 2 - Optimizing Matrix Transpose
 *
 * Name(s): Luke Verlangieri
 * NetID(s): lverl23
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[M][N], int B[N][M]);
 * and compute B = A^T.
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1 KiB direct mapped cache with a block size of 32 bytes.
 */

#include <stdio.h>
#include "support/cachelab.h"

int is_transpose(int M, int N, int A[M][N], int B[N][M]);

/*
 * transpose_submit - This is the transpose function that you will be graded
 *     on. Do not change the description string "Transpose submission", as the
 *     driver searches for that string to identify the transpose function to be
 *     graded.
 */

 #define SIZE_64 4 
 #define SIZE_32 8 
 #define INTS_PER_BLOCK 8

char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[M][N], int B[N][M]) {

    /* From Problem Statement 
     * 
     * s = 5 -> S = 32   Bytes 
     * k = 5 -> K = 32   Bytes
     * E = 1 -> C = 2048 Bytes
     */

    int row, column, k, l; // 4 
    int arr[INTS_PER_BLOCK]; // 8

    // On the 64 case, B will conflict itself for sizes > 4 
    if(M/INTS_PER_BLOCK >= INTS_PER_BLOCK){
        // Move in number of blocks through matrix
        for(row = 0; row < M; row += SIZE_64) {
            for(column = 0; column < N; column += SIZE_64) {

                // within size x size matrix, optimize by reading all
                // values in A, then writing, reduces conflict misses.
                for(k = 0; k < SIZE_64; k++) {
                    for(l = 0; l < SIZE_64; l++) {
                        arr[l] = A[row + k][column + l];
                    }

                    for(l = 0; l < SIZE_64; l++) {
                        B[column + l][row + k] = arr[l];
                    }
                }
            }
        }
    } else {
        for(row = 0; row < M; row += SIZE_32) {
            for(column = 0; column < N; column += SIZE_32) {
                for(k = 0; k < SIZE_32; k++) {
                    for(l = 0; l < SIZE_32; l++) {
                        arr[l] = A[row + k][column + l];
                    }

                    for(l = 0; l < SIZE_32; l++) {
                        B[column + l][row + k] = arr[l];
                    }
                }
            }
        }
    }
}


// You can define additional transpose functions below. We've defined a simple
// one below to help you get started.

char transpose_8_8_desc[] = " 8x8 (ints) Transpose submission";
void transpose_8_8(int M, int N, int A[M][N], int B[N][M]) {

    #define NUM_INTS_PER_BLOCK 32 / sizeof(int)

    /* From Problem Statement 
     * 
     * s = 5 -> S = 32   Bytes 
     * k = 5 -> K = 32   Bytes
     * E = 1 -> C = 2048 Bytes
     */

    int i, j, k, l;

    for(i = 0; i < M; i += 8) {
        for(j = 0; j < N; j += 8) {
            
            // Transpose K/2 x K/2 
            for(k = i; k < i + 8; k++) {
                for(l = j; l < j + 8; l++){
                    B[l][k] = A[k][l];
                }
            }

        }
    }
}

char transpose_4_4_desc[] = " 4 x 4 Transpose submission";
void transpose_4_4(int M, int N, int A[M][N], int B[N][M]) {

    #define NUM_INTS_PER_BLOCK 32 / sizeof(int)

    /* From Problem Statement 
     * 
     * s = 5 -> S = 32   Bytes 
     * k = 5 -> K = 32   Bytes
     * E = 1 -> C = 2048 Bytes
     */

    int i, j, k, l;

    for(i = 0; i < M; i += 4) {
        for(j = 0; j < N; j += 4) {
            
            // Transpose K/2 x K/2 
            for(k = i; k < i + 4; k++) {
                for(l = j; l < j + 4; l++){
                    B[l][k] = A[k][l];
                }
            }

        }
    }
}

char transpose_4_4_diag_desc[] = "4 x 4 w/ diag optimization";
void transpose_4_4_diag(int M, int N, int A[M][N], int B[N][M]) {

    #define NUM_INTS_PER_BLOCK 32 / sizeof(int)

    /* From Problem Statement 
     * 
     * s = 5 -> S = 32   Bytes 
     * k = 5 -> K = 32   Bytes
     * E = 1 -> C = 2048 Bytes
     */

    int i, j, k, l;

    for(int z = 0; z < M; z++){
        B[z][z] = A[z][z];
    }

    for(i = 0; i < M; i += NUM_INTS_PER_BLOCK/2) {
        for(j = 0; j < N; j += NUM_INTS_PER_BLOCK/2) {
            
            // Transpose K/2 x K/2 
            for(k = i; k < i + NUM_INTS_PER_BLOCK/2; k++) {
                for(l = j; l < j + NUM_INTS_PER_BLOCK/2; l++){
                    if(l != k){
                        B[l][k] = A[k][l];
                    }
                }
            }

        }
    }
}

/*
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[M][N], int B[N][M]) {
    int i, j, tmp;

    for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++) {
            tmp = A[i][j];
            B[j][i] = tmp;
        }
    }

}


/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions() {
    /* Register your solution function */
    registerTransFunction(transpose_submit, transpose_submit_desc);

    /* Register any additional transpose functions */
    registerTransFunction(trans, trans_desc);
    // registerTransFunction(transpose_8_8, transpose_8_8_desc);
    // registerTransFunction(transpose_4_4, transpose_4_4_desc);
    // registerTransFunction(transpose_4_4_diag, transpose_4_4_diag_desc);


}



/*
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[M][N], int B[N][M]) {
    int i, j;

    for (i = 0; i < M; i++) {
        for (j = 0; j < N; ++j) {
            if (A[i][j] != B[j][i]) {
                return 0;
            }
        }
    }
    return 1;
}
