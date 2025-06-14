/**
 * @file trans.c
 * @brief Contains various implementations of matrix transpose
 *
 * Each transpose function must have a prototype of the form:
 *   void trans(size_t M, size_t N, double A[N][M], double B[M][N],
 *              double tmp[TMPCOUNT]);
 *
 * All transpose functions take the following arguments:
 *
 *   @param[in]     M    Width of A, height of B
 *   @param[in]     N    Height of A, width of B
 *   @param[in]     A    Source matrix
 *   @param[out]    B    Destination matrix
 *   @param[in,out] tmp  Array that can store temporary double values
 *
 * A transpose function is evaluated by counting the number of hits and misses,
 * using the cache parameters and score computations described in the writeup.
 *
 * Programming restrictions:
 *   - No out-of-bounds references are allowed
 *   - No alterations may be made to the source array A
 *   - Data in tmp can be read or written
 *   - This file cannot contain any local or global doubles or arrays of doubles
 *   - You may not use unions, casting, global variables, or
 *     other tricks to hide array data in other forms of local or global memory.
 *
 * TODO: fill in your name and Andrew ID below.
 * @author Angela Wu <qianruw@andrew.cmu.edu>
 */

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>

#include "cachelab.h"
#define block_size 8
#define Bsize 64
#define Ssize 8

/**
 * @brief Checks if B is the transpose of A.
 *
 * You can call this function inside of an assertion, if you'd like to verify
 * the correctness of a transpose function.
 *
 * @param[in]     M    Width of A, height of B
 * @param[in]     N    Height of A, width of B
 * @param[in]     A    Source matrix
 * @param[out]    B    Destination matrix
 *
 * @return True if B is the transpose of A, and false otherwise.
 */
#ifndef NDEBUG
static bool is_transpose(size_t M, size_t N, double A[N][M], double B[M][N]) {
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < M; ++j) {
            if (A[i][j] != B[j][i]) {
                fprintf(stderr,
                        "Transpose incorrect.  Fails for B[%zd][%zd] = %.3f, "
                        "A[%zd][%zd] = %.3f\n",
                        j, i, B[j][i], i, j, A[i][j]);
                return false;
            }
        }
    }
    return true;
}
#endif

/*
 * You can define additional transpose functions here. We've defined
 * some simple ones below to help you get started, which you should
 * feel free to modify or delete.
 */

/**
 * @brief A simple baseline transpose function, not optimized for the cache.
 *
 * Note the use of asserts (defined in assert.h) that add checking code.
 * These asserts are disabled when measuring cycle counts (i.e. when running
 * the ./test-trans) to avoid affecting performance.
 */
static void trans_basic(size_t M, size_t N, double A[N][M], double B[M][N],
                        double tmp[TMPCOUNT]) {
    assert(M > 0);
    assert(N > 0);

    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < M; j++) {
            B[j][i] = A[i][j];
        }
    }

    assert(is_transpose(M, N, A, B));
}

/**
 * @brief A contrived example to illustrate the use of the temporary array.
 *
 * This function uses the first four elements of tmp as a 2x2 array with
 * row-major ordering.
 */
static void trans_tmp(size_t M, size_t N, double A[N][M], double B[M][N],
                      double tmp[TMPCOUNT]) {
    assert(M > 0);
    assert(N > 0);

    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < M; j++) {
            size_t di = i % 2;
            size_t dj = j % 2;
            tmp[2 * di + dj] = A[i][j];
            B[j][i] = tmp[2 * di + dj];
        }
    }

    assert(is_transpose(M, N, A, B));
}

void transpose_32(size_t M, size_t N, double A[N][M], double B[M][N],
                  double *tmp);
void transpose_1024(size_t M, size_t N, double A[N][M], double B[M][N],
                    double *tmp);
/**
 * @brief The solution transpose function that will be graded.
 *
 * You can call other transpose functions from here as you please.
 * It's OK to choose different functions based on array size, but
 * this function must be correct for all values of M and N.
 */
static void transpose_submit(size_t M, size_t N, double A[N][M], double B[M][N],
                             double tmp[TMPCOUNT]) {
    // if (M == N)
    //     trans_basic(M, N, A, B, tmp);
    // else
    //     trans_tmp(M, N, A, B, tmp);
    assert(M > 0);
    assert(N > 0);
    if (M == 32 && N == 32)
        transpose_32(M, N, A, B, tmp);
    else if (M == 1024 && N == 1024) {
        transpose_1024(M, N, A, B, tmp);
    } else
        trans_tmp(M, N, A, B, tmp);
}

/**
 * @brief Registers all transpose functions with the driver.
 *
 * At runtime, the driver will evaluate each function registered here, and
 * and summarize the performance of each. This is a handy way to experiment
 * with different transpose strategies.
 */
void registerFunctions(void) {
    // Register the solution function. Do not modify this line!
    registerTransFunction(transpose_submit, SUBMIT_DESCRIPTION);

    // Register any additional transpose functions
    registerTransFunction(trans_basic, "Basic transpose");
    registerTransFunction(trans_tmp, "Transpose using the temporary array");
}

/**
 * @brief Optimized transpose for 32x32 matrix using 8x8 blocking.
 *
 * @requires M == 32 && N == 32
 * @ensures For all i, j: B[j][i] == A[i][j]
 *
 * This function uses 8*8 blocking to improve spatial locality and reduce cache
 * misses. For non-diagonal blocks, data is transposed directly. For diagonal
 * blocks, we avoid writing diagonal elements early to reduce conflict misses.
 *
 * @param M Number of columns in A, must be 32
 * @param N Number of rows in A, must be 32
 * @param A Input matrix A[N][M], read-only
 * @param B Output matrix B[M][N], where result is written
 * @param tmp Temporary buffer (not used)
 */
void transpose_32(size_t M, size_t N, double A[N][M], double B[M][N],
                  double *tmp) {
    assert(M == 32 && N == 32);
    // 8*8 blocks A[i+k][j+0-7]
    for (size_t i = 0; i < N; i += block_size) {
        for (size_t j = 0; j < M; j += block_size) {
            if (i != j) {
                for (size_t k = 0; k < block_size; k++) {
                    for (size_t t = 0; t < block_size; t++) {
                        B[j + t][i + k] = A[i + k][j + t];
                    }
                }
            } else {
                for (size_t k = 0; k < block_size; k++) {
                    for (size_t t = 0; t < block_size; t++) {
                        if (k != t) {
                            B[j + t][i + k] = A[i + k][j + t];
                        }
                    }
                    B[j + k][i + k] = A[i + k][j + k];
                }
            }
        }
    }

    assert(is_transpose(M, N, A, B));
}

/**
 * @brief Optimized transpose for 1024x1024 matrix using multi-level blocking.
 *
 * This function uses coarse-grained (Bsize x Bsize) and fine-grained (Ssize x
 * Ssize) blocking to reduce conflict misses in direct-mapped cache. Diagonal
 * blocks are handled separately using `tmp` to avoid writing into B before the
 * entire subblock is processed, preventing self-conflict in the cache.
 *
 * @requires M == 1024 && N == 1024
 * @ensures For all i, j: B[j][i] == A[i][j]
 *
 * @param M Number of columns in A, must be 1024
 * @param N Number of rows in A, must be 1024
 * @param A Input matrix A[N][M], read-only
 * @param B Output matrix B[M][N], where result is written
 * @param tmp Temporary buffer of size at least Ssize, used for diagonal
 * buffering
 */
void transpose_1024(size_t M, size_t N, double A[N][M], double B[M][N],
                    double *tmp) {
    assert(M == 1024 && N == 1024);

    for (size_t i = 0; i < N; i += Bsize) {
        for (size_t j = 0; j < M; j += Bsize) {
            for (size_t ii = i; ii < i + Bsize; ii += Ssize) {
                for (size_t jj = j; jj < j + Bsize; jj += Ssize) {
                    if (ii != jj) {
                        // Non-diagonal blocks: direct transpose
                        for (size_t x = 0; x < Ssize; ++x) {
                            for (size_t y = 0; y < Ssize; ++y) {
                                B[jj + y][ii + x] = A[ii + x][jj + y];
                            }
                        }
                    } else {
                        // Diagonal block: use tmp to buffer diagonal entries
                        for (size_t x = 0; x < Ssize; ++x) {
                            tmp[x] = A[ii + x][ii + x];
                            for (size_t y = 0; y < Ssize; ++y) {
                                if (x != y)
                                    B[ii + y][ii + x] = A[ii + x][ii + y];
                            }
                        }
                        for (size_t x = 0; x < Ssize; ++x) {
                            B[ii + x][ii + x] = tmp[x]; // Restore diagonal
                        }
                    }
                }
            }
        }
    }

    assert(is_transpose(M, N, A, B));
}
