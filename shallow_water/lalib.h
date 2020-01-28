#ifndef LALIB_H
#define LALIB_H

/**
 * Linear algebra routines
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/**
 * Matrix
 */
typedef struct Matrix {
    int M;
    int N;
    double **vals;
} matrix_t;

/**
 * Print a matrix to stdout
 */
void mat_print(matrix_t mat) {
    for (int i = 0; i < mat.M; i++) {
        for (int j = 0; j < mat.N; j++) {
            printf("%.2e ", mat.vals[i][j]);
        }
        printf("\n");
    }
}

/**
 * Create a new M-by-N matrix
 */
matrix_t mat_new(int M, int N) {
    double **vals;
    vals = malloc(M * sizeof(double *));
    for (int i = 0; i < M; i++) {
        vals[i] = malloc(N * sizeof(double));
    }
    matrix_t mat = { M, N, vals };
    return mat;
}

/**
 * Create a matrix filled with zeros
 */
matrix_t mat_zeros(int M, int N) {
    matrix_t mat = mat_new(M, N);
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            mat.vals[i][j] = 0.0;
        }
    }
    return mat;
}

/**
 * Free matrix memory
 */
void mat_free(matrix_t mat) {
    free(mat.vals);
}

/**
 * Find maximum of matrix
 */
double mat_maximum(matrix_t mat) {
    double max = -INFINITY;
    for (int i = 0; i < mat.M; i++) {
        for (int j = 0; j < mat.N; j++) {
            if (mat.vals[i][j] > max) {
                max = mat.vals[i][j];
            }
        }
    }
    return max;
}

/**
 * Find minimum of matrix
 */
double mat_minimum(matrix_t mat) {
    double min = INFINITY;
    for (int i = 0; i < mat.M; i++) {
        for (int j = 0; j < mat.N; j++) {
            if (mat.vals[i][j] < min) {
                min = mat.vals[i][j];
            }
        }
    }
    return min;
}

/**
 * Compute minor matrix
 */
matrix_t la_minor_matrix(matrix_t mat, int row, int col) {

    matrix_t minor = mat_new(mat.M - 1, mat.N - 1);
    int idest = 0;
    for (int isrc = 0; isrc < mat.M; isrc++) {
        if (isrc == row) {
            continue;
        }
        int jdest = 0;
        for (int jsrc = 0; jsrc < mat.N; jsrc++) {
            if (jsrc == col) {
                continue;
            }
            minor.vals[idest][jdest] = mat.vals[isrc][jsrc];
            jdest++;
        }
        idest++;
    }
    return minor;
}

/**
 * Compute matrix determinant
 */
double la_determinant(matrix_t mat) {

    // Base cases
    if (mat.M == 1) {
        return mat.vals[0][0];
    }
    if (mat.M == 2) {
        return mat.vals[0][0] * mat.vals[1][1] - 
            mat.vals[0][1] * mat.vals[1][0];
    }

    // Otherwise, compute recursively
    double sign_fac = 1.0;
    double det = 0.0;
    int i = 0;
    matrix_t minor = mat_new(mat.M - 1, mat.N - 1);
    for (int j = 0; j < mat.N; j++) {
        matrix_t minor = la_minor_matrix(mat, i, j);
        det += sign_fac * mat.vals[i][j] * la_determinant(minor);
        sign_fac = -sign_fac;
        mat_free(minor);
    }
    return det;
}

/**
 * Compute cofactor matrix
 */
matrix_t la_cofactor_matrix(matrix_t mat) {
    
    matrix_t cofactor = mat_new(mat.M, mat.N);
    for (int i = 0; i < mat.M; i++) {
        for (int j = 0; j < mat.N; j++) {
            matrix_t minor = la_minor_matrix(mat, i, j);
            cofactor.vals[i][j] = pow(-1.0, (double) (i + j)) * 
                la_determinant(minor);
            mat_free(minor);
        }
    }
    return cofactor;
}

/** 
 * Compute matrix transpose
 */
matrix_t la_transpose(matrix_t mat) {

    matrix_t transpose = mat_new(mat.M, mat.N);
    for (int i = 0; i < mat.M; i++) {
        for (int j = 0; j < mat.N; j++) {
            transpose.vals[j][i] = mat.vals[i][j];
        }
    }
    return transpose;
}

/**
 * Compute matrix inverse
 */
matrix_t la_inverse(matrix_t mat) {
    
    double det = la_determinant(mat);
    matrix_t aux = la_cofactor_matrix(mat);
    matrix_t inverse = la_transpose(aux);
    for (int i = 0; i < inverse.M; i++) {
        for (int j = 0; j < inverse.N; j++) {
            inverse.vals[i][j] /= det;
        }
    }
    mat_free(aux);
    return inverse;
}

/**
 * Compute inner product of a row and column of two matrices
 */
double la_dot(matrix_t mat1, matrix_t mat2, int row1, int col2) {
    double dot = 0.0;
    for (int i = 0; i < mat1.N; i++) {
        dot += mat1.vals[row1][i] * mat2.vals[i][col2];
    }
    return dot;
}

/**
 * Compute matrix product
 */
matrix_t la_product(matrix_t mat1, matrix_t mat2) {

    matrix_t prod = mat_new(mat1.M, mat1.N);
    for (int i = 0; i < prod.M; i++) {
        for (int j = 0; j < prod.N; j++) {
            prod.vals[i][j] = la_dot(mat1, mat2, i, j);
        }
    }
    return prod;
}

#endif
