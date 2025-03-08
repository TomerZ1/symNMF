#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "symnmf.h"

/* Constants */
#define int_INFINITY 2147483647
#define EPSILON 1e-4
#define MAX_ITER 300
#define BETA 0.5

double** sym(double **data, int N, int d);
double** ddg(int N, double **A);
double** norm (int N, double **D, double **A);
void symnmf(double **curr_H, double **W, int N, int d);

void compute_root (double **D, int N, double *D_root);
void matrix_transpose(double **mat, double **transposed_mat, int N, int k);
void matrix_multiply(double **mat1, double **mat2, double **result, int n, int m, int p);
void matrix_subtract(double **mat1, double **mat2, double **result, int N, int M);
double frobenius_norm(double **mat, int N, int M);

void allocate_matrix_memory(double ***matrix, int N, int k);
void free_matrix_memory(double **matrix, int rows);

void get_matrix_dimensions(const char* filename, int* N, int* d);
double** read_file_to_matrix(const char* filename, int N, int d);
void print_matrix(double** matrix, int N, int d);

int main(int argc, char *argv[]) {
    char *goal, *filename;
    double **points, **A, **D, **W;
    int N = 0, d = 0;

    if (argc != 3) { /* Ensure correct number of arguments and extract */
        printf("An Error Has Occurred\n");
        return 1;
    }
    goal = argv[1]; 
    filename = argv[2];

    /* Determine matrix dimensions and read the file into matrix */
    get_matrix_dimensions(filename, &N, &d);
    points = read_file_to_matrix(filename, N, d);

    if (strcmp(goal, "sym") == 0){ /* Validate goal */
        A = sym(points, N, d); /* Compute similarity matrix */
        print_matrix(A, N, N); /* Print matrix */
    } 
    else if (strcmp(goal, "ddg") == 0){
        A = sym(points, N, d); /* Compute similarity matrix */
        D = ddg(N, A); /* Compute diagonal degree matrix */
        print_matrix(D, N, N); /* Print matrix */
        free_matrix_memory(D, N); /* Free memory */
    } 
    else if (strcmp(goal, "norm") == 0){
        A = sym(points, N, d); /* Compute similarity matrix */
        D = ddg(N, A); /* Compute diagonal degree matrix */
        W = norm(N, D, A); /* Compute norm matrix */
        print_matrix(W, N, N); /* Print matrix */
        free_matrix_memory(D, N); /* Free memory */
        free_matrix_memory(W, N); /* Free memory */
    }
    else { 
        printf("An Error Has Occurred\n");
        return 1;
    }
    free_matrix_memory(A, N); /* Free memory, if we are here we allocated at least A */
    free_matrix_memory(points, N); /* Free memory */
    return 0;
}

double** sym(double **data, int N, int d) {
    /**
    * @brief Computes the similarity matrix for a given set of points.
    *
    * This function calculates the similarity matrix A using the given data points.
    *
    * @param data A 2D array (N x d) containing the input points.
    * @param N The number of points (rows in the matrix).
    * @param d The number of dimensions (columns in the matrix).
    * @return A dynamically allocated 2D array (N x N) containing the similarity matrix.
    */
    int i, j, k;
    double sum, diff;
    double **A;

    /* Allocate memory for similarity matrix */
    allocate_matrix_memory(&A, N, N);

    /* Compute similarity matrix */
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            if (i == j) {
                A[i][j] = 0.0; /* No self-similarity */
            } else {
                sum = 0.0;
                for (k = 0; k < d; k++) {
                    diff = data[i][k] - data[j][k]; /* Compute difference */
                    sum += diff * diff; /* Compute sum of squares */
                }
                A[i][j] = exp(-sum / 2.0 ); /* Apply Gaussian similarity */
            }
        }
    }
    return A;
}

/* Compute diagonal degree matrix D for similarity matrix A */
double** ddg (int N, double **A){
    int i, j;
    double sum = 0.0;
    double **D;

    /* Allocate memory for diagonal degree matrix */
    allocate_matrix_memory(&D, N, N);

    /* Compute diagonal degree matrix */
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            sum += A[i][j];
            D[i][j] = 0; /* Initialize D to zeros */
        }
        D[i][i] = sum; /* Assign sum to diagonal */
        sum = 0;
    }

    return D;
}

/* Compute normalized matrix W for similarity matrix A and ddg matrix D */
double** norm (int N, double **D, double **A){
    int i, j;
    double *D_root;
    double **W;

    /* Allocate memory for norm matrix */
    allocate_matrix_memory(&W, N, N);

    /* Allocate memory for the root matrix */
    D_root = (double*)malloc(N * sizeof(double));

    /* Compute square root of diagonal matrix D */
    compute_root(D, N, D_root);

    /* Compute norm matrix */
    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++){
            W[i][j] = D_root[i] * A[i][j] * D_root[j];
        }
    }

    /* Free memory */
    free(D_root);

    return W;
}

/* Compute square root of diagonal matrix D */
void compute_root (double **D, int N, double *D_root){
    int i;
    for (i = 0; i < N; i++){ /* Compute square root of diagonal matrix D */
        D_root[i] = 1.0 / sqrt(D[i][i]);
    }
/* Store it in 1D array because this is just the diagonal */
}

/* Compute transpose of matrix */
void matrix_transpose(double **mat, double **transposed_mat, int N, int k){
    int i, j; /* transposed_mat is already initialized */
    for (i = 0; i < N; i++){
        for (j = 0; j < k; j++){
            transposed_mat[j][i] = mat[i][j];
        }
    }
}

/* Compute matrix multiplication */
void matrix_multiply(double **mat1, double **mat2, double **result, int n, int m, int p) {
    int i, j, k;
    for (i = 0; i < n; i++) {
        for (j = 0; j < p; j++) {
            result[i][j] = 0.0; 
            for (k = 0; k < m; k++) {
                result[i][j] += mat1[i][k] * mat2[k][j];  /* Compute dot product */
            }
        }
    }
}

/* Compute matrix subtraction */
void matrix_subtract(double **mat1, double **mat2, double **result, int N, int M){
    int i, j;
    for (i = 0; i < N; i++){
        for (j = 0; j < M; j++){
            result[i][j] = mat1[i][j] - mat2[i][j]; /* Compute difference */
        }
    }
}

/* Compute squared Frobenius norm of matrix */
double frobenius_norm(double **mat, int N, int M){
    int i, j;
    double sum = 0.0;
    for (i = 0; i < N; i++){ 
        for (j = 0; j < M; j++){
            sum += mat[i][j] * mat[i][j]; /* Compute sum of squares */
        }
    }
    return sum; /* Return squared norm */
}

/* Free memory allocated for matrix */
void free_matrix_memory(double **matrix, int rows) {
    int i;
    for (i = 0; i < rows; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

/* Allocate memory for matrix */
void allocate_matrix_memory(double ***matrix, int N, int k) {
    int i;
    *matrix = (double**)malloc(N * sizeof(double*));
    if (*matrix == NULL) {
        printf("Memory allocation failed\n");
        exit(1);
    }
    
    for (i = 0; i < N; i++) {
        (*matrix)[i] = (double*)malloc(k * sizeof(double));
        if ((*matrix)[i] == NULL) {
            printf("Memory allocation failed\n");
            exit(1);
        }
    }
}

/* Compute symmetric NMF */
void symnmf(double **curr_H, double **W, int N, int k){
    int i, j, iter;
    double **new_H, **WH, **HtH, **denomirator, **Ht, **margin_matrix; /* WH=W*H, HtH= H^t * H , Ht = H^t  */

    /* Allocate memory for matrices */
    allocate_matrix_memory(&new_H, N, k);
    allocate_matrix_memory(&WH, N, k);
    allocate_matrix_memory(&HtH, k, k);
    allocate_matrix_memory(&denomirator, N, k);
    allocate_matrix_memory(&Ht, k, N);
    allocate_matrix_memory(&margin_matrix, N, k);
    for (iter = 0; iter < MAX_ITER; iter++){  /* Update H step by step */
        matrix_transpose(curr_H, Ht, N, k);
        matrix_multiply(W, curr_H, WH, N, N, k); /* Nomirator */
        matrix_multiply(Ht, curr_H, HtH, k, N, k);
        matrix_multiply(curr_H, HtH, denomirator, N, k, k); /* Denomirator */
        for (i = 0; i < N; i++){ /* create new H */
            for (j = 0; j < k; j++){
                new_H[i][j] = curr_H[i][j] * (1 - BETA + BETA * (WH[i][j] / denomirator[i][j]));
            }
        }
        matrix_subtract(new_H, curr_H, margin_matrix, N, k); 
        for (i = 0; i < N; i++){  /* Update H */
            for (j = 0; j < k; j++){
                curr_H[i][j] = new_H[i][j];
            }
        }
        /* Check if the margin matrix is less than epsilon- convergence */
        if (frobenius_norm(margin_matrix, N, k) < EPSILON){
            break;
        }
    }
    /* Free memory */
    free_matrix_memory(new_H, N);
    free_matrix_memory(WH, N);
    free_matrix_memory(HtH, k);
    free_matrix_memory(denomirator, N);
    free_matrix_memory(Ht, k);
    free_matrix_memory(margin_matrix, N);
}

/* Function to determine N (rows) and d (columns) from file */
void get_matrix_dimensions(const char* filename, int* N, int* d) {
    FILE* file = fopen(filename, "r");
    char line[1024];
    char* token;

    if (!file) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    
    if (fgets(line, sizeof(line), file) == NULL) {
        printf("An Error Has Occurred\n");
        fclose(file);
        exit(1);
    }
    
    token = strtok(line, ",");
    while (token) {
        (*d)++;
        token = strtok(NULL, ",");
    }
    (*N) = 1;
    
    while (fgets(line, sizeof(line), file)) {
        (*N)++;
    }
    
    fclose(file);
}

/* Function to read file into a 2D matrix */
double** read_file_to_matrix(const char* filename, int N, int d) {
    FILE* file = fopen(filename, "r");
    char line[1024];
    double** matrix;
    int row = 0, col;
    char* token; 
    /* Allocate memory for matrix */
    if (!file) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    matrix = (double**)malloc(N * sizeof(double*));
    if (!matrix) {
        printf("An Error Has Occurred\n");
        fclose(file);
        exit(1);
    }
    for (row = 0; row < N; row++) {
        matrix[row] = (double*)malloc(d * sizeof(double));
        if (!matrix[row]) {
            printf("An Error Has Occurred\n");
            fclose(file);
            exit(1);
        }
    }
    /* Read file into matrix */
    row = 0;
    while (fgets(line, sizeof(line), file)) {
        col = 0;
        token = strtok(line, ",");
        while (token) {
            matrix[row][col] = atof(token);
            token = strtok(NULL, ",");
            col++;
        }
        row++;
    }
    fclose(file);
    return matrix;
}

/* Function to print a 2D matrix */
void print_matrix(double** matrix, int N, int d) {
    int i, j;
    printf("Matrix (%d x %d):\n", N, d);
    for (i = 0; i < N; i++) {
        for (j = 0; j < d; j++) {
            printf("%.4f", matrix[i][j]);
            if (j < d - 1) {
                printf(",");
            }
        }
        printf("\n");
    }
}