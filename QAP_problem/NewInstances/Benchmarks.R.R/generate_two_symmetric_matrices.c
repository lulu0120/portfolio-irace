#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 100

static void generate_symmetric_matrix(int matrix[N][N]) {
    for (int i = 0; i < N; i++) {
        matrix[i][i] = 0;  // zero diagonal
        for (int j = i + 1; j < N; j++) {
            int value = (rand() % 99) + 1;  // uniform integer in [1, 99]
            matrix[i][j] = value;
            matrix[j][i] = value;
        }
    }
}

static void print_matrix(FILE *out, int matrix[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(out, "%d", matrix[i][j]);
            if (j < N - 1) {
                fputc(' ', out);
            }
        }
        fputc('\n', out);
    }
}

int main(void) {
    int matrix1[N][N];
    int matrix2[N][N];

    srand((unsigned int)time(NULL));

    generate_symmetric_matrix(matrix1);
    generate_symmetric_matrix(matrix2);

    // First line: matrix dimensions
    printf("%d %d\n", N, N);

    // First matrix, then second matrix
    print_matrix(stdout, matrix1);
    print_matrix(stdout, matrix2);

    return 0;
}
