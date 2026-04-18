#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 100
#define INSTANCES 100

static void generate_symmetric_matrix(int matrix[N][N]) {
    for (int i = 0; i < N; i++) {
        matrix[i][i] = 0;
        for (int j = i + 1; j < N; j++) {
            int value = (rand() % 99) + 1;
            matrix[i][j] = value;
            matrix[j][i] = value;
        }
    }
}

static void print_matrix(FILE *out, int matrix[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(out, "%d", matrix[i][j]);
            if (j < N - 1) fprintf(out, " ");
        }
        fprintf(out, "\n");
    }
}

int main(void) {
    srand((unsigned int)time(NULL));

    for (int inst = 1; inst <= INSTANCES; inst++) {
        char filename[64];
        sprintf(filename, "RandomInstance.%d.dat", inst);

        FILE *f = fopen(filename, "w");
        if (!f) {
            perror("Error opening file");
            return 1;
        }

        int matrix1[N][N];
        int matrix2[N][N];

        generate_symmetric_matrix(matrix1);
        generate_symmetric_matrix(matrix2);

        // First line
        fprintf(f, "%d %d\n\n", N, N);

        // First matrix
        print_matrix(f, matrix1);

        // Empty line after first matrix
        fprintf(f, "\n");

        // Second matrix
        print_matrix(f, matrix2);

        fclose(f);
    }

    return 0;
}
